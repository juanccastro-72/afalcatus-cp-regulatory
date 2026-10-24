#!/usr/bin/env python3
"""Reproducible, exploratory regulatory-sequence inventory for annotated plastids.

This is a candidate-discovery workflow, NOT an experimentally validated promoter,
terminator, operon, or regulatory-RNA annotator. See METHODS.md and workbook README.
Internal coordinates are zero-based; exported sequence coordinates are one-based,
inclusive. Every hit includes an oriented genomic position path for exact recovery.
Mixed-strand trans-spliced products are represented by independent exon contexts;
they are not represented by their misleading genome-wide min/max bounding box.

Example (offline after downloading the supplied GenBank record):
python scripts/chloroplast_regulatory_inventory_improved.py \
    --genbank inputs/QB176266.1.gb --outdir . --null-reps 1000 --fold-null-reps 200
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import random
import re
import sys

import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import RNA

VERSION = '2.0.0'
COMP = str.maketrans('ACGTN', 'TGCAN')
SD_SEQUENCES = ('AGGAGG', 'GGAGG', 'AGGAG', 'GAGG', 'AGGA', 'GGAG')
PRODUCT_TYPES = {'CDS', 'rRNA', 'tRNA', 'tmRNA', 'misc_RNA'}


def rc(s):
    return s.translate(COMP)[::-1]


def near_pattern(consensus):
    """Overlapping exact/one-substitution matches (not insertion/deletion)."""
    variants = [consensus]
    for i, base in enumerate(consensus):
        variants.append(consensus[:i] + '[ACGT]' + consensus[i+1:])
    return re.compile('(?=(' + '|'.join(variants) + '))')


P35, P10 = near_pattern('TTGACA'), near_pattern('TATAAT')
OCT = re.compile('(?=(TATAATAT))')
HEX = re.compile('(?=(TATAAT))')
YRTA = re.compile('(?=([CT][AG]TA))')
GAA = re.compile('(?=(GAA))')


def near_hits(seq, pat, consensus):
    return [(m.start(), sum(a != b for a, b in zip(m.group(1), consensus)))
            for m in pat.finditer(seq) if set(m.group(1)) <= set('ACGT')]


def promoter_pairs(seq):
    """Six-base boxes; spacer = bases BETWEEN -35 and -10 boxes, 15..19."""
    h35 = near_hits(seq, P35, 'TTGACA')
    h10 = dict(near_hits(seq, P10, 'TATAAT'))
    hits = []
    for i, mm35 in h35:
        for spacer in range(15, 20):
            j = i + 6 + spacer
            if j in h10:
                mm10 = h10[j]
                score = 12 - mm35 - mm10 - 0.25 * abs(spacer - 17)
                hits.append((i, j + 6, j, spacer, mm35, mm10, score))
    return hits


def sd_hits(seq):
    """Spacing is the number of nucleotides AFTER motif and BEFORE start codon."""
    hits = []
    for motif in SD_SEQUENCES:
        for m in re.finditer('(?=(' + motif + '))', seq):
            a, b = m.start(), m.start() + len(motif)
            spacer = len(seq) - b
            if 3 <= spacer <= 15:
                hits.append((a, b, motif, spacer, float(len(motif))))
    # Do not count the shorter submotifs of a longer SD-like match again.
    return sorted([h for h in hits if not any(
        k[0] <= h[0] and k[1] >= h[1] and len(k[2]) > len(h[2])
        for k in hits)], key=lambda h: (h[0], h[1]))


def dinuc_shuffle(seq, rng):
    """Randomized Euler traversal: exact dinucleotide counts and fixed ends.

    This is NOT a uniform draw over all possible dinucleotide-preserving strings;
    empirical P-values are conditional, approximate screening diagnostics.
    """
    if len(seq) < 3:
        return seq
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    for v in edges.values():
        rng.shuffle(v)
    stack, trail = [seq[0]], []
    while stack:
        if edges[stack[-1]]:
            stack.append(edges[stack[-1]].pop())
        else:
            trail.append(stack.pop())
    return ''.join(trail[::-1])


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(1., adjusted)
    return out


def windows(genome, boundary0, strand, offsets, circular=True):
    """boundary0 = first/last transcribed base; offsets follow transcript direction."""
    positions = [boundary0 + strand * x for x in offsets]
    if circular:
        positions = [x % len(genome) for x in positions]
    else:
        positions = [x for x in positions if 0 <= x < len(genome)]
    seq = ''.join(genome[x] if strand == 1 else genome[x].translate(COMP) for x in positions)
    return seq, positions


def path_fields(positions, n):
    """Preserve transcription-order coordinates, including circular origin wraps."""
    if not positions:
        return {'genomic_start_1b': None, 'genomic_end_1b': None,
                'genomic_segments_1b': '', 'position_path_1b': '', 'wraps_origin': False}
    runs = [[positions[0]]]
    for x in positions[1:]:
        if abs(x - runs[-1][-1]) == 1:
            runs[-1].append(x)
        else:
            runs.append([x])
    return {'genomic_start_1b': min(positions)+1, 'genomic_end_1b': max(positions)+1,
            'genomic_segments_1b': ';'.join(f'{min(g)+1}..{max(g)+1}' for g in runs),
            'position_path_1b': ','.join(str(x+1) for x in positions),
            'wraps_origin': any(abs(a-b) == n-1 for a,b in zip(positions,positions[1:]))}


def make_annotations(record):
    genome = str(record.seq).upper()
    feature_rows, product_rows, anchor_rows, qc, translations = [], [], [], [], []
    exon_occupancy = [set() for _ in genome]
    cds_occupancy = np.zeros(len(genome), dtype=bool)
    for index, f in enumerate(record.features, 1):
        q = f.qualifiers
        parts = list(f.location.parts) if f.location is not None else []
        feature_rows.append({
            'feature_id': f'F{index:04d}', 'feature_type': f.type,
            'gene': q.get('gene',[''])[0], 'locus_tag': q.get('locus_tag',[''])[0],
            'product': q.get('product',[''])[0], 'location': str(f.location),
            'genomic_start_1b': int(f.location.start)+1 if parts else None,
            'genomic_end_1b': int(f.location.end) if parts else None,
            'strand': f.location.strand if parts else None,
            'extracted_length_bp': len(f) if parts else None, 'n_parts': len(parts),
            'trans_splicing': 'trans_splicing' in q, 'pseudo': 'pseudo' in q or 'pseudogene' in q,
            'qualifiers_json': json.dumps(q, ensure_ascii=False)})
        if f.type not in PRODUCT_TYPES or not parts:
            continue
        pid = f'P{len(product_rows)+1:04d}'
        gene = q.get('gene', q.get('locus_tag', [pid]))[0]
        strands = {p.strand for p in parts}
        mixed = len(strands) > 1
        split_contexts = mixed or 'trans_splicing' in q
        prod = {'product_id':pid, 'feature_id':f'F{index:04d}', 'gene':gene,
                'feature_type':f.type, 'product':q.get('product',[''])[0],
                'locus_tag':q.get('locus_tag',[''])[0], 'location':str(f.location),
                'n_parts':len(parts), 'spliced_length_bp':len(f),
                'mixed_strand':mixed, 'trans_splicing':'trans_splicing' in q,
                'pseudo':'pseudo' in q or 'pseudogene' in q,
                'codon_start':int(q.get('codon_start',['1'])[0]),
                'translation_table':int(q.get('transl_table',['11'])[0])}
        product_rows.append(prod)
        for p in parts:
            for x in range(int(p.start), int(p.end)):
                exon_occupancy[x].add(pid)
                if f.type == 'CDS': cds_occupancy[x] = True
        if mixed:
            qc.append({'severity':'WARNING','item':pid,'gene':gene,
                       'message':'Mixed-strand trans-spliced annotation: separate exon contexts; exon order taken from GenBank, not experimentally verified.'})
        if any(str(p.start).startswith(('<','>')) or str(p.end).startswith(('<','>')) for p in parts):
            qc.append({'severity':'WARNING','item':pid,'gene':gene,'message':'Fuzzy/partial feature boundary; boundary-based predictions uncertain.'})
        # For a cis-spliced transcript, location.parts are already biological order.
        contexts = [[p] for p in parts] if split_contexts else [parts]
        for k, context in enumerate(contexts):
            strand = context[0].strand
            if strand not in (1,-1):
                qc.append({'severity':'WARNING','item':pid,'gene':gene,'message':'Unknown strand: context skipped.'})
                continue
            first,last = context[0], context[-1]
            first0 = int(first.start) if strand == 1 else int(first.end)-1
            last0 = int(last.end)-1 if strand == 1 else int(last.start)
            aid = f'A{len(anchor_rows)+1:04d}'
            # Codon_start > 1 indicates a partial CDS; do not invent its RBS.
            sd_eligible = f.type == 'CDS' and k == 0 and prod['codon_start'] == 1 and not prod['pseudo']
            anchor_rows.append({**prod, 'anchor_id':aid,'strand':strand,
                'context_index':k+1,'context_kind':'trans_spliced_exon' if split_contexts else 'product_boundary',
                'context_start_1b':min(int(p.start) for p in context)+1,
                'context_end_1b':max(int(p.end) for p in context),
                'first_transcribed_base_1b':first0+1,'last_transcribed_base_1b':last0+1,
                'sd_eligible':sd_eligible})
        if f.type == 'CDS':
            seq = str(f.extract(record.seq)).upper()
            seq = seq[prod['codon_start']-1:]
            observed = q.get('translation',[''])[0]
            table = prod['translation_table']
            raw = str(Seq(seq[:len(seq)//3*3]).translate(table=table))
            normalized = raw.rstrip('*')
            if seq[:3] in ('ATG','GTG','TTG') and normalized:
                normalized = 'M'+normalized[1:]
            row = {'product_id':pid,'gene':gene,'translation_table':table,
                   'codon_start':prod['codon_start'],'coding_length_bp':len(seq),
                   'length_multiple_of_3':len(seq)%3 == 0,'start_codon':seq[:3],
                   'stop_codon':seq[-3:],'internal_stop_count':raw[:-1].count('*'),
                   'matches_annotation': normalized == observed,
                   'annotated_aa_length':len(observed),'computed_aa_length':len(normalized),
                   'annotated_translation':observed,'computed_translation':normalized}
            translations.append(row)
            if not row['matches_annotation'] or row['internal_stop_count'] or len(seq)%3:
                qc.append({'severity':'WARNING','item':pid,'gene':gene,
                           'message':'CDS translation/boundary check failed; inspect annotation before interpreting start-site motifs.'})
    return (pd.DataFrame(feature_rows), pd.DataFrame(product_rows), pd.DataFrame(anchor_rows),
            pd.DataFrame(qc,columns=['severity','item','gene','message']),
            pd.DataFrame(translations), exon_occupancy, cds_occupancy)


def hairpins(seq, stem_min=4, stem_max=12, loop_min=3, loop_max=10):
    """Exact antiparallel Watson-Crick hairpins, maximal outer stem per loop.

    GU pairs, bulges, pseudoknots and imperfect stems are NOT enumerated. Whole
    downstream-window MFE structures are also exported to partially expose this
    blind spot, but are not interpreted as separate regulatory-element calls.
    """
    found = {}
    n = len(seq)
    for looplen in range(loop_min,loop_max+1):
        for a in range(n - 2*stem_min - looplen + 1):
            l = a + stem_min
            r = l + looplen
            if seq[a:l] != rc(seq[r:r+stem_min]) or set(seq[a:l]) - set('ACGT'):
                continue
            x,y,s = a,r+stem_min,stem_min
            while s < stem_max and x > 0 and y < n and seq[x-1] == seq[y].translate(COMP) and seq[x-1] in 'ACGT':
                x -= 1; y += 1; s += 1
            key = (l,r)
            h = (x,y,s,looplen,l,r)
            if key not in found or s > found[key][2]: found[key] = h
    return sorted(found.values())


def fold_hairpin(seq, hit, md):
    a,b,stem,loop,l,r = hit
    core = seq[a:b].replace('T','U')
    intended = '('*stem + '.'*loop + ')'*stem
    fc = RNA.fold_compound(core,md)
    structure,mfe = fc.mfe()
    energy = fc.eval_structure(intended)
    pt = RNA.ptable(structure)
    retained = sum(pt[i+1] == len(core)-i for i in range(stem))/stem
    return core,intended,structure,float(mfe),float(energy),retained


def annotate_hit(anchor, category, seq, positions, a, b, occupancy, cds_mask, n):
    pos = positions[a:b]
    other = set().union(*(occupancy[x] for x in pos)) - {anchor['product_id']} if pos else set()
    return {'anchor_id':anchor['anchor_id'],'product_id':anchor['product_id'],
            'gene':anchor['gene'],'feature_type':anchor['feature_type'],
            'strand':anchor['strand'],'context_kind':anchor['context_kind'],
            'element_type':category,'sequence_5to3':seq[a:b], 'length_bp':b-a,
            **path_fields(pos,n),'overlap_CDS_bp':int(sum(cds_mask[x] for x in pos)),
            'overlap_other_products':';'.join(sorted(other)),
            'evidence':'computational sequence candidate; not experimentally validated'}


def analyze(record, anchors, occupancy, cds_mask, args):
    genome = str(record.seq).upper(); n=len(genome)
    circular = record.annotations.get('topology','linear').lower() == 'circular'
    md = RNA.md(); md.temperature = args.temperature
    rng = random.Random(args.seed)
    promoter_rows, raw_rows, sd_rows, hp_rows, window_rows, null_rows, full_folds = [],[],[],[],[],[],[]
    for ix, anchor in enumerate(anchors.to_dict('records')):
        first0,last0,strand = anchor['first_transcribed_base_1b']-1,anchor['last_transcribed_base_1b']-1,anchor['strand']
        up,uppos = windows(genome,first0,strand,range(-args.upstream,0),circular)
        down,downpos = windows(genome,last0,strand,range(1,args.downstream+1),circular)
        short,shortpos = windows(genome,first0,strand,range(-35,0),circular)
        for label,s,p in [('upstream',up,uppos),('downstream',down,downpos),('SD_window',short,shortpos)]:
            window_rows.append({'anchor_id':anchor['anchor_id'],'product_id':anchor['product_id'],
                'gene':anchor['gene'],'strand':strand,'window':label,'sequence_5to3':s,
                'length_bp':len(s),'GC_percent':100*(s.count('G')+s.count('C'))/len(s) if s else 0,
                **path_fields(p,n)})
        pairs = promoter_pairs(up)
        octamers = list(OCT.finditer(up)); hexamers = list(HEX.finditer(up))
        for a,b,j,spacer,mm35,mm10,score in pairs:
            row = annotate_hit(anchor,'PEP_paired_boxes',up,uppos,a,b,occupancy,cds_mask,n)
            row.update({'minus35_sequence':up[a:a+6],'minus10_sequence':up[j:j+6],
                'minus35_mismatches':mm35,'minus10_mismatches':mm10,'spacer_bp':spacer,
                'heuristic_score':score,'relative_start_to_boundary':a-len(up),
                'relative_end_to_boundary':b-1-len(up),
                'priority':'sequence_prioritized' if mm35+mm10 <= 1 and 16 <= spacer <= 18 and row['overlap_CDS_bp']==0 else 'exploratory',
                'statistical_test':'PEP_max_score'})
            promoter_rows.append(row)
        # Preserve isolated algal promoter motifs: -35 is not universally required.
        for matches,label,width,test in [(octamers,'algal_TATAATAT',8,'octamer_count'),(hexamers,'isolated_TATAAT',6,'hexamer_count')]:
            for m in matches:
                a,b=m.start(),m.start()+width
                row=annotate_hit(anchor,label,up,uppos,a,b,occupancy,cds_mask,n)
                row.update({'relative_start_to_boundary':a-len(up),'relative_end_to_boundary':b-1-len(up),
                            'priority':'exploratory','heuristic_score':width,'statistical_test':test})
                promoter_rows.append(row)
        # Diagnostic short words are explicitly EXCLUDED from the candidate inventory.
        for name,pat,width in [('PEP_minus35_near',P35,6),('PEP_minus10_near',P10,6),('YRTA_exploratory_NOT_NEP_call',YRTA,4),('GAA_background_NOT_promoter',GAA,3)]:
            for m in pat.finditer(up):
                row=annotate_hit(anchor,name,up,uppos,m.start(),m.start()+width,occupancy,cds_mask,n)
                row['relative_start_to_boundary']=m.start()-len(up)
                raw_rows.append(row)
        sd = sd_hits(short) if anchor['sd_eligible'] else []
        for a,b,motif,spacer,score in sd:
            row=annotate_hit(anchor,'SD_like',short,shortpos,a,b,occupancy,cds_mask,n)
            row.update({'motif':motif,'spacer_to_start_codon_bp':spacer,
                'motif_start_relative_to_CDS':a-len(short),'heuristic_score':score,
                'priority':'exploratory','statistical_test':'SD_max_length'})
            sd_rows.append(row)
        hs=[]
        for h in hairpins(down):
            core,intended,structure,mfe,energy,retained=fold_hairpin(down,h,md)
            if energy > -2: continue
            a,b,stem,loop,l,r=h
            row=annotate_hit(anchor,'downstream_hairpin',down,downpos,a,b,occupancy,cds_mask,n)
            tail=down[b:b+8]
            row.update({'left_arm':down[a:l],'loop':down[l:r],'right_arm':down[r:b],
                'stem_pairs':stem,'loop_length':loop,'MFE_kcal_mol':mfe,'intended_energy_kcal_mol':energy,
                'intended_structure':intended,'MFE_structure':structure,'MFE_intended_pair_fraction':retained,
                'tail_8nt':tail,'tail_T_count':tail.count('T'),
                'tail_T_fraction':tail.count('T')/len(tail) if tail else 0,
                'distance_after_boundary_bp':a+1,'relative_end_after_boundary_bp':b,
                'priority':'sequence_prioritized' if energy<=-5 and stem>=6 and retained>=.75 and row['overlap_CDS_bp']==0 else 'exploratory',
                'interpretation':'3prime structure/processing/stability candidate; not a proven transcription terminator',
                'statistical_test':'hairpin_min_intended_energy','selected_nonredundant':False})
            hs.append(row)
        # Keep all maximal alternatives; provide a deterministic non-overlapping subset.
        occupied=set()
        for row in sorted(hs,key=lambda x:(x['intended_energy_kcal_mol'],-x['stem_pairs'],x['distance_after_boundary_bp'])):
            span=set(range(row['distance_after_boundary_bp'],row['relative_end_after_boundary_bp']+1))
            if not span & occupied:
                row['selected_nonredundant']=True; occupied |= span
        hp_rows.extend(hs)
        if down:
            structure, mfe=RNA.fold_compound(down.replace('T','U'),md).mfe()
            full_folds.append({'anchor_id':anchor['anchor_id'],'gene':anchor['gene'],
                'RNA_sequence':down.replace('T','U'),'MFE_structure':structure,'MFE_kcal_mol':float(mfe),
                'temperature_C':args.temperature,'interpretation':'whole-window fold, not a regulatory-element annotation'})
        observed={'PEP_max_score':max([x[-1] for x in pairs],default=0.),
                  'octamer_count':len(octamers),'hexamer_count':len(hexamers)}
        if anchor['sd_eligible']: observed['SD_max_length']=max([x[-1] for x in sd],default=0.)
        nulls={k:[] for k in observed}
        for _ in range(args.null_reps):
            su=dinuc_shuffle(up,rng)
            nulls['PEP_max_score'].append(max([x[-1] for x in promoter_pairs(su)],default=0.))
            nulls['octamer_count'].append(sum(1 for _ in OCT.finditer(su)))
            nulls['hexamer_count'].append(sum(1 for _ in HEX.finditer(su)))
            if anchor['sd_eligible']:
                ss=dinuc_shuffle(short,rng)
                nulls['SD_max_length'].append(max([x[-1] for x in sd_hits(ss)],default=0.))
        for test,obs in observed.items():
            vals=np.array(nulls[test],dtype=float)
            null_rows.append({'anchor_id':anchor['anchor_id'],'product_id':anchor['product_id'],'gene':anchor['gene'],
                'test':test,'observed_statistic':obs,'null_mean':float(vals.mean()) if len(vals) else np.nan,
                'null_sd':float(vals.std(ddof=1)) if len(vals)>1 else np.nan,
                'null_replicates':args.null_reps,
                'empirical_P':(1+int((vals>=obs).sum()))/(len(vals)+1) if len(vals) else np.nan,
                'null_model':'randomized Euler dinucleotide-preserving shuffle; approximate, nonuniform',
                'test_scope':'per-anchor maximum or count; NOT a per-hit P-value'})
        if args.fold_null_reps:
            obs=min([x['intended_energy_kcal_mol'] for x in hs],default=0.)
            vals=[]
            for _ in range(args.fold_null_reps):
                ds=dinuc_shuffle(down,rng)
                energies=[fold_hairpin(ds,h,md)[4] for h in hairpins(ds)]
                vals.append(min([x for x in energies if x<=-2],default=0.))
            null_rows.append({'anchor_id':anchor['anchor_id'],'product_id':anchor['product_id'],'gene':anchor['gene'],
                'test':'hairpin_min_intended_energy','observed_statistic':obs,
                'null_mean':float(np.mean(vals)),'null_sd':float(np.std(vals,ddof=1)) if len(vals)>1 else np.nan,
                'null_replicates':args.fold_null_reps,
                'empirical_P':(1+sum(x<=obs for x in vals))/(len(vals)+1),
                'null_model':'randomized Euler dinucleotide-preserving shuffle; approximate, nonuniform',
                'test_scope':'per-anchor strongest enumerated perfect-stem structure; NOT a per-hit P-value'})
        if (ix+1)%10==0 or ix==0:
            print(f'Analyzed context {ix+1}/{len(anchors)}: {anchor["gene"]}',flush=True)
    null_df=pd.DataFrame(null_rows)
    null_df['BH_q_within_test_family']=np.nan
    for _,g in null_df.groupby('test'):
        null_df.loc[g.index,'BH_q_within_test_family']=bh_adjust(g.empirical_P)
    frames=[]
    for rows,prefix in [(promoter_rows,'PR'),(sd_rows,'SD'),(hp_rows,'HP')]:
        df=pd.DataFrame(rows)
        if not df.empty:
            df.insert(0,'association_id',[f'{prefix}{i+1:05d}' for i in range(len(df))])
            st=null_df.rename(columns={'test':'statistical_test','empirical_P':'anchor_empirical_P',
                'BH_q_within_test_family':'anchor_BH_q'})
            df=df.merge(st[['anchor_id','statistical_test','anchor_empirical_P','anchor_BH_q']],
                        how='left',on=['anchor_id','statistical_test'],validate='many_to_one')
        else:
            # Explicit schemas keep zero-hit categories exportable and plottable.
            common=['association_id','anchor_id','product_id','gene','feature_type','strand',
                'context_kind','element_type','sequence_5to3','length_bp','genomic_start_1b',
                'genomic_end_1b','genomic_segments_1b','position_path_1b','wraps_origin',
                'overlap_CDS_bp','overlap_other_products','evidence','priority','statistical_test',
                'anchor_empirical_P','anchor_BH_q']
            extra={'PR':['relative_start_to_boundary','relative_end_to_boundary','spacer_bp','heuristic_score'],
                   'SD':['spacer_to_start_codon_bp','motif','heuristic_score'],
                   'HP':['selected_nonredundant','stem_pairs','MFE_kcal_mol','intended_energy_kcal_mol',
                         'distance_after_boundary_bp','tail_T_fraction']}[prefix]
            df=pd.DataFrame(columns=common+extra)
        frames.append(df)
    return dict(zip(['Promoters','SD_like','Hairpins'],frames)) | {
        'Raw_motif_diagnostics':pd.DataFrame(raw_rows),'Sequence_windows':pd.DataFrame(window_rows),
        'Null_tests':null_df,'Downstream_MFE_folds':pd.DataFrame(full_folds)}


def neighborhoods(anchors, n, gapmax=200, circular=True):
    """Conservative physical-neighbor graph; opposite strands block connections.

    Overlaps are retained down to -50 bp; larger overlaps are not joined. Exon
    contexts of trans-spliced genes are separate nodes, never a genome-wide box.
    """
    rows=anchors.sort_values(['context_start_1b','context_end_1b']).to_dict('records')
    edges=[]; parent={r['anchor_id']:r['anchor_id'] for r in rows}
    def root(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    for i in range(len(rows) if circular else max(0,len(rows)-1)):
        a,b=rows[i],rows[(i+1)%len(rows)]
        wrap=i==len(rows)-1
        gap=b['context_start_1b']+(n if wrap else 0)-a['context_end_1b']-1
        joined=a['strand']==b['strand'] and -50<=gap<=gapmax
        if joined: parent[root(b['anchor_id'])]=root(a['anchor_id'])
        edges.append({'left_anchor':a['anchor_id'],'right_anchor':b['anchor_id'],
            'left_gene':a['gene'],'right_gene':b['gene'],'intergenic_gap_bp':gap,
            'same_strand':a['strand']==b['strand'],'crosses_origin':wrap,'joined':joined})
    groups=defaultdict(list)
    for r in rows:groups[root(r['anchor_id'])].append(r)
    output=[]
    for i,g in enumerate(groups.values(),1):
        ids={r['anchor_id'] for r in g}
        ingroup=[e for e in edges if e['joined'] and e['left_anchor'] in ids and e['right_anchor'] in ids]
        # Establish circular genomic order, then reverse for minus transcription.
        wrap=any(e['crosses_origin'] for e in ingroup)
        if wrap and len(g)>1:
            cuts=[j for j in range(len(g)) if not any(e['left_anchor']==g[j-1]['anchor_id'] and e['right_anchor']==g[j]['anchor_id'] for e in ingroup)]
            if cuts: j=cuts[0];g=g[j:]+g[:j]
        if g[0]['strand']==-1:g=g[::-1]
        output.append({'neighborhood_id':f'NB{i:03d}','n_contexts':len(g),
            'n_distinct_products':len({r['product_id'] for r in g}),
            'anchor_ids_transcription_order':';'.join(r['anchor_id'] for r in g),
            'genes_transcription_order':';'.join(r['gene'] for r in g),
            'strand':g[0]['strand'],'crosses_origin':wrap,
            'max_intergenic_gap_bp':max([e['intergenic_gap_bp'] for e in ingroup],default=None),
            'interpretation':'candidate co-transcription neighborhood' if len(g)>1 else 'singleton context; NOT an operon',
            'evidence':f'physical adjacency, same strand, -50..{gapmax} bp gap; no RNA evidence'})
    return pd.DataFrame(output),pd.DataFrame(edges)


def unique_sites(inventory):
    if inventory.empty:return inventory.copy()
    rows=[]
    for _,g in inventory.groupby(['element_type','strand','position_path_1b'],sort=False):
        r=g.iloc[0].to_dict()
        for col in ['association_id','anchor_id','product_id','gene']:
            r[col+'_all']=';'.join(dict.fromkeys(g[col].astype(str)))
        r['n_anchor_associations']=len(g)
        rows.append(r)
    out=pd.DataFrame(rows)
    out.insert(0,'site_id',[f'SITE{i+1:05d}' for i in range(len(out))])
    return out


def gene_summary(products, anchors, data):
    out=products.copy()
    count=lambda frame: frame.groupby('product_id').size() if len(frame) else pd.Series(dtype=int)
    out['n_contexts']=out.product_id.map(count(anchors)).fillna(0).astype(int)
    for name,key in [('PEP_pairs','PEP_paired_boxes'),('algal_octamers','algal_TATAATAT'),('isolated_TATAAT','isolated_TATAAT')]:
        df=data['Promoters'];sub=df[df.element_type==key] if len(df) else df
        out[name]=out.product_id.map(count(sub)).fillna(0).astype(int)
    for name,key in [('SD_like_sites','SD_like'),('hairpin_alternatives','Hairpins')]:
        out[name]=out.product_id.map(count(data[key])).fillna(0).astype(int)
    hp=data['Hairpins']
    sub=hp[hp.selected_nonredundant] if len(hp) else hp
    out['nonredundant_hairpins']=out.product_id.map(count(sub)).fillna(0).astype(int)
    if len(hp):out['best_hairpin_intended_energy']=out.product_id.map(hp.groupby('product_id').intended_energy_kcal_mol.min())
    out['interpretation']='motif absence is not evidence of absent transcription/translation/regulation'
    return out


def genome_regions(record):
    rows=[]
    for f in record.features:
        if f.type not in {'misc_feature','repeat_region'}:continue
        seq=str(f.extract(record.seq)).upper()
        rows.append({'feature_type':f.type,'description':';'.join(f.qualifiers.get('note',[])),
            'start_1b':int(f.location.start)+1,'end_1b':int(f.location.end),'length_bp':len(f),
            'GC_percent':100*(seq.count('G')+seq.count('C'))/len(seq)})
    return pd.DataFrame(rows)


def coordinate_qc(inventory, genome):
    rows=[]
    for r in inventory.to_dict('records'):
        pos=[int(x)-1 for x in r['position_path_1b'].split(',')]
        recovered=''.join(genome[x] if r['strand']==1 else genome[x].translate(COMP) for x in pos)
        rows.append({'association_id':r['association_id'],'sequence_length_match':len(pos)==r['length_bp'],
            'coordinates_in_range':all(0<=x<len(genome) for x in pos),
            'sequence_round_trip_pass':recovered==r['sequence_5to3']})
    return pd.DataFrame(rows)


def export_bed_gff(inventory, record, outdir):
    """Export unique category-specific sites, splitting circular-origin hits."""
    bed=[];gff=['##gff-version 3',f'##sequence-region {record.id} 1 {len(record)}']
    for r in unique_sites(inventory).to_dict('records'):
        segments=r['genomic_segments_1b'].split(';')
        for j,seg in enumerate(segments,1):
            a,b=map(int,seg.split('..'));strand='+' if r['strand']==1 else '-'
            name=f'{r["site_id"]}.{j}'
            bed.append(f'{record.id}\t{a-1}\t{b}\t{name}\t0\t{strand}')
            products=r['product_id_all'].replace(';',',')
            attrs=f'ID={name};Name={r["element_type"]};Note=computational_candidate_not_validated;associated_products={products}'
            gff.append(f'{record.id}\tchloroplast_inventory\tsequence_feature\t{a}\t{b}\t.\t{strand}\t.\t{attrs}')
    (outdir/'regulatory_candidates.bed').write_text('\n'.join(bed)+'\n')
    (outdir/'regulatory_candidates.gff3').write_text('\n'.join(gff)+'\n')


def run(args):
    out=Path(args.outdir);results=out/'results';figures=out/'figures';logs=out/'logs'
    for p in [results,figures,logs]:p.mkdir(parents=True,exist_ok=True)
    if args.genbank:
        gb=Path(args.genbank);record=SeqIO.read(gb,'genbank')
    else:
        if not args.email or args.email.endswith('@example.com'):
            raise ValueError('Use --genbank for offline input or supply a real --email for NCBI Entrez.')
        Entrez.email=args.email;Entrez.tool='chloroplast_regulatory_inventory'
        if args.api_key:Entrez.api_key=args.api_key
        with Entrez.efetch(db='nuccore',id=args.accession,rettype='gbwithparts',retmode='text') as h:
            text=h.read()
        (out/'inputs').mkdir(exist_ok=True);gb=out/'inputs'/f'{args.accession}.gb';gb.write_text(text)
        record=SeqIO.read(gb,'genbank')
    if len(record)==0:raise ValueError('Empty genome')
    f,p,a,qc,tr,occ,cdsmask=make_annotations(record)
    if a.empty:
        raise ValueError('No oriented CDS/RNA product features found. Supply an annotated GenBank file, not an unannotated FASTA.')
    data=analyze(record,a,occ,cdsmask,args)
    inventory=pd.concat([x for x in [data['Promoters'],data['SD_like'],data['Hairpins']] if len(x)],ignore_index=True,sort=False) if any(len(data[k]) for k in ['Promoters','SD_like','Hairpins']) else data['Promoters'].copy()
    nb,edges=neighborhoods(a,len(record),args.operon_gap,record.annotations.get('topology')=='circular')
    sites=unique_sites(inventory)
    coord=coordinate_qc(inventory,str(record.seq).upper())
    if len(coord) and not coord.iloc[:,1:].all().all():raise AssertionError('Coordinate round-trip failed')
    genome=str(record.seq).upper()
    region=genome_regions(record)
    repeats=[f for f in record.features if f.type=='repeat_region' and 'inverted' in f.qualifiers.get('rpt_type',[])]
    ir_identity=None
    if len(repeats)==2:
        x=str(repeats[0].extract(record.seq));y=rc(str(repeats[1].extract(record.seq)))
        ir_identity=sum(i==j for i,j in zip(x,y))/max(len(x),len(y))*100
    versions={x:importlib.metadata.version(x) for x in ['biopython','ViennaRNA','numpy','pandas','matplotlib','xlsxwriter','openpyxl']}
    meta={'accession':record.id,'description':record.description,'organism':record.annotations.get('organism',''),
        'genome_length_bp':len(record),'topology':record.annotations.get('topology',''),
        'record_date':record.annotations.get('date',''),'GC_percent':100*(genome.count('G')+genome.count('C'))/len(genome),
        'ambiguous_bases':sum(b not in 'ACGT' for b in genome),'n_annotated_features':len(record.features),
        'n_product_features':len(p),'n_distinct_gene_names':p.gene.nunique(),'n_boundary_contexts':len(a),
        'n_mixed_strand_products':int(p.mixed_strand.sum()),'IR_reverse_complement_identity_percent':ir_identity,
        'n_inventory_associations':len(inventory),'n_category_specific_unique_sites':len(sites),
        'analysis_UTC':datetime.now(timezone.utc).isoformat(),'pipeline_version':VERSION,
        'input_GenBank_sha256':hashlib.sha256(gb.read_bytes()).hexdigest(),
        'sequence_sha256':hashlib.sha256(genome.encode()).hexdigest(),'python':platform.python_version(),
        'software_versions':json.dumps(versions),'parameters':json.dumps(vars(args),default=str)}
    data.update({'Metadata':pd.DataFrame([{'field':k,'value':v} for k,v in meta.items()]),
        'Features':f,'Products':p,'Boundary_contexts':a,'Gene_summary':gene_summary(p,a,data),
        'Genome_regions':region,'Neighborhoods':nb,'Adjacency_edges':edges,
        'Inventory_associations':inventory,'Unique_sites':sites,'Annotation_QC':qc,
        'Translation_QC':tr,'Coordinate_QC':coord})
    from inventory_reporting import add_summary, create_workbook, create_figures
    add_summary(data)
    if (out/'inputs'/'baseline_original_summary.json').exists():
        base=json.loads((out/'inputs'/'baseline_original_summary.json').read_text())
        data['Original_baseline']=pd.DataFrame([{'metric':k,'original_value':v} for k,v in base.items()])
    for name,df in data.items():df.to_csv(results/f'{name}.tsv',sep='\t',index=False)
    SeqIO.write(record,results/f'{record.id}.fasta','fasta')
    export_bed_gff(inventory,record,results)
    (logs/'run_manifest.json').write_text(json.dumps(meta,indent=2,default=str))
    create_figures(data,record,figures,args.temperature)
    create_workbook(data,out,results/f'{record.id}_regulatory_inventory.xlsx')
    print(json.dumps({**{k:meta[k] for k in ['accession','genome_length_bp','GC_percent','n_product_features','n_boundary_contexts','n_inventory_associations','n_category_specific_unique_sites']},
        'PEP_pairs':int((data['Promoters'].element_type=='PEP_paired_boxes').sum()) if len(data['Promoters']) else 0,
        'SD_like':len(data['SD_like']),'hairpin_alternatives':len(data['Hairpins']),
        'multicontext_neighborhoods':int((nb.n_contexts>1).sum()),
        'q_lt_0_05_by_test':data['Null_tests'].assign(significant=lambda d:d.BH_q_within_test_family<.05).groupby('test').significant.sum().to_dict()},indent=2),flush=True)
    return data


def parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--genbank',help='Local annotated GenBank file (preferred for reproducibility).')
    p.add_argument('--accession',default='QB176266.1');p.add_argument('--email');p.add_argument('--api-key')
    p.add_argument('--outdir',default='chloroplast_analysis')
    p.add_argument('--upstream',type=int,default=300);p.add_argument('--downstream',type=int,default=120)
    p.add_argument('--operon-gap',type=int,default=200)
    p.add_argument('--null-reps',type=int,default=1000);p.add_argument('--fold-null-reps',type=int,default=200)
    p.add_argument('--temperature',type=float,default=25.0);p.add_argument('--seed',type=int,default=20261023)
    return p


if __name__=='__main__':
    args=parser().parse_args()
    if args.null_reps < 1 or args.fold_null_reps < 0:raise ValueError('null-reps >=1 and fold-null-reps >=0 required')
    run(args)
