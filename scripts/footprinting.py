#!/usr/bin/env python3
"""De novo motif discovery / phylogenetic footprinting of upstream non-coding regions.

1. Per-gene footprinting (MEME, ZOOPS, given strand, 2nd-order Markov
   background estimated from all upstream non-coding sequences of the set):
   orthologous upstream non-coding stretches (<=300 nt, >=30 nt) of key genes.
2. Pooled discriminative discovery (STREME, shuffled controls preserving
   3-mer composition) of the proximal -100..-1 non-coding stretch, by clade.
3. TOMTOM comparison of every discovered motif against a reference library of
   experimentally validated chloroplast elements (validated_elements_curated.tsv).
"""
import argparse
import pickle
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

import os as _os
import shutil as _shutil
# External tools (MAFFT, trimAl, IQ-TREE, MEME Suite, BLAST+, ARAGORN): set CPREG_BIN to their bin
# directory, or put them on PATH (the directory containing 'mafft' is then used).
BIN = Path(_os.environ.get('CPREG_BIN') or (Path(_shutil.which('mafft')).parent if _shutil.which('mafft') else '/usr/local/bin'))
KEY_GENES = ['psbA', 'psbD', 'psbB', 'psbC', 'psbE', 'psbH', 'psbK', 'rbcL', 'atpA', 'atpB', 'atpE', 'atpH', 'atpI',
             'petA', 'petB', 'petD', 'petG', 'psaB', 'psaC', 'psaJ', 'rps4', 'rps7', 'rps12', 'rpl20', 'rpl23', 'rrn16',
             'tufA', 'chlL', 'chlN', 'ccsA', 'clpP', 'ycf3', 'ycf4', 'rpoA', 'ftsH', 'cemA']
VALIDATED = {  # id: (sequence, element, curated row)
    'CR_-10_octamer': ('TATAATAT', 'Chlamydomonas -10 octamer', 'CE001'),
    'PEP_-35': ('TTGACA', 'sigma70 -35 (16S rRNA promoter)', 'CE002'),
    'PEP_-10': ('TATAAT', 'sigma70 -10 hexamer', 'CE002'),
    'SD_GGAGG': ('GGAGG', 'Shine-Dalgarno consensus', 'CE029'),
    'psbD_PRB1': ('GGAG', 'psbD PRB1 SD-like', 'CE018'),
    'rbcL_stab': ('ATTTCCGGAC', 'rbcL 5UTR stability element +38/+47', 'CE042'),
    'atpB_stab': ('ATAAGCGTTAGT', 'atpB 5UTR stability element +31/+42', 'CE043'),
}


def meme_motif_file(path, motifs, bg=None):
    lines = ['MEME version 5', '', 'ALPHABET= ACGT', '', 'strands: +', '',
             'Background letter frequencies', 'A 0.35 C 0.15 G 0.15 T 0.35', '']
    for mid, (seq, name, src) in motifs.items():
        lines.append(f'MOTIF {mid} {src}')
        lines.append(f'letter-probability matrix: alength= 4 w= {len(seq)} nsites= 20 E= 0')
        for c in seq:
            p = [0.91 if c == b else 0.03 for b in 'ACGT']
            lines.append(' '.join(f'{x:.3f}' for x in p))
        lines.append('')
    Path(path).write_text('\n'.join(lines))


def parse_meme(xml):
    out = []
    root = ET.parse(xml).getroot()
    seqs = {s.get('id'): (s.get('name'), int(s.get('length'))) for s in root.iter('sequence')}
    for m in root.iter('motif'):
        sites = [(seqs[s.get('sequence_id')][0], int(s.get('position')) - seqs[s.get('sequence_id')][1])
                 for s in m.iter('contributing_site')]
        rel = [x for _, x in sites]
        af = [x for n, x in sites if n == 'QB176266.1']
        out.append(dict(motif=m.get('id'), consensus=m.get('name'), width=int(m.get('width')),
                        sites=int(m.get('sites')), E_value=float(m.get('e_value')), llr=float(m.get('llr')),
                        median_site_pos_rel_start=float(np.median(rel)) if rel else np.nan,
                        IQR_site_pos=(f'{np.percentile(rel,25):.0f}..{np.percentile(rel,75):.0f}' if rel else ''),
                        Afalcatus_site_pos=';'.join(map(str, af)), site_genomes=';'.join(n for n, _ in sites)))
    return out


def parse_streme(xml):
    out = []
    root = ET.parse(xml).getroot()
    for m in root.iter('motif'):
        out.append(dict(motif=m.get('id'), consensus=m.get('id').split('-', 1)[-1], width=int(m.get('width')),
                        sites=float(m.get('train_pos_count', 0)), E_value=float(m.get('test_evalue', m.get('evalue', 'nan'))),
                        score_threshold=m.get('score_threshold')))
    return out


def tomtom(query, db, oc):
    subprocess.run([str(BIN/'tomtom'), '-no-ssc', '-oc', str(oc), '-min-overlap', '4', '-dist', 'pearson',
                    '-evalue', '-thresh', '10', '-verbosity', '1', str(query), str(db)], check=True)
    f = Path(oc)/'tomtom.tsv'
    if not f.exists():
        return pd.DataFrame()
    return pd.read_csv(f, sep='\t', comment='#')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..'); ap.add_argument('--tag', default='main')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; W = P/'results/motifs'; W.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t'); man = man[man.included_in_comparative_set]
    clade = dict(zip(man.insdc_accession, man.clade)); fam = dict(zip(man.insdc_accession, man.family))
    an = []
    for acc in man.insdc_accession:
        r = pickle.load(open(P/f'results/scan_cache/{a.tag}/{acc}.pkl', 'rb'))
        an.append(r['anchors'][['genome', 'anchor_id', 'gene', 'feature_type', 'nc_upstream']])
    an = pd.concat(an)
    an = an.drop_duplicates(['genome', 'gene'])  # IR copies / exon contexts: first context per gene
    allseq = W/'all_upstream_nc.fa'
    allseq.write_text(''.join(f'>{r.genome}|{r.gene}\n{r.nc_upstream}\n' for r in an.itertuples() if len(r.nc_upstream) >= 30))
    bg = W/'background_markov2.txt'
    subprocess.run([str(BIN/'fasta-get-markov'), '-m', '2', '-norc', str(allseq), str(bg)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    db = W/'validated_reference_motifs.meme'
    meme_motif_file(db, VALIDATED)
    rows, tt = [], []
    # 1. per-gene footprinting
    for gene in KEY_GENES:
        g = an[(an.gene == gene) & (an.nc_upstream.str.len() >= 30)]
        if len(g) < 15:
            continue
        fa = W/f'fp_{gene}.fa'
        fa.write_text(''.join(f'>{r.genome}\n{r.nc_upstream}\n' for r in g.itertuples()))
        oc = W/f'meme_{gene}'
        subprocess.run([str(BIN/'meme'), str(fa), '-dna', '-oc', str(oc), '-mod', 'zoops', '-nmotifs', '3',
                        '-minw', '6', '-maxw', '20', '-bfile', str(bg), '-nostatus', '-searchsize', '0',
                        '-p', '1'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for m in parse_meme(oc/'meme.xml'):
            rows.append(dict(analysis='per_gene_MEME', group=gene, n_sequences=len(g),
                             includes_Afalcatus='QB176266.1' in set(g.genome), **m))
        t = tomtom(oc/'meme.txt', db, W/f'tomtom_{gene}')
        if len(t):
            tt.append(t.assign(group=gene, analysis='per_gene_MEME'))
    # 2. pooled STREME by group
    groups = {
        'A_falcatus': ['QB176266.1'],
        'Selenastraceae': [x for x in man.insdc_accession if fam[x] == 'Selenastraceae'],
        'Sphaeropleales': [x for x in man.insdc_accession if clade[x] == 'Sphaeropleales'],
        'Chlamydomonadales': [x for x in man.insdc_accession if clade[x] == 'Chlamydomonadales'],
        'C_reinhardtii': ['FJ423446'],
    }
    for gname, accs in groups.items():
        g = an[an.genome.isin(accs) & (an.nc_upstream.str.len() >= 30)]
        fa = W/f'streme_{gname}.fa'
        fa.write_text(''.join(f'>{r.genome}|{r.gene}\n{r.nc_upstream[-100:]}\n' for r in g.itertuples()))
        oc = W/f'streme_{gname}'
        subprocess.run([str(BIN/'streme'), '--p', str(fa), '--dna', '--oc', str(oc), '--nmotifs', '8',
                        '--minw', '6', '--maxw', '15', '--order', '2', '--seed', '20261023', '--verbosity', '1'],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for m in parse_streme(oc/'streme.xml'):
            rows.append(dict(analysis='pooled_STREME_-100..-1', group=gname, n_sequences=len(g),
                             includes_Afalcatus='QB176266.1' in accs, **m))
        t = tomtom(oc/'streme.txt', db, W/f'tomtom_streme_{gname}')
        if len(t):
            tt.append(t.assign(group=gname, analysis='pooled_STREME_-100..-1'))
    R = pd.DataFrame(rows); R.to_csv(T/'motif_discovery.tsv', sep='\t', index=False)
    TT = pd.concat(tt) if tt else pd.DataFrame(); TT.to_csv(T/'motif_tomtom_vs_validated.tsv', sep='\t', index=False)
    print(R.to_string()); print(TT.head(40).to_string())


if __name__ == '__main__':
    main()
