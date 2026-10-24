#!/usr/bin/env python3
"""Benchmark the uniform scan against experimentally validated C. reinhardtii elements.

Validated elements (data/validated_elements_curated.tsv, primary literature)
are mapped onto the C. reinhardtii plastome FJ423446 in transcript orientation.
Where an exactly reported 5'UTR element has a reported TSS-relative position
(Anthonisen et al. 2001: rbcL +38..+47, atpB +31..+42) or a 5'UTR length is
reported (psbD 74-nt leader, Nickelsen et al. 1994), the TSS is inferred and the
scan's TATAATAT calls are evaluated for the expected -10 placement
(first octamer base at TSS-15..TSS-5). Other classes are evaluated as
detected / not detected by the uniform scan parameters (SD-like words, 3' window
MFE z-score and BH q, perfect-stem hairpin energy). The same element sequences
and positional rules are then projected onto every other genome.
"""
import argparse
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import binomtest

import chloroplast_regulatory_inventory_improved as v2

TSS_ANCHORS = {  # gene: (sequence in transcript DNA orientation, TSS-relative position of first base, source)
    'rbcL': ('ATTTCCGGAC', 38, 'CE042 Anthonisen et al. 2001 RNA 7:1024'),
    'atpB': ('ATAAGCGTTAGT', 31, 'CE043 Anthonisen et al. 2001 RNA 7:1024'),
}
LEADER = {'psbD': (74, 'CE016 Nickelsen et al. 1994 EMBO J (74-nt psbD leader)')}
THREE_PRIME = {  # gene: curated element ids reporting 3' IR function
    'atpB': 'CE050/CE051/CE053/CE055', 'rbcL': 'CE052/CE056', 'psaB': 'CE052/CE054',
    'petD': 'CE053', 'petA': 'CE053'}
SD_VALIDATED = {'psbD': ('GGAG', -14, 'CE018 Nickelsen et al. 1999 (PRB1 GGAG, -14/-11)'),
                'psbA': ('GGAG*', None, 'CE009 Mayfield et al. 1994 (SD element required)')}
SD_DISPENSABLE = {'atpB', 'atpE', 'rps4', 'rps7'}  # Fargo et al. 1998: SD-like present but dispensable


def near(a, b, mm):
    return len(a) == len(b) and sum(x != y for x, y in zip(a, b)) <= mm


def find_mm(seq, motif, mm):
    return [i for i in range(len(seq) - len(motif) + 1) if near(seq[i:i+len(motif)], motif, mm)]


def upstream(rec, gene, L=600):
    g = str(rec.seq).upper(); n = len(g)
    out = []
    for f in rec.features:
        if f.type in ('CDS', 'rRNA') and f.qualifiers['gene'][0] == gene and 'exclude_from_boundary_analysis' not in f.qualifiers:
            parts = list(f.location.parts)
            p = parts[0]
            first0 = int(p.start) if p.strand == 1 else int(p.end) - 1
            s, pos = v2.windows(g, first0, p.strand, range(-L, 0), True)
            out.append((s, pos, p.strand))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..'); ap.add_argument('--tag', default='main')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t'); man = man[man.included_in_comparative_set]
    cr = SeqIO.read(P/'data/standardized/FJ423446.gb', 'genbank')
    scan = pickle.load(open(P/f'results/scan_cache/{a.tag}/FJ423446.pkl', 'rb'))
    an = scan['anchors']
    rows = []
    # 1. Promoter -10 octamer placement relative to inferred TSS
    tss = {}
    for gene, (el, rel, src) in TSS_ANCHORS.items():
        for s, pos, strand in upstream(cr, gene)[:1]:
            hits = find_mm(s, el, 0)
            if hits:
                tss[gene] = (hits[0] - len(s)) - (rel - 1)
                rows.append(dict(element_class='TSS_inference', gene=gene, validated_source=src,
                                 observation=f'{el} found at {hits[0]-len(s)} nt from start codon; inferred TSS at {tss[gene]}',
                                 detected=True))
    for gene, (leader, src) in LEADER.items():
        tss[gene] = -leader
        rows.append(dict(element_class='TSS_inference', gene=gene, validated_source=src,
                         observation=f'TSS inferred at {-leader} from reported leader length', detected=True))
    for gene, t in tss.items():
        s, pos, strand = upstream(cr, gene)[0]
        octs = [m.start() - len(s) for m in re.finditer('(?=TATAATAT)', s)]
        ok = [o for o in octs if t - 15 <= o <= t - 5]
        rows.append(dict(element_class='promoter_-10_octamer', gene=gene,
                         validated_source='CE001 Klein et al. 1992 PNAS 89:3453 (-10 TATAATAT)',
                         observation=f'octamers at {octs}; expected first base in [{t-15},{t-5}] (TSS-15..TSS-5); '
                                     f'TSS-relative: {[o - t for o in octs]}',
                         detected=bool(ok), within_300nt_scan_window=bool(ok) and min(ok) >= -300))
    # 16S rRNA sigma70-type promoter (-35/-10) -- lenient pair search
    s, pos, strand = upstream(cr, 'rrn16')[0]
    w = s[-300:]
    lenient = []
    for i in range(len(w) - 6):
        m35 = sum(x != y for x, y in zip(w[i:i+6], 'TTGACA'))
        if m35 > 2:
            continue
        for sp in range(15, 20):
            j = i + 6 + sp
            if j + 6 <= len(w):
                m10 = sum(x != y for x, y in zip(w[j:j+6], 'TATAAT'))
                if m10 <= 2:
                    lenient.append((m35 + m10, abs(sp - 17), i - 300, w[i:i+6], sp, w[j:j+6]))
    lenient.sort()
    pairs = sorted(v2.promoter_pairs(s[-300:]), key=lambda h: -h[-1])
    rows.append(dict(element_class='promoter_-35/-10', gene='rrn16', validated_source='CE002 Klein et al. 1992',
                     observation=(f'best pair {s[-300:][pairs[0][0]:pairs[0][0]+6]}..{s[-300:][pairs[0][2]:pairs[0][2]+6]} '
                                  f'spacer {pairs[0][3]} mm {pairs[0][4]+pairs[0][5]} at {pairs[0][0]-300}') if pairs else 'no strict pair (<=1 mismatch per box)',
                     detected=bool(pairs and pairs[0][4]+pairs[0][5] <= 1),
                     detected_lenient=bool(lenient and lenient[0][0] <= 3),
                     lenient_best=(f'{lenient[0][3]}-[{lenient[0][4]}]-{lenient[0][5]} at {lenient[0][2]} ({lenient[0][0]} mm)' if lenient else '')))
    # 2. SD-like
    for gene, (motif, relpos, src) in SD_VALIDATED.items():
        r = an[an.gene == gene].iloc[0]
        up = upstream(cr, gene, 60)[0][0]
        far = [m.start() - 60 for m in re.finditer('(?=GGAG)', up)]
        rows.append(dict(element_class='SD_like_required', gene=gene, validated_source=src,
                         observation=f'scan SD calls (3-15 nt spacer rule): {r.SD_words or "none"}; GGAG within -60..-1 at {far}',
                         detected=bool(r.SD_like), detected_lenient=bool(far)))
    for gene in sorted(SD_DISPENSABLE):
        r = an[an.gene == gene]
        if len(r):
            r = r.iloc[0]
            rows.append(dict(element_class='SD_like_dispensable_info', gene=gene, validated_source='CE029 Fargo et al. 1998 MGG 257:271',
                             observation=f'scan SD calls: {r.SD_words or "none"}', detected=bool(r.SD_like)))
    # 3. 3' structure
    for gene, src in THREE_PRIME.items():
        r = an[an.gene == gene].iloc[0]
        rows.append(dict(element_class='3prime_IR', gene=gene, validated_source=src,
                         observation=f'down_z={r.down_z:.2f}; q={r.down_q_BH:.3f}; best perfect-stem dG={r.best_hairpin_energy}',
                         detected=bool(r.down_q_BH < .05), detected_lenient=bool(r.down_z <= -2 or r.best_hairpin_energy <= -10)))
    B = pd.DataFrame(rows)
    B.to_csv(T/'benchmark_Creinhardtii.tsv', sep='\t', index=False)
    summ = []
    for cls, g in B[~B.element_class.isin(['TSS_inference', 'SD_like_dispensable_info'])].groupby('element_class'):
        k, nn = int(g.detected.sum()), len(g)
        ci = binomtest(k, nn).proportion_ci(method='exact') if nn else None
        kl = int(g.get('detected_lenient', pd.Series(dtype=bool)).fillna(g.detected).astype(bool).sum()) if 'detected_lenient' in g else k
        summ.append(dict(element_class=cls, n_validated=nn, n_detected=k, recall=k/nn, n_detected_lenient=kl,
                         CI95_low=ci.low, CI95_high=ci.high))
    pd.DataFrame(summ).to_csv(T/'benchmark_summary.tsv', sep='\t', index=False)
    print(B.to_string()); print(pd.DataFrame(summ).to_string())
    # 4. Projection of validated sequence elements onto every genome
    proj = []
    for r in man.itertuples():
        rec = SeqIO.read(P/'data/standardized'/f'{r.insdc_accession}.gb', 'genbank')
        for gene, (el, rel, src) in list(TSS_ANCHORS.items()):
            for s, pos, strand in upstream(rec, gene)[:1]:
                h0, h1, h2 = find_mm(s, el, 0), find_mm(s, el, 1), find_mm(s, el, 2)
                proj.append(dict(genome=r.insdc_accession, species=r.species, clade=r.clade, gene=gene, element=el,
                                 exact=len(h0), le1mm=len(h1), le2mm=len(h2),
                                 positions_le2mm=';'.join(str(x - len(s)) for x in h2)))
        for gene in ('psbD',):
            for s, pos, strand in upstream(rec, gene)[:1]:
                w = s[-20:-5]
                proj.append(dict(genome=r.insdc_accession, species=r.species, clade=r.clade, gene=gene,
                                 element='GGAG (PRB1, -20..-6)', exact=w.count('GGAG'), le1mm=np.nan, le2mm=np.nan,
                                 positions_le2mm=';'.join(str(m.start() - 20) for m in re.finditer('(?=GGAG)', w))))
    pd.DataFrame(proj).to_csv(T/'validated_element_projection.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
