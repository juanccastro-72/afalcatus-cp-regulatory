#!/usr/bin/env python3
"""Window-length / masking sensitivity of genome-level promoter-motif enrichment.

For upstream windows of 150, 300 and 500 nt, with (a) non-coding masking (the
contiguous non-coding stretch only, as in the main scan) and (b) no masking
(full window, may include coding DNA), observed genome-wide TATAATAT, TATAAT,
TTGACA and strict PEP-pair counts are compared with 200 dinucleotide-shuffled
replicates of the same windows.
"""
import argparse
import random
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

import chloroplast_regulatory_inventory_improved as v2
import cpreg_scan as cs

WINDOWS = (150, 300, 500)


def one(job):
    acc, path, nrep = job
    rec = SeqIO.read(path, 'genbank')
    g = str(rec.seq).upper()
    _, _, anchors, _, _, occ, _ = cs.prepare(rec)
    noncoding = np.array([len(x) == 0 for x in occ])
    out = []
    for L in WINDOWS:
        for mask in ('noncoding', 'unmasked'):
            rng = random.Random(f'20261023-{acc}-{L}-{mask}')
            seqs = []
            for a in anchors.to_dict('records'):
                up, pos = v2.windows(g, a['first_transcribed_base_1b']-1, a['strand'], range(-L, 0), True)
                if mask == 'noncoding':
                    k = 0
                    for x in reversed(pos):
                        if not noncoding[x]:
                            break
                        k += 1
                    up = up[len(up)-k:]
                seqs.append(up)
            def stat(ss):
                return np.array([sum(cs.count(cs.OCT, s) for s in ss), sum(cs.count(cs.HEX, s) for s in ss),
                                 sum(cs.count(cs.M35, s) for s in ss), sum(len(cs.strong_pairs(s)) for s in ss)])
            o = stat(seqs)
            null = np.array([stat([v2.dinuc_shuffle(s, rng) if len(s) > 2 else s for s in seqs]) for _ in range(nrep)])
            for i, m in enumerate(('TATAATAT', 'TATAAT', 'TTGACA', 'PEP_pair')):
                mu = null[:, i].mean()
                out.append(dict(genome=acc, window=L, masking=mask, motif=m, observed=int(o[i]), expected=mu,
                                OE=o[i]/mu if mu else np.nan, total_bp=sum(map(len, seqs)),
                                P_enriched=(1 + (null[:, i] >= o[i]).sum())/(nrep+1),
                                P_depleted=(1 + (null[:, i] <= o[i]).sum())/(nrep+1)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..'); ap.add_argument('--nrep', type=int, default=200)
    ap.add_argument('--threads', type=int, default=2)
    a = ap.parse_args(); P = Path(a.project)
    man = pd.read_csv(P/'results/tables/genome_manifest.tsv', sep='\t'); man = man[man.included_in_comparative_set]
    jobs = [(acc, P/'data/standardized'/f'{acc}.gb', a.nrep) for acc in man.insdc_accession]
    rows = []
    with Pool(a.threads) as pool:
        for r in pool.imap_unordered(one, jobs):
            rows.extend(r); print(r[0]['genome'], flush=True)
    pd.DataFrame(rows).to_csv(P/'results/tables/promoter_window_sensitivity.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
