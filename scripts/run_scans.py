#!/usr/bin/env python3
"""Run the uniform regulatory scan on every standardised plastome (parallel)."""
import argparse
import pickle
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from Bio import SeqIO

import cpreg_scan as cs

FULL_HAIRPINS = {'QB176266.1', 'FJ423446'}


def work(job):
    acc, path, out, a = job
    f = out / f'{acc}.pkl'
    if f.exists():
        return acc, 'cached'
    t = time.time()
    rec = SeqIO.read(path, 'genbank')
    res = cs.scan_genome(rec, acc, nrep=a.nrep, nrep_fold=a.nrep_fold, nrep_start=a.nrep_start,
                         temperature=a.temperature, seed=a.seed, keep_all_hairpins=acc in FULL_HAIRPINS)
    pickle.dump(res, open(f, 'wb'))
    return acc, f'{time.time()-t:.0f}s'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project', default='..')
    ap.add_argument('--nrep', type=int, default=200)
    ap.add_argument('--nrep-fold', type=int, default=100)
    ap.add_argument('--nrep-start', type=int, default=50)
    ap.add_argument('--temperature', type=float, default=25.0)
    ap.add_argument('--seed', type=int, default=20261023)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--tag', default='main')
    a = ap.parse_args()
    P = Path(a.project)
    out = P / 'results' / 'scan_cache' / a.tag; out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(P/'results/tables/genome_manifest.tsv', sep='\t')
    man = man[man.included_in_comparative_set]
    # focal and benchmark genomes first
    order = sorted(man.insdc_accession, key=lambda x: (x not in FULL_HAIRPINS, x))
    jobs = [(acc, P/'data/standardized'/f'{acc}.gb', out, a) for acc in order]
    with Pool(a.threads) as pool:
        for acc, msg in pool.imap_unordered(work, jobs):
            print(acc, msg, flush=True)


if __name__ == '__main__':
    main()
