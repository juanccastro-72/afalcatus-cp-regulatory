#!/usr/bin/env python3
"""Align RNA-seq reads to the A. falcatus plastome (QB176266.1) with minimap2 (mappy, sr
preset) and compute strand-specific per-base coverage and 5'-end counts. One .npz per run.

Usage: python rnaseq_plastome_coverage.py --genome QB176266.1.gb --reads reads/ --out cov/
       [--max-reads N]   # cap reads per mate (0 = all)   [--workers K]
Requires: mappy, biopython, numpy. Tolerant of truncated (range-downloaded) gzip FASTQ.
"""
import argparse, glob, gzip, os
from pathlib import Path
from multiprocessing import Pool
import numpy as np, mappy as mp
from Bio import SeqIO

FASTA = '_plastome.fa'

def fq(path, cap):
    """Yield read sequences; stop cleanly if the gzip stream is truncated."""
    op = gzip.open if path.endswith('.gz') else open
    n = 0
    try:
        with op(path, 'rt') as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().strip(); fh.readline(); fh.readline()
                if not s:
                    break
                yield s
                n += 1
                if cap and n >= cap:
                    break
    except (EOFError, OSError):
        return  # truncated final gzip member -> use what we read

def process_run(args):
    run, reads_dir, out_dir, cap, genome_id, n = args
    al = mp.Aligner(FASTA, preset='sr')
    covp = np.zeros(n); covm = np.zeros(n); fivep = np.zeros(n); fivem = np.zeros(n)
    total = mapped = 0
    for mate in ('1', '2'):
        fp = f'{reads_dir}/{run}_{mate}.fastq.gz'
        if not os.path.exists(fp):
            continue
        for seq in fq(fp, cap):
            total += 1
            best = None
            for h in al.map(seq):
                if best is None or h.mlen > best.mlen:
                    best = h
            if not best:
                continue
            mapped += 1
            s, e = best.r_st, best.r_en
            if best.strand == 1:
                covp[s:e] += 1; fivep[s] += 1
            else:
                covm[s:e] += 1; fivem[e-1] += 1
    np.savez_compressed(f'{out_dir}/{run}.npz', covp=covp, covm=covm, fivep=fivep, fivem=fivem,
                        total=total, mapped=mapped, genome=genome_id, length=n)
    return f'{run}: {mapped}/{total} reads mapped to plastome ({100*mapped/max(total,1):.3f}%)'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genome', required=True); ap.add_argument('--reads', required=True)
    ap.add_argument('--out', default='cov'); ap.add_argument('--max-reads', type=int, default=0)
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args(); Path(a.out).mkdir(parents=True, exist_ok=True)
    rec = SeqIO.read(a.genome, 'genbank'); n = len(rec.seq)
    Path(FASTA).write_text(f'>{rec.id}\n{rec.seq}\n')
    runs = sorted(set(os.path.basename(f).split('_')[0] for f in glob.glob(f'{a.reads}/*_1.fastq.gz')))
    jobs = [(run, a.reads, a.out, a.max_reads, rec.id, n) for run in runs]
    with Pool(min(a.workers, len(jobs))) as pool:
        for msg in pool.imap_unordered(process_run, jobs):
            print(msg, flush=True)

if __name__ == '__main__':
    main()
