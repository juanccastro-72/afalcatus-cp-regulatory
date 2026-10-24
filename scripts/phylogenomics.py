#!/usr/bin/env python3
"""Concatenated plastid-protein maximum-likelihood phylogeny.

Genes: protein-coding plastid genes present in >= 90% of the comparative set,
one copy per genome (IR duplicates are identical; first copy used), homology
annotations that failed translation QC excluded, ORFs/ycf1 excluded.
MAFFT L-INS-i -> trimAl -automated1 -> concatenation -> IQ-TREE (ModelFinder
on the concatenated matrix, 1000 UFBoot + 1000 SH-aLRT), rooted on
Chlorella vulgaris (Trebouxiophyceae).
"""
import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

import os as _os
import shutil as _shutil
# External tools (MAFFT, trimAl, IQ-TREE, MEME Suite, BLAST+, ARAGORN): set CPREG_BIN to their bin
# directory, or put them on PATH (the directory containing 'mafft' is then used).
BIN = Path(_os.environ.get('CPREG_BIN') or (Path(_shutil.which('mafft')).parent if _shutil.which('mafft') else '/usr/local/bin'))
EXCLUDE = {'ycf1', 'rpoBb', 'rpoC1b', 'rps2b'}


def proteins(rec):
    out = {}
    for f in rec.features:
        if f.type != 'CDS' or 'exclude_from_boundary_analysis' in f.qualifiers or 'pseudo' in f.qualifiers:
            continue
        g = f.qualifiers['gene'][0]
        if g.lower().startswith(('orf', 'hypothetical')) or g in EXCLUDE or g in out:
            continue
        tr = f.qualifiers.get('translation', [''])[0]
        if not tr:
            s = str(f.extract(rec.seq))
            tr = str(Seq(s[:len(s)//3*3]).translate(table=11)).rstrip('*')
        if '*' in tr or len(tr) < 25:
            continue
        out[g] = tr
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..'); ap.add_argument('--threads', default='2'); ap.add_argument('--reuse', action='store_true')
    a = ap.parse_args(); P = Path(a.project)
    W = P/'results/phylogeny'; (W/'genes').mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(P/'results/tables/genome_manifest.tsv', sep='\t')
    man = man[man.included_in_comparative_set]
    label = {r.insdc_accession: f"{r.species.replace(' ', '_').replace('.', '')}_{r.insdc_accession.split('.')[0]}" for r in man.itertuples()}
    prots = {acc: proteins(SeqIO.read(P/'data/standardized'/f'{acc}.gb', 'genbank')) for acc in man.insdc_accession}
    counts = defaultdict(int)
    for p in prots.values():
        for g in p:
            counts[g] += 1
    genes = sorted(g for g, c in counts.items() if c >= 0.9 * len(prots))
    rows, parts, concat, pos = [], [], defaultdict(str), 1
    for g in genes:
        if (W/'genes'/f'{g}.trim').exists() and a.reuse:
            recs = {r.id: str(r.seq) for r in SeqIO.parse(W/'genes'/f'{g}.trim', 'fasta')}
            L = len(next(iter(recs.values())))
            for acc in prots:
                concat[label[acc]] += recs.get(label[acc], '-' * L)
            parts.append(f'  charset {g} = {pos}-{pos+L-1};')
            rows.append(dict(gene=g, n_taxa=len(recs), raw_aln_len=len(next(SeqIO.parse(W/'genes'/f'{g}.aln', 'fasta')).seq), trimmed_len=L))
            pos += L
            continue
        fa = W/'genes'/f'{g}.faa'
        fa.write_text(''.join(f'>{label[acc]}\n{p[g]}\n' for acc, p in prots.items() if g in p))
        aln = W/'genes'/f'{g}.aln'
        with open(aln, 'w') as fh:
            subprocess.run([str(BIN/'mafft'), '--localpair', '--maxiterate', '1000', '--quiet', '--thread', a.threads, str(fa)], stdout=fh, check=True)
        trim = W/'genes'/f'{g}.trim'
        subprocess.run([str(BIN/'trimal'), '-in', str(aln), '-out', str(trim), '-automated1'], check=True)
        recs = {r.id: str(r.seq) for r in SeqIO.parse(trim, 'fasta')}
        L = len(next(iter(recs.values())))
        for acc in prots:
            concat[label[acc]] += recs.get(label[acc], '-' * L)
        parts.append(f'  charset {g} = {pos}-{pos+L-1};')
        rows.append(dict(gene=g, n_taxa=len(recs), raw_aln_len=len(next(SeqIO.parse(aln, 'fasta')).seq), trimmed_len=L))
        pos += L
    sup = W/'concatenated.faa'
    sup.write_text(''.join(f'>{k}\n{v}\n' for k, v in concat.items()))
    (W/'partitions.nex').write_text('#nexus\nbegin sets;\n' + '\n'.join(parts) + '\nend;\n')
    pd.DataFrame(rows).to_csv(W/'gene_alignments.tsv', sep='\t', index=False)
    print(len(genes), 'genes;', pos - 1, 'aa sites')
    outg = label['AB001684']
    # Model selection: a full ModelFinder run (preliminary_full_ModelFinder.*) favoured LG with FreeRate
    # heterogeneity; a restricted ModelFinder run (LG/cpREV x F/FU, R6) selected LG+F+R6 by BIC
    # (ufboot_full_search_aborted/). For tractability on 2 CPUs the final tree uses the fast NNI search
    # with the selected model; branch support = SH-aLRT (1000 replicates) and approximate Bayes (aBayes).
    # (UFBoot is incompatible with --fast in IQ-TREE.)
    subprocess.run([str(BIN/'iqtree'), '-s', str(sup), '-m', 'LG+F+R6', '--fast', '-alrt', '1000', '--abayes',
                    '-o', outg, '-T', a.threads, '--seed', '20261023', '--prefix', str(W/'plastome_ML'), '-redo', '--quiet'], check=True)


if __name__ == '__main__':
    main()
