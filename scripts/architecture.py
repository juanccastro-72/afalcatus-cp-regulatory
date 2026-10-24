#!/usr/bin/env python3
"""Plastome architecture, gene content and gene-order (co-orientation) conservation."""
import argparse
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

import os as _os
import shutil as _shutil
# External tools (MAFFT, trimAl, IQ-TREE, MEME Suite, BLAST+, ARAGORN): set CPREG_BIN to their bin
# directory, or put them on PATH (the directory containing 'mafft' is then used).
BIN = Path(_os.environ.get('CPREG_BIN') or (Path(_shutil.which('mafft')).parent if _shutil.which('mafft') else '/usr/local/bin'))


def inverted_repeat(seq):
    with tempfile.TemporaryDirectory() as td:
        f = Path(td)/'g.fa'; f.write_text(f'>g\n{seq}\n')
        out = subprocess.run([str(BIN/'blastn'), '-query', str(f), '-subject', str(f), '-strand', 'minus',
                              '-evalue', '1e-50', '-outfmt', '6 qstart qend sstart send length pident', '-dust', 'no'],
                             stdout=subprocess.PIPE, text=True, check=True).stdout
    best = None
    for line in out.splitlines():
        qs, qe, ss, se, L, pid = line.split('\t')
        qs, qe, ss, se, L, pid = int(qs), int(qe), int(ss), int(se), int(L), float(pid)
        if qs < se and L >= 1000 and pid >= 95 and (best is None or L > best[4]):
            best = (qs, qe, ss, se, L, pid)
    return best


def gene_list(rec):
    out = []
    for f in rec.features:
        if f.type not in ('CDS', 'rRNA') or 'exclude_from_boundary_analysis' in f.qualifiers:
            continue
        g = f.qualifiers['gene'][0]
        if g.lower().startswith(('orf', 'hypothetical')) or ' ' in g:
            continue
        parts = list(f.location.parts)
        for p in parts if len({x.strand for x in parts}) > 1 else [f.location]:
            out.append((int(p.start), int(p.end), p.strand, g))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t'); man = man[man.included_in_comparative_set]
    arch, content, orders = [], defaultdict(dict), {}
    for r in man.itertuples():
        rec = SeqIO.read(P/'data/standardized'/f'{r.insdc_accession}.gb', 'genbank')
        s = str(rec.seq).upper(); n = len(s)
        cov = np.zeros(n, bool); introns = 0; trans = 0; genes = set()
        for f in rec.features:
            if f.type in ('CDS', 'rRNA', 'tRNA'):
                for p in f.location.parts:
                    cov[int(p.start):int(p.end)] = True
                if f.type == 'CDS':
                    parts = list(f.location.parts)
                    if len({p.strand for p in parts}) > 1 or 'trans_splicing' in f.qualifiers:
                        trans += 1
                    elif len(parts) > 1:
                        introns += len(parts) - 1
                g = f.qualifiers['gene'][0]
                if f.type in ('CDS', 'rRNA') and not g.lower().startswith(('orf', 'hypothetical')) and ' ' not in g:
                    content[r.insdc_accession][g] = content[r.insdc_accession].get(g, 0) + 1
                    genes.add(g)
                if f.type == 'tRNA':
                    genes.add(g)
        ir = inverted_repeat(s)
        irl = ir[4] if ir else 0
        nc = s_nc = ''.join(c for c, x in zip(s, cov) if not x)
        arch.append(dict(genome=r.insdc_accession, species=r.species, family=r.family, clade=r.clade,
                         length_bp=n, GC=100*(s.count('G')+s.count('C'))/n,
                         IR_length_bp=irl, IR_identity=ir[5] if ir else np.nan,
                         IR_coords=f'{ir[0]}-{ir[1]} / {ir[3]}-{ir[2]}' if ir else 'not detected',
                         coding_fraction=cov.mean(), intergenic_fraction=1-cov.mean(),
                         intergenic_GC=100*(nc.count('G')+nc.count('C'))/max(len(nc), 1),
                         n_distinct_genes=len(genes), n_cis_introns_CDS=introns, n_trans_spliced_CDS=trans))
        orders[r.insdc_accession] = gene_list(rec)
    A = pd.DataFrame(arch); A.to_csv(T/'plastome_architecture.tsv', sep='\t', index=False)
    C = pd.DataFrame(content).fillna(0).astype(int).T
    C.index.name = 'genome'; C.to_csv(T/'gene_content_matrix.tsv', sep='\t')
    # co-oriented adjacency in transcription order
    def pairs(gl, n):
        out = {}
        for i in range(len(gl)):
            x, y = gl[i], gl[(i+1) % len(gl)]
            if x[2] == y[2] and x[3] != y[3]:
                gap = (y[0] - x[1]) % n if x[2] == 1 else (y[0] - x[1]) % n
                key = (x[3], y[3]) if x[2] == 1 else (y[3], x[3])
                out[key] = gap
        return out
    lens = dict(zip(A.genome, A.length_bp))
    P_all = {g: pairs(o, lens[g]) for g, o in orders.items()}
    focal = 'QB176266.1'
    rows = []
    for key, gap in P_all[focal].items():
        cons = [g for g in P_all if g != focal and key in P_all[g]]
        rows.append(dict(upstream_gene=key[0], downstream_gene=key[1], focal_gap_bp=gap,
                         n_genomes_conserved=len(cons), frac_conserved=len(cons)/(len(P_all)-1),
                         conserved_in=';'.join(cons)))
    G = pd.DataFrame(rows).sort_values('frac_conserved', ascending=False)
    G.to_csv(T/'gene_pair_coorientation_conservation.tsv', sep='\t', index=False)
    # conserved clusters in A. falcatus (chains of pairs conserved in >= 50%)
    ok = {(r.upstream_gene, r.downstream_gene) for r in G.itertuples() if r.frac_conserved >= .5}
    nxt = {u: d for u, d in ok}; prv = {d: u for u, d in ok}
    clusters = []
    for u in nxt:
        if u in prv:
            continue
        chain = [u]
        while chain[-1] in nxt and nxt[chain[-1]] not in chain:
            chain.append(nxt[chain[-1]])
        clusters.append(chain)
    pd.DataFrame([dict(cluster_id=f'CL{i+1:02d}', n_genes=len(c), genes='-'.join(c)) for i, c in
                  enumerate(sorted(clusters, key=len, reverse=True))]).to_csv(T/'conserved_cooriented_clusters.tsv', sep='\t', index=False)
    print(A[['species', 'length_bp', 'IR_length_bp', 'intergenic_fraction', 'n_distinct_genes']].to_string())
    print(G.head(20).to_string())


if __name__ == '__main__':
    main()
