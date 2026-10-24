#!/usr/bin/env python3
"""Site-level QC of per-gene MEME footprints.

For every MEME site: genome, position relative to the first transcribed base,
whether the site could lie inside an upstream in-frame ORF extension (no
in-frame stop codon between the site and the annotated start codon; such a
site might be coding if the start codon were mis-annotated), and the best
TOMTOM match to the validated reference library.
"""
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

STOPS = ('TAA', 'TAG', 'TGA')


def extension_len(up300):
    """nt upstream of the start codon before the first in-frame stop (CDS frame)."""
    L = len(up300)
    for i in range(L - 3, -1, -3):
        if up300[i:i+3] in STOPS:
            return L - (i + 3)
    return L


def pwm(xml, motif_id):
    root = ET.parse(xml).getroot()
    for m in root.iter('motif'):
        if m.get('id') == motif_id:
            rows = [[float(v.text) for v in a.iter('value')] for a in m.iter('alphabet_array')]
            return pd.DataFrame(rows[1:] if len(rows) > int(m.get('width')) else rows, columns=list('ACGT'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; W = P/'results/motifs'
    A = pd.read_csv(T/'all_boundary_contexts.tsv.gz', sep='\t')
    A1 = A.drop_duplicates(['genome', 'gene']).set_index(['genome', 'gene'])
    D = pd.read_csv(T/'motif_discovery.tsv', sep='\t')
    tt = pd.read_csv(T/'motif_tomtom_vs_validated.tsv', sep='\t')
    rows, srows = [], []
    for r in D[D.analysis == 'per_gene_MEME'].itertuples():
        xml = W/f'meme_{r.group}'/'meme.xml'
        root = ET.parse(xml).getroot()
        seqs = {s.get('id'): (s.get('name'), int(s.get('length'))) for s in root.iter('sequence')}
        for m in root.iter('motif'):
            if m.get('id') != r.motif:
                continue
            w = int(m.get('width'))
            inext = []
            for s in m.iter('contributing_site'):
                g, L = seqs[s.get('sequence_id')]
                pos = int(s.get('position')) - L
                ctx = A1.loc[(g, r.group)] if (g, r.group) in A1.index else None
                ext = extension_len(ctx.upstream300) if (ctx is not None and ctx.feature_type == 'CDS') else 0
                coding_risk = (-pos) <= ext
                inext.append(coding_risk)
                site = ''.join(l.get('letter_id') for l in s.iter('letter_ref'))
                srows.append(dict(gene=r.group, motif=r.motif, consensus=r.consensus, genome=g, site_pos_rel_start=pos,
                                  site=site, inside_possible_ORF_extension=coding_risk))
        best = tt[(tt.group == r.group) & (tt.Query_ID == r.consensus)].sort_values('p-value').head(1)
        rows.append(dict(gene=r.group, motif=r.motif, consensus=r.consensus, width=r.width, sites=r.sites,
                         n_sequences=r.n_sequences, E_value=r.E_value, median_pos=r.median_site_pos_rel_start,
                         IQR_pos=r.IQR_site_pos, Afalcatus_pos=r.Afalcatus_site_pos,
                         site_fraction=r.sites / r.n_sequences,
                         frac_sites_possible_ORF_extension=float(np.mean(inext)) if inext else np.nan,
                         best_validated_match=best.Target_ID.iloc[0] if len(best) else '',
                         match_p=best['p-value'].iloc[0] if len(best) else np.nan,
                         match_q=best['q-value'].iloc[0] if len(best) else np.nan))
    R = pd.DataFrame(rows)
    R['positional_footprint'] = (R.E_value < 1e-10) & (R.site_fraction >= .6) & \
        R.IQR_pos.str.split('..', regex=False).apply(lambda x: abs(float(x[1]) - float(x[0])) <= 15 if isinstance(x, list) and len(x) == 2 else False)
    R.to_csv(T/'footprint_motifs_QC.tsv', sep='\t', index=False)
    pd.DataFrame(srows).to_csv(T/'footprint_sites.tsv.gz', sep='\t', index=False)
    print(R[R.positional_footprint][['gene', 'consensus', 'sites', 'n_sequences', 'E_value', 'median_pos', 'IQR_pos', 'Afalcatus_pos',
                                      'frac_sites_possible_ORF_extension', 'best_validated_match', 'match_q']].to_string())
    print(len(R), R.positional_footprint.sum())


if __name__ == '__main__':
    main()
