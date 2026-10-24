#!/usr/bin/env python3
"""Standardise all plastome records and audit them against Table S1.

Outputs data/standardized/<acc>.gb and results/tables/genome_manifest.tsv.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from cpreg_genomes import load_record, standardise, homology_annotate, n_products

HOMOLOGY_REFS = {  # unannotated WGS scaffolds -> closest annotated relatives
    'CM002678': ['QB176266.1', 'OM329943', 'KT199250'],
    'CM031383': ['MF423120', 'DQ396875', 'OR502664'],
}
DUPLICATE_OF = {'CM007919': 'KX756229'}
SUPPLEMENT = {'OM329943': ['QB176266.1', 'KT199250', 'CM002678']}  # incomplete submitter annotation


def flag_failed(rec):
    for f in rec.features:
        if f.type == 'CDS' and f.qualifiers.get('inference', [''])[0] == 'homology_transfer_v3':
            note = f.qualifiers['note'][0]
            if 'pseudo' in f.qualifiers or 'internal_stops=0' not in note:
                f.qualifiers['pseudo'] = ['']
                f.qualifiers['exclude_from_boundary_analysis'] = ['homology ORF failed translation QC (intron/split gene?)']
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project', default='.')
    a = ap.parse_args()
    P = Path(a.project)
    man = pd.read_csv(P/'data/genome_manifest_input.tsv', sep='\t')
    s1 = pd.read_excel(P/'data/Table_S1.xlsx', sheet_name=0)
    s1.columns = [c.strip() for c in s1.columns]
    out = P/'data/standardized'; out.mkdir(parents=True, exist_ok=True)
    (P/'results/tables').mkdir(parents=True, exist_ok=True)
    recs, rows = {}, []
    for r in man.itertuples():
        acc = r.insdc_accession
        f = P/'data/genbank'/(f'{acc}.gb' if acc == 'QB176266.1' else f'{acc}.embl')
        raw = f.read_bytes()
        rec = standardise(load_record(f), acc)
        recs[acc] = rec
        rows.append(dict(species=r.species, family=r.family, clade=r.order_clade, role=r.role,
                         tableS1_accession=r.tableS1_accession, insdc_accession=acc, file=f.name,
                         file_sha256=hashlib.sha256(raw).hexdigest(),
                         seq_sha256=hashlib.sha256(str(rec.seq).upper().encode()).hexdigest(),
                         length_bp=len(rec.seq), topology=rec.annotations.get('topology'),
                         record_date=rec.annotations.get('date'),
                         organism=rec.annotations.get('organism'),
                         n_CDS_original=sum(1 for x in rec.features if x.type == 'CDS')))
    for acc, refs in HOMOLOGY_REFS.items():
        recs[acc] = flag_failed(homology_annotate(recs[acc], [recs[x] for x in refs], f'refs={",".join(refs)}'))
    for acc, refs in SUPPLEMENT.items():
        extra = flag_failed(homology_annotate(recs[acc], [recs[x] for x in refs], f'supplement refs={",".join(refs)}'))
        have = {f.qualifiers['gene'][0] for f in recs[acc].features if f.type in ('CDS', 'rRNA')}
        add = [f for f in extra.features if f.type in ('CDS', 'rRNA') and f.qualifiers['gene'][0] not in have]
        recs[acc].features = sorted(recs[acc].features + add, key=lambda f: (f.type != 'source', int(f.location.start)))
        print(acc, 'supplemented', len(add), sorted(f.qualifiers['gene'][0] for f in add))
    for row in rows:
        acc = row['insdc_accession']; rec = recs[acc]
        s = str(rec.seq).upper()
        row['GC_percent'] = round(100*(s.count('G')+s.count('C'))/len(s), 3)
        row['ambiguous_bases'] = len(s) - sum(s.count(x) for x in 'ACGT')
        row['annotation_source'] = ('homology_transfer_v3' if acc in HOMOLOGY_REFS else
                                    'INSDC submitter + homology supplement' if acc in SUPPLEMENT else 'INSDC submitter')
        row['n_CDS_failed_QC_excluded'] = sum(1 for x in rec.features if 'exclude_from_boundary_analysis' in x.qualifiers)
        row['n_CDS_used'] = sum(1 for x in rec.features if x.type == 'CDS')
        row['n_rRNA_used'] = sum(1 for x in rec.features if x.type == 'rRNA')
        row['n_tRNA_used'] = sum(1 for x in rec.features if x.type == 'tRNA')
        row['included_in_comparative_set'] = acc not in DUPLICATE_OF
        row['exclusion_reason'] = (f'sequence identical to {DUPLICATE_OF[acc]} (same assembly, unannotated CON record)'
                                   if acc in DUPLICATE_OF else '')
        # Table S1 audit
        t = s1[s1['GenBank Accession Code'].astype(str).str.strip() == row['tableS1_accession']]
        if len(t):
            t = t.iloc[0]
            row['tableS1_taxid'] = t['Taxonomy']; row['tableS1_size_Mb'] = t['Size (Mb)']; row['tableS1_GC'] = t['GC (%)']
            row['size_diff_vs_tableS1_bp'] = int(round(len(s) - 1e6*float(t['Size (Mb)'])))
            row['GC_diff_vs_tableS1'] = round(row['GC_percent'] - float(t['GC (%)']), 2)
            row['record_taxid'] = next((x.split(':')[1] for f in rec.features if f.type == 'source'
                                        for x in f.qualifiers.get('db_xref', []) if x.startswith('taxon:')), '')
            row['taxid_matches_tableS1'] = str(row['record_taxid']) == str(t['Taxonomy'])
        if row['included_in_comparative_set']:
            SeqIO.write(rec, out/f'{acc}.gb', 'genbank')
    df = pd.DataFrame(rows)
    df.to_csv(P/'results/tables/genome_manifest.tsv', sep='\t', index=False)
    print(df[['species', 'insdc_accession', 'length_bp', 'GC_percent', 'n_CDS_used', 'n_rRNA_used',
              'n_tRNA_used', 'size_diff_vs_tableS1_bp', 'GC_diff_vs_tableS1', 'taxid_matches_tableS1']].to_string())


if __name__ == '__main__':
    main()
