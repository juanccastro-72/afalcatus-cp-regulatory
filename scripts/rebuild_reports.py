#!/usr/bin/env python3
"""Rebuild Excel/figures from saved TSVs without repeating randomized controls."""
import argparse
import json
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from inventory_reporting import add_summary,create_workbook,create_figures

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project',default='.',help='Directory containing inputs/results/figures/logs.')
args=p.parse_args();root=Path(args.project)
data={}
for f in sorted((root/'results').glob('*.tsv')):
    if f.stem in {'Deliverable_QC'}:continue
    try:data[f.stem]=pd.read_csv(f,sep='\t')
    except pd.errors.EmptyDataError:data[f.stem]=pd.DataFrame()
record=SeqIO.read(root/'inputs'/'QB176266.1.gb','genbank')
add_summary(data)
data['Summary'].to_csv(root/'results'/'Summary.tsv',sep='\t',index=False)
meta=dict(zip(data['Metadata'].field,data['Metadata'].value))
temperature=json.loads(meta['parameters']).get('temperature',25.)
create_figures(data,record,root/'figures',temperature)
create_workbook(data,root,root/'results'/f'{record.id}_regulatory_inventory.xlsx')
print('Rebuilt workbook and figures from existing analysis results.')
