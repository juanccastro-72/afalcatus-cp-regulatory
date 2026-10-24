#!/usr/bin/env python3
"""Verify workbook data, genome coordinate paths, JPEG DPI and vector SVGs."""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
from PIL import Image
from openpyxl import load_workbook
from Bio import SeqIO

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project',default='.')
args=p.parse_args();root=Path(args.project)
rows=[]
def check(name,passed,details):
    rows.append({'check':name,'status':'PASS' if passed else 'FAIL','details':str(details)})
    if not passed:print('FAIL:',name,details)

results=root/'results';record=SeqIO.read(root/'inputs'/'QB176266.1.gb','genbank')
genome=str(record.seq).upper();comp=str.maketrans('ACGTN','TGCAN')
inv=pd.read_csv(results/'Inventory_associations.tsv',sep='\t')
hp=pd.read_csv(results/'Hairpins.tsv',sep='\t')
tr=pd.read_csv(results/'Translation_QC.tsv',sep='\t')
all_paths=True
for _,r in inv.iterrows():
    positions=[int(x)-1 for x in r.position_path_1b.split(',')]
    seq=''.join(genome[x] if r.strand==1 else genome[x].translate(comp) for x in positions)
    all_paths &= seq==r.sequence_5to3 and len(seq)==r.length_bp
check('all_candidate_coordinate_round_trips',all_paths,f'{len(inv)} candidate/context associations')
check('all_hairpins_antiparallel',all(r.left_arm==r.right_arm.translate(comp)[::-1] for _,r in hp.iterrows()),f'{len(hp)} hairpin alternatives')
check('MFE_not_above_intended_energy',bool((hp.MFE_kcal_mol<=hp.intended_energy_kcal_mol+.02).all()),'0.02 kcal/mol floating-point tolerance')
check('all_CDS_match_annotated_proteins',bool(tr.matches_annotation.all()),f'{len(tr)} CDS')
check('all_CDS_no_internal_stops',bool((tr.internal_stop_count==0).all()),'Annotated translation table and codon_start used')

wb=load_workbook(results/f'{record.id}_regulatory_inventory.xlsx',read_only=True,data_only=True)
required=['README','Summary','Gene_summary','Promoters','SD_like','Hairpins','Inventory_associations','Unique_sites','Null_tests','References','Code_audit','Unit_tests']
check('workbook_required_sheets',all(x in wb.sheetnames for x in required),','.join(wb.sheetnames))
for sheet,file in [('Gene_summary','Gene_summary'),('Promoters','Promoters'),('SD_like','SD_like'),('Hairpins','Hairpins'),('Inventory_associations','Inventory_associations'),('Unique_sites','Unique_sites'),('Null_tests','Null_tests')]:
    table=pd.read_csv(results/f'{file}.tsv',sep='\t')
    check(f'workbook_{sheet}_rows',wb[sheet].max_row-1==len(table),f'Excel={wb[sheet].max_row-1}; TSV={len(table)}')
wb.close()
manifest=pd.read_csv(results/'Figure_manifest.tsv',sep='\t')
for _,r in manifest.iterrows():
    with Image.open(root/r.JPEG) as im:
        dpi=im.info.get('dpi',(0,0))
        check(r.figure+'_JPEG_600dpi',all(abs(float(x)-600)<.5 for x in dpi),f'{im.size[0]}x{im.size[1]} pixels; dpi={dpi}')
    svg=ET.parse(root/r.SVG).getroot()
    tags=[x.tag.split('}')[-1] for x in svg.iter()]
    check(r.figure+'_SVG_vector_only','image' not in tags and 'text' not in tags and 'script' not in tags,
          f'SVG {svg.get("version")}; images={tags.count("image")}; live text={tags.count("text")}; paths={tags.count("path")}')
    check(r.figure+'_SVG_has_viewBox','viewBox' in svg.attrib,svg.attrib.get('viewBox','missing'))
for filename,bed in [('regulatory_candidates.bed',True),('regulatory_candidates.gff3',False)]:
    good=True;count=0
    for line in (results/filename).read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t');a,b=(int(f[1]),int(f[2])) if bed else (int(f[3])-1,int(f[4]))
        good &= 0<=a<b<=len(record);count+=1
    check(filename+'_bounds',good,f'{count} intervals; length={len(record)}')
tests=pd.read_csv(root/'logs'/'unit_test_results.tsv',sep='\t')
check('regression_tests_pass',bool((tests.status=='PASS').all()),f'{len(tests)} tests')
rows.append({'check':'CorelDRAW_application_import','status':'NOT_TESTED','details':'Proprietary application unavailable; SVG 1.1 XML/vector/no-font/no-raster checks passed where indicated.'})
qc=pd.DataFrame(rows)
qc.to_csv(root/'logs'/'deliverable_qc.tsv',sep='\t',index=False)
qc.to_csv(results/'Deliverable_QC.tsv',sep='\t',index=False)
print(qc.groupby('status').size().to_string())
if (qc.status=='FAIL').any():raise SystemExit(1)
