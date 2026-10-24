#!/usr/bin/env python3
"""Reproduce original logic offline without repeating the original NCBI requests."""
import importlib.util
import json
from pathlib import Path
from Bio import SeqIO

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('original',root/'inputs'/'original_chloroplast_regulatory_inventory.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
r=SeqIO.read(root/'inputs'/'QB176266.1.gb','genbank')
f=old.extract_features(r);t=old.identify_transcription_units(f)
p=old.predict_promoters(r,t);s=old.predict_shine_dalgarno(r,t);h=old.predict_terminators(r,t);o=old.predict_operons(t)
summary={'original_features_extracted':len(f),'original_candidate_TUs':len(t),
    'original_promoter_motif_associations':len(p),'original_SD_associations':len(s),
    'original_hairpin_associations':len(h),'original_operon_groups_including_singletons':len(o),
    'original_multigene_groups':int((o.n_genes>1).sum()),
    'original_singleton_groups_mislabeled_operons':int((o.n_genes==1).sum()),
    'original_total_associations':len(p)+len(s)+len(h),
    'original_promoter_motif_counts':p.motif_name.value_counts().to_json()}
(root/'inputs'/'baseline_original_summary.json').write_text(json.dumps(summary,indent=2))
for name,df in [('promoters',p),('SD',s),('hairpins',h),('operons',o)]:
    df.to_csv(root/'inputs'/f'original_baseline_{name}.tsv',sep='\t',index=False)
print(json.dumps(summary,indent=2))
