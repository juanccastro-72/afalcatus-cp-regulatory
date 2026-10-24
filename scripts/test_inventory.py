#!/usr/bin/env python3
"""Focused regression tests; run from project root: python scripts/test_inventory.py."""
from collections import Counter
import random
import unittest
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
sys.path.insert(0,str(Path(__file__).resolve().parent))
import chloroplast_regulatory_inventory_improved as m


class InventoryTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(m.rc('AGTC'),'GACT')

    def test_true_antiparallel_stem(self):
        self.assertTrue(any(h[:4]==(0,11,4,3) for h in m.hairpins('AGTCAAAGACT')))

    def test_parallel_complement_is_not_stem(self):
        self.assertFalse(any(h[:4]==(0,11,4,3) for h in m.hairpins('AGTCAAATCAG')))

    def test_circular_plus_window(self):
        s,p=m.windows('ACGTACGT',1,1,range(-3,0))
        self.assertEqual(p,[6,7,0]);self.assertEqual(s,'GTA')

    def test_circular_minus_window(self):
        s,p=m.windows('ACGTACGT',6,-1,range(-3,0))
        self.assertEqual(p,[1,0,7]);self.assertEqual(s,'GTA')

    def test_linear_clips_not_wraps(self):
        s,p=m.windows('ACGTACGT',1,1,range(-3,0),False)
        self.assertEqual(s,'A');self.assertEqual(p,[0])

    def test_wrapped_coordinate_segments(self):
        x=m.path_fields([6,7,0,1],8)
        self.assertTrue(x['wraps_origin']);self.assertEqual(x['genomic_segments_1b'],'7..8;1..2')

    def test_spaced_promoter_pair(self):
        h=m.promoter_pairs('TTGACA'+'C'*17+'TATAAT')
        self.assertTrue(any(x[0]==0 and x[3]==17 and x[4:6]==(0,0) for x in h))

    def test_wrong_spacer_rejected(self):
        self.assertFalse(m.promoter_pairs('TTGACA'+'C'*10+'TATAAT'))

    def test_sd_spacing_after_motif(self):
        seq='C'*8+'AGGAGG'+'C'*7
        h=m.sd_hits(seq);self.assertEqual(len(h),1)
        self.assertEqual(h[0][3],7);self.assertEqual(h[0][2],'AGGAGG')

    def test_sd_too_close_rejected(self):
        self.assertFalse(m.sd_hits('CCCCGAGGCC'))

    def test_overlapping_motif_matches(self):
        self.assertEqual([x.start() for x in m.OCT.finditer('TATAATATATAATAT')],[0,7])

    def test_shuffle_preserves_dinucleotides_and_ends(self):
        seq='AATATGGCCAATGCGTAACCGGTTATATA';rng=random.Random(7)
        for _ in range(50):
            shuffled=m.dinuc_shuffle(seq,rng)
            self.assertEqual(Counter(zip(seq,seq[1:])),Counter(zip(shuffled,shuffled[1:])))
            self.assertEqual((seq[0],seq[-1]),(shuffled[0],shuffled[-1]))

    def test_bh_known_values(self):
        np.testing.assert_allclose(m.bh_adjust([.01,.04,.03,.002]),[.02,.04,.04,.008])

    def test_empty_annotation_safe(self):
        # All empty exporters are tested separately through a synthetic run.
        self.assertEqual(m.hairpins('AAAA'),[]);self.assertEqual(m.sd_hits('AAAA'),[])

    def test_zero_hit_categories_have_schemas(self):
        rec=SeqRecord(Seq('A'*500),id='zero')
        rec.annotations={'topology':'circular'}
        rec.features=[SeqFeature(SimpleLocation(200,250,1),type='rRNA',qualifiers={'gene':['x']})]
        _,p,a,_,_,occ,mask=m.make_annotations(rec)
        args=SimpleNamespace(seed=7,temperature=25,upstream=50,downstream=30,null_reps=2,fold_null_reps=0)
        d=m.analyze(rec,a,occ,mask,args)
        for key in ['Promoters','SD_like','Hairpins']:
            self.assertTrue(d[key].empty);self.assertIn('element_type',d[key].columns)
        self.assertEqual(m.gene_summary(p,a,d).PEP_pairs.iloc[0],0)

    def test_mixed_strand_exon_contexts(self):
        rec=SeqRecord(Seq('A'*200),id='test')
        rec.features=[SeqFeature(CompoundLocation([SimpleLocation(100,112,1),SimpleLocation(10,31,-1)]),type='CDS',qualifiers={'gene':['x'],'trans_splicing':[''],'translation':['KKKKKKKKKKK']})]
        _,_,a,_,_,_,_=m.make_annotations(rec)
        self.assertEqual(len(a),2)
        self.assertEqual(a.first_transcribed_base_1b.tolist(),[101,31])
        self.assertEqual(a.sd_eligible.tolist(),[True,False])

    def test_minus_cis_spliced_boundaries(self):
        rec=SeqRecord(Seq('A'*100),id='test')
        rec.features=[SeqFeature(CompoundLocation([SimpleLocation(70,90,-1),SimpleLocation(10,30,-1)]),type='rRNA',qualifiers={'gene':['r']})]
        _,_,a,_,_,_,_=m.make_annotations(rec)
        self.assertEqual(a.first_transcribed_base_1b.tolist(),[90])
        self.assertEqual(a.last_transcribed_base_1b.tolist(),[11])

    def test_origin_neighborhood_and_minus_order(self):
        a=pd.DataFrame([{'anchor_id':'a','product_id':'p1','gene':'left','strand':-1,'context_start_1b':5,'context_end_1b':15},
                        {'anchor_id':'b','product_id':'p2','gene':'right','strand':-1,'context_start_1b':80,'context_end_1b':95}])
        nb,e=m.neighborhoods(a,100,gapmax=20,circular=True)
        self.assertEqual(len(nb),1);self.assertTrue(nb.iloc[0].crosses_origin)
        self.assertEqual(nb.iloc[0].genes_transcription_order,'left;right')

    def test_actual_genome_product_contexts(self):
        path=Path(__file__).resolve().parents[1]/'inputs'/'QB176266.1.gb'
        if not path.exists():self.skipTest('Supplied genome not present')
        rec=SeqIO.read(path,'genbank');_,p,a,_,tr,_,_=m.make_annotations(rec)
        self.assertEqual(len(p),110);self.assertEqual(len(a),113)
        self.assertEqual(set(p[p.mixed_strand].gene),{'psaA','rpl32'})
        self.assertEqual(int(a.sd_eligible.sum()),74)


if __name__=='__main__':
    class RecordingResult(unittest.TextTestResult):
        rows=[]
        def addSuccess(self,test):
            super().addSuccess(test);self.rows.append({'test':str(test),'status':'PASS','details':''})
        def addFailure(self,test,err):
            super().addFailure(test,err);self.rows.append({'test':str(test),'status':'FAIL','details':str(err[1])})
        def addError(self,test,err):
            super().addError(test,err);self.rows.append({'test':str(test),'status':'ERROR','details':str(err[1])})
        def addSkip(self,test,reason):
            super().addSkip(test,reason);self.rows.append({'test':str(test),'status':'SKIP','details':reason})
    result=unittest.TextTestRunner(verbosity=2,resultclass=RecordingResult).run(unittest.defaultTestLoader.loadTestsFromTestCase(InventoryTests))
    path=Path(__file__).resolve().parents[1]/'logs';path.mkdir(exist_ok=True)
    pd.DataFrame(result.rows).to_csv(path/'unit_test_results.tsv',sep='\t',index=False)
    sys.exit(not result.wasSuccessful())
