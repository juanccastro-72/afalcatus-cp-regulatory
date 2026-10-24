#!/usr/bin/env python3
"""Regression tests for the comparative (v3) modules. Run: python test_comparative.py"""
import random
from collections import Counter

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import chloroplast_regulatory_inventory_improved as v2
import cpreg_genomes as cg
import cpreg_scan as cs
from motif_qc import extension_len

results = []


def check(name, cond):
    results.append((name, bool(cond)))
    print(('PASS ' if cond else 'FAIL ') + name)


# 1 gene-name normalisation
f = SeqFeature(FeatureLocation(0, 10, 1), type='CDS', qualifiers={'gene': ['rplB']})
check('rplB -> rpl2', cg.normalise_name(f)[1] == 'rpl2')
f = SeqFeature(FeatureLocation(0, 10, 1), type='rRNA', qualifiers={'gene': ['rrs']})
check('rrs -> rrn16', cg.normalise_name(f)[1] == 'rrn16')
f = SeqFeature(FeatureLocation(0, 10, 1), type='rRNA', qualifiers={'product': ['23S ribosomal RNA']})
check('23S product -> rrn23', cg.normalise_name(f)[1] == 'rrn23')
f = SeqFeature(FeatureLocation(0, 10, 1), type='CDS', qualifiers={'gene': ['cbbL']})
check('cbbL -> rbcL', cg.normalise_name(f)[1] == 'rbcL')

# 2 dinucleotide shuffle preserves dinucleotides and ends
rng = random.Random(1)
s = ''.join(rng.choice('AAATTTGC') for _ in range(300))
t = v2.dinuc_shuffle(s, rng)
check('shuffle preserves dinucleotides', Counter(zip(s, s[1:])) == Counter(zip(t, t[1:])) and s[0] == t[0] and s[-1] == t[-1])

# 3 octamer / PEP rules
check('octamer overlapping count', cs.count(cs.OCT, 'TATAATATAATAT') == 2)
seq = 'TTGACA' + 'A' * 17 + 'TATAAT'
check('strict PEP pair 17-nt spacer', len(cs.strong_pairs(seq)) == 1)
check('PEP pair 19-nt spacer rejected by strict rule', len(cs.strong_pairs('TTGACA' + 'A' * 19 + 'TATAAT')) == 0)

# 4 SD spacer measured after motif end
short = 'A' * 20 + 'GGAGG' + 'A' * 10
hits = v2.sd_hits(short)
check('SD spacer = bases after motif (longest word AGGAGG kept)', any(h[2] == 'AGGAGG' and h[3] == 10 for h in hits))

# 5 ORF-extension length
check('extension stops at in-frame stop', extension_len('AAAAAA' + 'TAA' + 'AAAAAA') == 6)
check('no in-frame stop -> full length', extension_len('AAAAAAAAA') == 9)

# 6 scan on synthetic circular record: non-coding masking and wrap-around
g = ('A' * 50 + 'TATAATAT' + 'A' * 42) + 'ATG' + 'GCT' * 30 + 'TAA' + 'GCGCGCAAAAGCGCGC' + 'T' * 80
rec = SeqRecord(Seq(g), id='syn', annotations={'topology': 'circular', 'molecule_type': 'DNA'},
                features=[SeqFeature(FeatureLocation(100, 196, 1), type='CDS',
                                     qualifiers={'gene': ['x'], 'codon_start': ['1'], 'transl_table': ['11'],
                                                 'translation': [str(Seq(g[100:196]).translate(table=11)).rstrip('*')]})])
res = cs.scan_genome(rec, 'syn', nrep=10, nrep_fold=5, nrep_start=5)
a = res['anchors'].iloc[0]
check('octamer found at -50 relative to start', a.TATAATAT_pos == '-50')
check('upstream window wraps circular origin (300 nt)', len(a.upstream300) == 300)
check('non-coding stretch stops at the CDS after wrapping (196 nt)', a.upstream_noncoding_len == 196)
check('start codon recorded', a.start_codon == 'ATG')

# 7 antiparallel hairpin
hp = v2.hairpins('AAAAGCGCAAATTTGCGCAAAA')
check('antiparallel stem detected', any(h[2] >= 4 for h in hp))

ok = all(x for _, x in results)
print(f'{sum(x for _, x in results)}/{len(results)} tests passed')
raise SystemExit(0 if ok else 1)
