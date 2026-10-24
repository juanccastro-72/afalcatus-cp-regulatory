#!/usr/bin/env python3
"""Genome loading, gene-name normalisation and homology-based re-annotation.

Part of the comparative chloroplast regulatory-element workflow (v3).

* INSDC records are read from GenBank (.gb) or ENA EMBL (.embl) flat files.
  ENA "CON" records carry a CO line before the expanded SQ block, which
  Biopython rejects; CO lines are removed before parsing (sequence unchanged).
* Gene names are normalised to plastid nomenclature (rplB -> rpl2, cbbL -> rbcL,
  rrs/16S -> rrn16, ...). The original name is always kept.
* Records without usable annotation (WGS chloroplast scaffolds) are annotated by
  homology transfer: tblastn (proteins) / blastn (rRNAs) using annotated
  Sphaeropleales references, with start/stop extension in frame, plus ARAGORN
  for tRNAs. These features are flagged `inference=homology_transfer_v3` and
  are QC-checked by translation (no internal stops; complete ORF).
"""
from __future__ import annotations

import io
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import os as _os
import shutil as _shutil
# External tools (MAFFT, trimAl, IQ-TREE, MEME Suite, BLAST+, ARAGORN): set CPREG_BIN to their bin
# directory, or put them on PATH (the directory containing 'mafft' is then used).
BIN = Path(_os.environ.get('CPREG_BIN') or (Path(_shutil.which('mafft')).parent if _shutil.which('mafft') else '/usr/local/bin'))

SYN = {
    'rplb': 'rpl2', 'rple': 'rpl5', 'rpll': 'rpl12', 'rpln': 'rpl14', 'rplp': 'rpl16',
    'rplt': 'rpl20', 'rplw': 'rpl23', 'rpmj': 'rpl36', 'rpsb': 'rps2', 'rpsc': 'rps3',
    'rpsd': 'rps4', 'rpsg': 'rps7', 'rpsh': 'rps8', 'rpsi': 'rps9', 'rpsk': 'rps11',
    'rpsl': 'rps12', 'rpsn': 'rps14', 'rpsr': 'rps18', 'rpss': 'rps19',
    'cbbl': 'rbcL', 'tuf': 'tufA', 'pafi': 'ycf3', 'pafii': 'ycf4', 'pbf1': 'psbN',
    'psb30': 'ycf12', 'clpp1': 'clpP', 'psbd1': 'psbD', 'rpoc': 'rpoC1',
    'rpob-1': 'rpoBa', 'rpob-2': 'rpoBb', 'rpoc1a': 'rpoC1', 'rpoc1b': 'rpoC1b',
    'rps2a': 'rps2', 'rps2b': 'rps2b',
}
PRODUCT_SYN = [
    (r'30S ribosomal protein S3', 'rps3'), (r'alpha subunit of RNA polymerase', 'rpoA'),
    (r'^beta subunit of RNA polymerase', 'rpoB'), (r"subunit beta'", 'rpoC1'),
    (r'hotochlorophyllide reductase subunit N', 'chlN'),
    (r'protochlorophyllide reductase subunit B', 'chlB'),
    (r'chloroplast protein RF1', 'ycf1'), (r'clp protease', 'clpP'),
]
RRNA = [(r'(^|[^0-9])(16S|rrn16|rrs|rrnS|rns|small subunit)', 'rrn16'),
        (r'(^|[^0-9])(23S|rrn23|rrl|rrnL|large subunit)', 'rrn23'),
        (r'(^|[^0-9])(5S|rrn5|rrf)', 'rrn5'),
        (r'rrn3', 'rrn3'), (r'rrn7', 'rrn7')]

CANONICAL = None  # filled lazily: set of plastid gene names with canonical case


def load_record(path):
    path = Path(path)
    if path.suffix in ('.gb', '.gbk', '.genbank'):
        return SeqIO.read(path, 'genbank')
    text = ''.join(l for l in open(path) if not l.startswith('CO '))
    return SeqIO.read(io.StringIO(text), 'embl')


def normalise_name(feature):
    q = feature.qualifiers
    raw = (q.get('gene') or q.get('locus_tag') or q.get('product') or ['?'])[0].strip()
    product = q.get('product', [''])[0]
    if feature.type == 'rRNA' or raw.lower() in ('rrs', 'rrl', 'rrf'):
        for pat, name in RRNA:
            if re.search(pat, raw, re.I) or re.search(pat, product, re.I):
                return raw, name
        return raw, raw
    if feature.type == 'tRNA':
        return raw, raw
    key = raw.lower()
    if key in SYN:
        return raw, SYN[key]
    if raw == '?' or ' ' in raw:
        for pat, name in PRODUCT_SYN:
            if re.search(pat, raw) or re.search(pat, product):
                return raw, name
    if re.fullmatch(r'(accd|acc[d])', key):
        return raw, 'accD'
    return raw, raw


def standardise(record, accession):
    """Return a copy with `gene_norm` and `gene_raw` qualifiers on product features."""
    for f in record.features:
        if f.type in ('CDS', 'rRNA', 'tRNA'):
            raw, norm = normalise_name(f)
            # rRNA mis-typed as CDS in some submissions
            if f.type == 'CDS' and norm in ('rrn16', 'rrn23', 'rrn5'):
                f.type = 'rRNA'
            f.qualifiers['gene_raw'] = [raw]
            f.qualifiers['gene'] = [norm]
            f.qualifiers.setdefault('note', []).append(f'v3 normalised gene name from {raw!r}')
    record.annotations.setdefault('molecule_type', 'DNA')
    record.id = accession
    return record


def n_products(record):
    return sum(1 for f in record.features if f.type == 'CDS' and f.qualifiers.get('gene', [''])[0] not in ('rrn16', 'rrn23'))


# ----------------------------------------------------------------------------
# Homology transfer for unannotated plastome scaffolds
# ----------------------------------------------------------------------------
STARTS = ('ATG', 'GTG', 'TTG')
STOPS = ('TAA', 'TAG', 'TGA')


def _run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def reference_sets(ref_records):
    prot, rrna = {}, {}
    for rec in ref_records:
        for f in rec.features:
            g = f.qualifiers.get('gene', [''])[0]
            if f.type == 'CDS' and g and not g.lower().startswith('orf') and g not in prot:
                tr = f.qualifiers.get('translation', [''])[0]
                if len(tr) >= 25:
                    prot[g] = tr
            if f.type == 'rRNA' and g in ('rrn16', 'rrn23', 'rrn5') and g not in rrna:
                rrna[g] = str(f.extract(rec.seq))
    return prot, rrna


def _orf_from_hit(genome2, strand, s, e, qstart, qlen, n):
    """Refine a tblastn HSP into a complete ORF on a doubled circular sequence.

    s,e: 0-based half-open coordinates of the HSP on the + strand of genome2.
    Returns (start0, end0) on genome2 + strand (half-open) or None.
    """
    seq = genome2 if strand == 1 else str(Seq(genome2).reverse_complement())
    L = len(seq)
    if strand == -1:
        s, e = L - e, L - s
    # extend downstream to stop
    j = e
    while j + 3 <= L and seq[j-3:j] not in STOPS:
        j += 3
        if j - e > 3000:
            return None
    # extend upstream: predicted start position from query offset
    target = s - 3 * (qstart - 1)
    best, i = None, s
    while i - 3 >= 0 and seq[i-3:i] not in STOPS and s - i < 3000:
        i -= 3
        if seq[i:i+3] in STARTS:
            score = abs(i - target) + (0 if seq[i:i+3] == 'ATG' else 30)
            if best is None or score < best[0]:
                best = (score, i)
    if seq[s:s+3] in STARTS:
        score = abs(s - target) + (0 if seq[s:s+3] == 'ATG' else 30)
        if best is None or score < best[0]:
            best = (score, s)
    if best is None:
        return None
    a = best[1]
    if strand == 1:
        return a, j
    return L - j, L - a


def homology_annotate(record, ref_records, label):
    genome = str(record.seq).upper()
    n = len(genome)
    doubled = genome + genome[:20000]
    prot, rrna = reference_sets(ref_records)
    feats = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / 'g.fa').write_text(f'>g\n{doubled}\n')
        (td / 'p.fa').write_text(''.join(f'>{g}\n{s}\n' for g, s in prot.items()))
        (td / 'r.fa').write_text(''.join(f'>{g}\n{s}\n' for g, s in rrna.items()))
        _run([str(BIN/'makeblastdb'), '-in', str(td/'g.fa'), '-dbtype', 'nucl', '-out', str(td/'g')])
        fmt = '6 qseqid sstart send evalue bitscore qstart qend qlen pident'
        _run([str(BIN/'tblastn'), '-query', str(td/'p.fa'), '-db', str(td/'g'), '-db_gencode', '11',
              '-evalue', '1e-8', '-outfmt', fmt, '-out', str(td/'p.tsv'), '-seg', 'no'])
        _run([str(BIN/'blastn'), '-query', str(td/'r.fa'), '-db', str(td/'g'), '-evalue', '1e-20',
              '-outfmt', fmt, '-out', str(td/'r.tsv')])
        hits = defaultdict(list)
        for line in open(td/'p.tsv'):
            q, ss, se, ev, bs, qs, qe, ql, pid = line.split('\t')
            hits[q].append((float(bs), int(ss), int(se), int(qs), int(qe), int(ql), float(pid)))
        used = []
        for g, hs in hits.items():
            hs.sort(reverse=True)
            # cluster HSPs belonging to best locus (same strand, within 6 kb)
            bs, ss, se, qs, qe, ql, pid = hs[0]
            strand = 1 if se >= ss else -1
            lo, hi = min(ss, se) - 1, max(ss, se)
            locus = [h for h in hs if (1 if h[2] >= h[1] else -1) == strand and
                     min(h[1], h[2]) - 1 < hi + 6000 and max(h[1], h[2]) > lo - 6000]
            cover = sum(h[4] - h[3] + 1 for h in locus) / ql
            if cover < 0.4:
                continue
            # choose HSP containing query N-terminus (smallest qstart) and C-terminus
            nterm = min(locus, key=lambda h: h[3])
            cterm = max(locus, key=lambda h: h[4])
            if strand == 1:
                s, e = min(nterm[1], nterm[2]) - 1, max(cterm[1], cterm[2])
            else:
                s, e = min(cterm[1], cterm[2]) - 1, max(nterm[1], nterm[2])
            # single-exon ORF refinement when one HSP covers most of the protein
            single = len(locus) == 1 or (cterm is nterm)
            orf = _orf_from_hit(doubled, strand, s, e, nterm[3], ql, n) if single else None
            if orf:
                a, b = orf
                note = 'single-HSP ORF extended to start/stop codons'
            else:
                if e <= s:  # HSP order inconsistent (e.g. trans-split gene): use locus span
                    s = min(min(h[1], h[2]) for h in locus) - 1
                    e = max(max(h[1], h[2]) for h in locus)
                a, b = s, e
                note = 'multi-HSP locus (possible intron/split gene): boundaries approximate'
            a %= n
            b = a + (b - orf[0] if orf else e - s)
            if any(abs(a - u) < 30 and g2 == g for g2, u in used):
                continue
            used.append((g, a))
            loc = FeatureLocation(a, min(b, n), strand) if b <= n else CompoundLocation(
                [FeatureLocation(a, n, strand), FeatureLocation(0, b - n, strand)] if strand == 1 else
                [FeatureLocation(0, b - n, strand), FeatureLocation(a, n, strand)])
            fq = {'gene': [g], 'gene_raw': [g], 'inference': ['homology_transfer_v3'],
                  'note': [f'{label}: tblastn vs reference protein; {note}; cover={cover:.2f}; pid={pid:.1f}'],
                  'transl_table': ['11'], 'codon_start': ['1']}
            f = SeqFeature(loc, type='CDS', qualifiers=fq)
            if orf:
                cds = str(f.extract(record.seq)).upper()
                aa = str(Seq(cds[:len(cds)//3*3]).translate(table=11))
                fq['translation'] = [('M' + aa[1:]).rstrip('*')]
                fq['note'][0] += f'; internal_stops={aa[:-1].count("*")}'
            else:
                fq['pseudo'] = ['']  # do not treat as SD-eligible
            feats.append(f)
        for line in open(td/'r.tsv'):
            q, ss, se, ev, bs, qs, qe, ql, pid = line.split('\t')
            ss, se, ql = int(ss), int(se), int(ql)
            if abs(se - ss) + 1 < 0.6 * ql:
                continue
            strand = 1 if se >= ss else -1
            a = (min(ss, se) - 1) % n
            b = a + abs(se - ss) + 1
            if b > n or any(abs(a - u) < 50 and g2 == q for g2, u in used):
                continue
            used.append((q, a))
            feats.append(SeqFeature(FeatureLocation(a, b, strand), type='rRNA', qualifiers={
                'gene': [q], 'gene_raw': [q], 'inference': ['homology_transfer_v3'],
                'note': [f'{label}: blastn vs reference rRNA; pid={float(pid):.1f}; approximate ends']}))
        # tRNAs with ARAGORN (circular, genetic code 11)
        (td/'c.fa').write_text(f'>g\n{genome}\n')
        out = subprocess.run([str(BIN/'aragorn'), '-t', '-gc11', '-c', '-w', str(td/'c.fa')],
                             stdout=subprocess.PIPE, text=True).stdout
        for line in out.splitlines():
            m = re.match(r'\s*\d+\s+tRNA-(\w+)\s+(c?)\[(\d+),(\d+)\]\s+\S+\s+\((\w+)\)', line)
            if m:
                aa, comp, s1, e1, ac = m.groups()
                s1, e1 = int(s1), int(e1)
                if e1 < s1 or e1 > n:
                    continue
                strand = -1 if comp else 1
                feats.append(SeqFeature(FeatureLocation(s1 - 1, e1, strand), type='tRNA', qualifiers={
                    'gene': [f'trn{aa[0] if len(aa)==1 else aa}-{ac.upper().replace("T","U")}'],
                    'inference': ['ARAGORN -t -gc11 -c'], 'note': [f'{label}: ARAGORN tRNA prediction']}))
    feats.sort(key=lambda f: int(f.location.start))
    new = SeqRecord(record.seq, id=record.id, name=record.name, description=record.description,
                    annotations=dict(record.annotations), features=[x for x in record.features if x.type == 'source'] + feats)
    new.annotations['molecule_type'] = 'DNA'
    return new
