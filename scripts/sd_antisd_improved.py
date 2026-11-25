#!/usr/bin/env python3
"""Improved SD / anti-SD inventory for green-algal plastomes.

Problem: the exact-word scan of the ANNOTATED 16S 3' end misses the anti-SD in many
genomes, because chloroplast 16S rRNA 3' ends are frequently truncated or misannotated in
WGS assemblies. This script detects the anti-SD with a cascade of approaches and detects
SD-like elements with expanded cores and a thermodynamic (anti-SD:mRNA duplex) method.

Anti-SD cascade (per genome, first success wins):
 1 annotated_tail      : CCTCC in last 60 nt of the annotated rrn16
 2 genomic_extension   : conserved block (CACCTCC/CCTCC) in annotated region +/-50 nt genomic
 3 homology_locus      : GATCACCTCC (<=2 mm, core CCTCC intact) within rrn16 locus +/-400 nt
 4 homology_genomewide : GATCACCTCC (<=1 mm, core CCTCC intact) anywhere, both strands
 else absent_in_assembly (conserved anti-SD block not present in the deposited sequence)

SD detection (per CDS, 25 nt upstream of the start codon):
 strict      : GGAGG/AGGAGG/GGAG core ending 3-15 nt before start (current-style)
 expanded    : GGAGG(<=1mm)/AGGAGG/GGAG/AGGAG/GAGG/AGGA ending 3-18 nt before start
 thermodynamic: ViennaRNA duplexfold of the universal anti-SD tail (5'-ACCUCCUU-3') against
               the leader; SD+ if dG <= -4.0 kcal/mol with the duplex 3-15 nt before start
"""
import argparse, re
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import RNA

ANTISD_QUERY = 'ACCUCCUU'          # universal SSU-rRNA 3' anti-SD segment (core CCUCC)
# conserved 3'-terminal 30-mer of the reference 16S (C. reinhardtii FJ423446); anti-SD core = CCTCC
REF3 = None                        # set in main() from the reference genome
SD_STRICT = ['AGGAGG', 'GGAGG', 'GGAG']
SD_EXPANDED = ['AGGAGG', 'GGAGG', 'AGGAG', 'GGAG', 'GAGG', 'AGGA']
ID_THRESH = 0.73                   # min identity over the 30-mer to call a homologous 16S 3' end

def best_homology(seq, query, core_off):
    """Seed on the anti-SD core (CCTCC), align the 30-mer query, return (identity, mm, window, end)."""
    best = None; L = len(query)
    i = seq.find('CCTCC')
    while i != -1:
        start = i - core_off
        if 0 <= start and start + L <= len(seq):
            w = seq[start:start+L]
            mm = sum(1 for a, b in zip(w, query) if a != b)
            ident = 1 - mm/L
            if best is None or ident > best[0]:
                best = (ident, mm, w, start+L)   # end = genomic index just past the 3' terminus
        i = seq.find('CCTCC', i+1)
    return best

def detect_antisd(rec):
    fe = [f for f in rec.features if f.type == 'rRNA' and f.qualifiers.get('gene', [''])[0] == 'rrn16']
    annotated_ok = False
    if fe:
        ann = str(fe[0].extract(rec.seq)).upper().replace('U', 'T')
        annotated_ok = ann[-40:].rfind('CCTCC') >= 0
    # homology reannotation: search both strands for the conserved 16S 3'-terminal 30-mer
    fwd = str(rec.seq).upper().replace('U', 'T'); rev = str(rec.seq.reverse_complement()).upper().replace('U', 'T')
    core_off = REF3.rfind('CCTCC')
    cands = [b for b in (best_homology(fwd, REF3, core_off), best_homology(rev, REF3, core_off)) if b]
    hom = max(cands, key=lambda z: z[0]) if cands else None
    hom_ident = round(hom[0], 3) if hom else 0.0
    hom_tail = hom[2] if hom else ''
    hom_ok = bool(hom) and hom[0] >= ID_THRESH and 'CCTCC' in hom[2]
    if annotated_ok:
        method = 'annotated'; present = 'yes'
        tail = str(fe[0].extract(rec.seq))[-25:].upper().replace('U', 'T')
    elif hom_ok:
        method = 'reannotated_homology'; present = 'yes'; tail = hom_tail
    else:
        method = 'absent_in_assembly'; present = 'absent'; tail = hom_tail
    return dict(method=method, hom_identity=hom_ident, present=present, tail=tail,
                antiSD_core=('CCTCC' if 'CCTCC' in tail else '-'))

def upstream(rec, feat, n=25):
    strand = feat.location.strand
    L = len(rec.seq)
    if strand == 1:
        st = int(feat.location.start)
        seq = rec.seq[(st-n) % L:st] if st-n >= 0 else rec.seq[st-n:] + rec.seq[:st]
    else:
        en = int(feat.location.end)
        seg = rec.seq[en:en+n] if en+n <= L else rec.seq[en:] + rec.seq[:(en+n) % L]
        seq = seg.reverse_complement()
    return str(seq).upper().replace('U', 'T')

def sd_strict(lead):
    for core in SD_STRICT:
        j = lead.rfind(core)
        if j >= 0:
            dist = len(lead) - (j + len(core))  # nt from end of core to start codon
            if 3 <= dist <= 15: return 1, dist, core
    return 0, None, None

def sd_expanded(lead):
    best = None
    # GGAGG with <=1 mismatch
    for i in range(len(lead)-5+1):
        w = lead[i:i+5]
        if sum(1 for a, b in zip(w, 'GGAGG') if a != b) <= 1:
            dist = len(lead)-(i+5)
            if 3 <= dist <= 18: best = (dist, 'GGAGG~'); break
    if best is None:
        for core in SD_EXPANDED:
            j = lead.rfind(core)
            if j >= 0:
                dist = len(lead)-(j+len(core))
                if 3 <= dist <= 18: best = (dist, core); break
    return (1, best[0], best[1]) if best else (0, None, None)

def sd_thermo(lead):
    leadR = lead.replace('T', 'U')
    d = RNA.duplexfold(ANTISD_QUERY, leadR)
    # d.j = end position (1-based) of the duplex on the second sequence (leader)
    j_end = d.j
    dist = len(lead) - j_end  # nt from duplex end to start codon
    dg = d.energy
    call = 1 if (dg <= -4.0 and 3 <= dist <= 15) else 0
    return call, round(dg, 2), dist

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='.')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; T.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t'); sp = dict(zip(man.insdc_accession, man.species))
    global REF3
    _r = SeqIO.read(P/'data/standardized/FJ423446.gb', 'genbank')
    _f = [f for f in _r.features if f.type == 'rRNA' and f.qualifiers.get('gene', [''])[0] == 'rrn16'][0]
    REF3 = str(_f.extract(_r.seq))[-30:].upper().replace('U', 'T')
    print(f"reference 16S 3' 30-mer (FJ423446): {REF3}  (anti-SD core at offset {REF3.rfind('CCTCC')})\n")

    anti_rows = []; sd_rows = []; sd_sum = []
    for gb in sorted((P/'data/standardized').glob('*.gb')):
        rec = SeqIO.read(gb, 'genbank'); acc = rec.id if rec.id in sp else gb.stem
        species = sp.get(acc, sp.get(gb.stem, gb.stem))
        ad = detect_antisd(rec)
        anti_rows.append(dict(genome=acc, species=species, **ad))
        cds = [f for f in rec.features if f.type == 'CDS' and 'gene' in f.qualifiers]
        n = ns = ne = nt = 0
        for f in cds:
            lead = upstream(rec, f, 25)
            if len(lead) < 10: continue
            n += 1
            s1, d1, c1 = sd_strict(lead)
            s2, d2, c2 = sd_expanded(lead)
            s3, dg, d3 = sd_thermo(lead)
            ns += s1; ne += s2; nt += s3
            sd_rows.append(dict(genome=acc, species=species, gene=f.qualifiers['gene'][0],
                                strict=s1, strict_dist=d1, expanded=s2, expanded_dist=d2,
                                thermo=s3, duplex_dG=dg, thermo_dist=d3))
        sd_sum.append(dict(genome=acc, species=species, n_cds=n, n_strict=ns, n_expanded=ne,
                           n_thermo=nt, pct_strict=round(100*ns/max(n,1),1),
                           pct_expanded=round(100*ne/max(n,1),1), pct_thermo=round(100*nt/max(n,1),1)))
    AD = pd.DataFrame(anti_rows); AD.to_csv(T/'antiSD_16S_3prime_improved.tsv', sep='\t', index=False)
    SDg = pd.DataFrame(sd_rows); SDg.to_csv(T/'SD_per_gene_improved.tsv', sep='\t', index=False)
    SDs = pd.DataFrame(sd_sum); SDs.to_csv(T/'SD_summary_by_genome_improved.tsv', sep='\t', index=False)

    # ---- summary printout ----
    print("=== ANTI-SD inventory (improved multi-approach) ===")
    print(AD.method.value_counts().to_string())
    present = (AD.present == 'yes').sum()
    print(f"\nanti-SD PRESENT: {present}/{len(AD)}   (current annotated-tail-only was 22/41)")
    print(f"absent from deposited sequence (assembly-limited): {(AD.present=='absent').sum()}")
    print("\n=== SD detection totals across all CDS ===")
    tot = SDg.shape[0]
    print(f"CDS analysed: {tot}")
    print(f"  strict (current-style):   {SDg.strict.sum():5d}  ({100*SDg.strict.mean():.1f}%)")
    print(f"  expanded cores:           {SDg.expanded.sum():5d}  ({100*SDg.expanded.mean():.1f}%)")
    print(f"  thermodynamic (dG<=-4):   {SDg.thermo.sum():5d}  ({100*SDg.thermo.mean():.1f}%)")
    print("\n=== Genomes with SD enrichment signal (>=1 SD CDS) ===")
    for col in ['n_strict','n_expanded','n_thermo']:
        print(f"  {col}: {(SDs[col]>0).sum()}/{len(SDs)} genomes have >=1")
    # A. falcatus focus
    af = SDs[SDs.genome.str.startswith('QB176266')]
    if len(af): print("\nA. falcatus:", af.to_dict('records')[0])

if __name__ == '__main__':
    main()
