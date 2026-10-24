#!/usr/bin/env python3
"""Uniform per-genome regulatory-signal scan (comparative workflow v3).

For every CDS/rRNA boundary context of a standardised plastome record this
module measures, with identical parameters in every genome:

  upstream (transcript-oriented, relative to the first transcribed base):
    * contiguous non-coding stretch length (no CDS/rRNA/tRNA exon on either strand)
    * Chlamydomonas-type -10 octamer TATAATAT (exact), TATAAT hexamer (exact),
      -35 TTGACA (exact), and paired PEP-like -35/-10 boxes (<=1 mismatch each,
      15-19 nt spacer) -- scanned ONLY inside the non-coding stretch (<=300 nt)
    * SD-like GGAG-family words, 3-15 nt spacer before the start codon (CDS only)
    * start-codon-region folding (-30..+30; ViennaRNA MFE) and dinucleotide-
      shuffle z-score (translation-initiation accessibility)
  downstream:
    * whole-window MFE of +1..+120 and its dinucleotide-shuffle z-score
    * strongest perfect antiparallel stem-loop (v2 enumerator, Turner 2004)

Null model: randomised Euler-path dinucleotide-preserving shuffles of each
window (approximate, non-uniform; see v2 METHODS). The SAME shuffles of every
anchor are summed into genome-level replicate totals, giving genome-level
observed/expected ratios and empirical P-values for motif enrichment and for
positional profiles. Background control windows are drawn from distal
non-coding DNA (>300 nt from any annotated start) to estimate positional
specificity.
"""
from __future__ import annotations

import random
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import RNA

import chloroplast_regulatory_inventory_improved as v2

UP, DOWN, SDW, STARTW = 300, 120, 35, 30
OCT = re.compile('(?=TATAATAT)')
HEX = re.compile('(?=TATAAT)')
M35 = re.compile('(?=TTGACA)')
SD_WORDS = ('GGAGG', 'AGGAG', 'GGAG', 'GAGG')
BINS = np.arange(-300, 1, 25)
KMER_K = 8


def count(pat, s):
    return sum(1 for _ in pat.finditer(s))


def positions(pat, s, offset):
    return [m.start() + offset for m in pat.finditer(s)]


def strong_pairs(s):
    """PEP-like pairs with <=1 total mismatch and 16-18 nt spacer."""
    return [h for h in v2.promoter_pairs(s) if h[4] + h[5] <= 1 and 16 <= h[3] <= 18]


def sd_present(short):
    return any(len(h[2]) >= 4 and h[2] in SD_WORDS + ('AGGAGG', 'AGGA', 'GGAG') for h in v2.sd_hits(short))


def mfe(seq, md):
    if not seq:
        return np.nan, ''
    fc = RNA.fold_compound(seq.replace('T', 'U'), md)
    st, e = fc.mfe()
    return float(e), st


def kmers(s, k=KMER_K):
    return Counter(s[i:i+k] for i in range(len(s)-k+1) if set(s[i:i+k]) <= set('ACGT'))


def prepare(record):
    features, products, anchors, qc, translations, occ, cds = v2.make_annotations(record)
    excluded = {f'F{i:04d}' for i, f in enumerate(record.features, 1)
                if 'exclude_from_boundary_analysis' in f.qualifiers}
    anchors = anchors[anchors.feature_type.isin(['CDS', 'rRNA']) & ~anchors.feature_id.isin(excluded)].copy()
    return features, products, anchors, qc, translations, occ, cds


def scan_genome(record, acc, nrep=200, nrep_fold=100, nrep_start=50, temperature=25.0, seed=20261023,
                keep_all_hairpins=False):
    rng = random.Random(f'{seed}-{acc}')
    md = RNA.md(); md.temperature = temperature
    genome = str(record.seq).upper(); n = len(genome)
    circular = record.annotations.get('topology', 'linear').lower() == 'circular'
    features, products, anchors, qc, translations, occ, cds = prepare(record)
    noncoding = np.array([len(x) == 0 for x in occ])
    starts0 = set()
    rows, hits, hps = [], [], []
    rep = {k: np.zeros(nrep) for k in ('oct', 'hex', 'm35', 'pep', 'sd', 'oct_bg', 'hex_bg')}
    obs = Counter()
    prof_obs = np.zeros(len(BINS)-1); prof_null = np.zeros((nrep, len(BINS)-1))
    km_obs, km_null = Counter(), Counter()
    total_nc = 0
    for a in anchors.to_dict('records'):
        first0, last0, strand = a['first_transcribed_base_1b']-1, a['last_transcribed_base_1b']-1, a['strand']
        starts0.add(first0)
        up, uppos = v2.windows(genome, first0, strand, range(-UP, 0), circular)
        # contiguous non-coding stretch immediately upstream
        L = 0
        for x in reversed(uppos):
            if not noncoding[x]:
                break
            L += 1
        nc = up[len(up)-L:]; off = -L
        total_nc += L
        down, downpos = v2.windows(genome, last0, strand, range(1, DOWN+1), circular)
        Ld = 0
        for x in downpos:
            if not noncoding[x]:
                break
            Ld += 1
        o_oct, o_hex, o_m35 = positions(OCT, nc, off), positions(HEX, nc, off), positions(M35, nc, off)
        pairs = strong_pairs(nc)
        is_cds = a['feature_type'] == 'CDS'
        short, _ = v2.windows(genome, first0, strand, range(-SDW, 0), circular)
        sdh = v2.sd_hits(short) if (is_cds and a['sd_eligible']) else []
        sd_ok = bool(sdh)
        obs['oct'] += len(o_oct); obs['hex'] += len(o_hex); obs['m35'] += len(o_m35)
        obs['pep'] += len(pairs); obs['sd'] += int(sd_ok)
        for p in o_oct:
            b = np.searchsorted(BINS, p, side='right') - 1
            prof_obs[b] += 1
        km_obs.update(kmers(nc[-100:]))
        # shuffles of the non-coding stretch and of the SD window
        n_oct_null = np.zeros(nrep)
        for r in range(nrep):
            s = v2.dinuc_shuffle(nc, rng) if L > 2 else nc
            po = positions(OCT, s, off)
            n_oct_null[r] = len(po)
            rep['oct'][r] += len(po); rep['hex'][r] += count(HEX, s); rep['m35'][r] += count(M35, s)
            rep['pep'][r] += len(strong_pairs(s))
            for p in po:
                prof_null[r, np.searchsorted(BINS, p, side='right') - 1] += 1
            if r < 20:
                km_null.update(kmers(s[-100:]))
            if is_cds and a['sd_eligible']:
                rep['sd'][r] += int(bool(v2.sd_hits(v2.dinuc_shuffle(short, rng))))
        # folding: start region and downstream window
        st_e = st_z = st_p = np.nan
        if is_cds:
            sreg, _ = v2.windows(genome, first0, strand, range(-STARTW, STARTW+1), circular)
            st_e, _ = mfe(sreg, md)
            nulls = [mfe(v2.dinuc_shuffle(sreg, rng), md)[0] for _ in range(nrep_start)]
            st_z = (st_e - np.mean(nulls)) / (np.std(nulls, ddof=1) or np.nan)
            st_p = (1 + sum(x <= st_e for x in nulls)) / (nrep_start + 1)  # lower = more structured
            st_p_open = (1 + sum(x >= st_e for x in nulls)) / (nrep_start + 1)
        d_e, d_st = mfe(down, md)
        dn = [mfe(v2.dinuc_shuffle(down, rng), md)[0] for _ in range(nrep_fold)]
        d_z = (d_e - np.mean(dn)) / (np.std(dn, ddof=1) or np.nan)
        d_p = (1 + sum(x <= d_e for x in dn)) / (nrep_fold + 1)
        best = None; allh = []
        for h in v2.hairpins(down):
            core, intended, structure, mf, energy, retained = v2.fold_hairpin(down, h, md)
            if energy > -2:
                continue
            x0, y0, stem, loop, l, rr = h
            tail = down[y0:y0+8]
            rec = dict(anchor_id=a['anchor_id'], gene=a['gene'], start_after_stop=x0+1, end_after_stop=y0,
                       stem_pairs=stem, loop_nt=loop, left_arm=down[x0:l], loop=down[l:rr], right_arm=down[rr:y0],
                       intended_energy=energy, core_MFE=mf, MFE_pair_retention=retained,
                       tail8=tail, tail_T=tail.count('T'),
                       overlaps_coding=bool(any(not noncoding[p] for p in downpos[x0:y0])),
                       **v2.path_fields(downpos[x0:y0], n))
            allh.append(rec)
            if best is None or energy < best['intended_energy']:
                best = rec
        if keep_all_hairpins:
            hps.extend(allh)
        elif best:
            hps.append(best)
        for kind, plist, w in (('TATAATAT', o_oct, 8), ('TATAAT', o_hex, 6), ('TTGACA', o_m35, 6)):
            for p in plist:
                i = len(up) + p
                hits.append(dict(anchor_id=a['anchor_id'], gene=a['gene'], feature_type=a['feature_type'], motif=kind,
                                 rel_to_first_base=p, sequence=up[i:i+w], **v2.path_fields(uppos[i:i+w], n)))
        for (i, j_end, j, spacer, mm35, mm10, score) in pairs:
            ii = len(up) - L + i
            hits.append(dict(anchor_id=a['anchor_id'], gene=a['gene'], feature_type=a['feature_type'], motif='PEP_pair',
                             rel_to_first_base=off + i, sequence=up[ii:ii + (j_end - i)], spacer=spacer,
                             mismatches=mm35 + mm10, minus35=nc[i:i+6], minus10=nc[j:j+6],
                             **v2.path_fields(uppos[ii:ii + (j_end - i)], n)))
        for (s0, e0, motif, spacer, score) in sdh:
            hits.append(dict(anchor_id=a['anchor_id'], gene=a['gene'], feature_type='CDS', motif='SD_like',
                             rel_to_first_base=s0 - len(short), sequence=motif, spacer=spacer))
        rows.append(dict(genome=acc, anchor_id=a['anchor_id'], product_id=a['product_id'], gene=a['gene'],
                         feature_type=a['feature_type'], strand=strand, context_kind=a['context_kind'],
                         first_base_1b=first0+1, last_base_1b=last0+1, start_codon=genome[first0:first0+3] if strand == 1 else v2.rc(genome[first0-2:first0+1]) if is_cds else '',
                         upstream_noncoding_len=L, downstream_noncoding_len=Ld,
                         upstream300=up, nc_upstream=nc,
                         n_TATAATAT=len(o_oct), TATAATAT_pos=';'.join(map(str, o_oct)),
                         n_TATAAT=len(o_hex), TATAAT_pos=';'.join(map(str, o_hex)),
                         n_TTGACA=len(o_m35), n_PEP_pairs=len(pairs),
                         best_PEP=(f'{nc[pairs[0][0]:pairs[0][0]+6]}-{pairs[0][3]}-{nc[pairs[0][2]:pairs[0][2]+6]}@{off+pairs[0][0]}' if pairs else ''),
                         octamer_null_mean=float(n_oct_null.mean()),
                         octamer_P=(1 + int((n_oct_null >= len(o_oct)).sum())) / (nrep + 1) if o_oct else 1.0,
                         SD_like=sd_ok, SD_words=';'.join(f'{h[2]}@{h[0]-len(short)}' for h in sdh),
                         start_region_MFE=st_e, start_region_z=st_z, start_P_structured=st_p,
                         start_P_open=(st_p_open if is_cds else np.nan),
                         down120=down, down_MFE=d_e, down_MFE_structure=d_st, down_z=d_z, down_P=d_p,
                         best_hairpin_energy=best['intended_energy'] if best else np.nan,
                         best_hairpin_stem=best['stem_pairs'] if best else np.nan,
                         best_hairpin_start=best['start_after_stop'] if best else np.nan,
                         best_hairpin_tailT=best['tail_T'] if best else np.nan))
    df = pd.DataFrame(rows)
    if len(df):
        df['down_q_BH'] = v2.bh_adjust(df.down_P.values)
        cds_mask = df.feature_type == 'CDS'
        df.loc[cds_mask, 'start_open_q_BH'] = v2.bh_adjust(df.loc[cds_mask, 'start_P_open'].values)
        df['octamer_q_BH'] = v2.bh_adjust(df.octamer_P.values)
    # distal non-coding background (>300 nt from any boundary start) for specificity
    near = np.zeros(n, dtype=bool)
    for s0 in starts0:
        for d in range(-UP, UP+1):
            near[(s0 + d) % n] = True
    bgmask = noncoding & ~near
    runs, i = [], 0
    while i < n:
        if bgmask[i]:
            j = i
            while j < n and bgmask[j]:
                j += 1
            runs.append((i, j)); i = j
        else:
            i += 1
    bg_len = sum(j - i for i, j in runs)
    bg_oct = sum(count(OCT, genome[i:j]) + count(OCT, v2.rc(genome[i:j])) for i, j in runs)
    bg_hex = sum(count(HEX, genome[i:j]) + count(HEX, v2.rc(genome[i:j])) for i, j in runs)
    bg_rows = []
    long_runs = [(i, j) for i, j in runs if j - i >= DOWN]
    for k in range(min(60, len(long_runs) * 3)):
        i, j = long_runs[rng.randrange(len(long_runs))]
        s0 = rng.randrange(i, j - DOWN + 1)
        w = genome[s0:s0+DOWN]
        if rng.random() < .5:
            w = v2.rc(w)
        e, _ = mfe(w, md)
        nn = [mfe(v2.dinuc_shuffle(w, rng), md)[0] for _ in range(nrep_fold // 2)]
        bg_rows.append(dict(genome=acc, bg_window=k+1, start_1b=s0+1, MFE=e,
                            z=(e - np.mean(nn)) / (np.std(nn, ddof=1) or np.nan)))
    summ = dict(genome=acc, n_anchors=len(df), n_CDS_anchors=int((df.feature_type == 'CDS').sum()) if len(df) else 0,
                upstream_nc_bp=total_nc, distal_nc_bp=bg_len, distal_TATAATAT=bg_oct, distal_TATAAT=bg_hex,
                distal_TATAATAT_per_kb=1000*bg_oct/(2*bg_len) if bg_len else np.nan,
                upstream_TATAATAT_per_kb=1000*obs['oct']/total_nc if total_nc else np.nan,
                genome_GC=100*(genome.count('G')+genome.count('C'))/n,
                upstream_nc_GC=100*sum(r['nc_upstream'].count('G')+r['nc_upstream'].count('C') for r in rows)/max(total_nc, 1))
    for k in ('oct', 'hex', 'm35', 'pep', 'sd'):
        mu = rep[k].mean()
        summ[f'{k}_observed'] = obs[k]; summ[f'{k}_expected'] = mu
        summ[f'{k}_OE'] = obs[k] / mu if mu else np.nan
        summ[f'{k}_P_enriched'] = (1 + int((rep[k] >= obs[k]).sum())) / (nrep + 1)
        summ[f'{k}_P_depleted'] = (1 + int((rep[k] <= obs[k]).sum())) / (nrep + 1)
    if len(df):
        summ['down_z_median'] = float(df.down_z.median()); summ['down_frac_z_le_m2'] = float((df.down_z <= -2).mean())
        summ['down_frac_q05'] = float((df.down_q_BH < .05).mean())
        c = df[df.feature_type == 'CDS']
        summ['start_z_median'] = float(c.start_region_z.median())
        summ['start_frac_z_ge_1'] = float((c.start_region_z >= 1).mean())
        summ['SD_fraction_CDS'] = float(c.SD_like.mean())
    bgdf = pd.DataFrame(bg_rows)
    if len(bgdf):
        summ['bg_down_z_median'] = float(bgdf.z.median()); summ['bg_frac_z_le_m2'] = float((bgdf.z <= -2).mean())
    prof = pd.DataFrame(dict(genome=acc, bin_start=BINS[:-1], bin_end=BINS[1:]-1, observed=prof_obs,
                             expected=prof_null.mean(0), P_enriched=[(1+(prof_null[:, b] >= prof_obs[b]).sum())/(nrep+1)
                                                                       for b in range(len(prof_obs))]))
    kmer_rows = []
    scale = nrep / min(nrep, 20)
    for kmer, o in km_obs.items():
        e = km_null.get(kmer, 0) / min(nrep, 20)
        kmer_rows.append(dict(genome=acc, kmer=kmer, observed=o, expected=e, OE=(o + .5) / (e + .5)))
    km = pd.DataFrame(kmer_rows).sort_values('OE', ascending=False)
    return dict(anchors=df, hits=pd.DataFrame(hits), hairpins=pd.DataFrame(hps), summary=pd.DataFrame([summ]),
                profile=prof, kmers=km, background=bgdf, qc=qc.assign(genome=acc),
                translations=translations.assign(genome=acc), products=products.assign(genome=acc))
