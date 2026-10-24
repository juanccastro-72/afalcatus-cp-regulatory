#!/usr/bin/env python3
"""Aggregate per-genome scans into comparative tables and statistics."""
import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import chloroplast_regulatory_inventory_improved as v2

FOCAL = 'QB176266.1'
CORE = ['psbA', 'psbB', 'psbC', 'psbD', 'psbE', 'psbH', 'psbK', 'psbT', 'psaA', 'psaB', 'psaC', 'psaJ', 'petA', 'petB', 'petD',
        'petG', 'petL', 'atpA', 'atpB', 'atpE', 'atpF', 'atpH', 'atpI', 'rbcL', 'ccsA', 'cemA', 'chlB', 'chlL', 'chlN', 'clpP',
        'ftsH', 'infA', 'rpl2', 'rpl5', 'rpl14', 'rpl16', 'rpl20', 'rpl23', 'rpl36', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8',
        'rps9', 'rps11', 'rps12', 'rps14', 'rps18', 'rps19', 'rpoA', 'rpoBa', 'rpoC1', 'rpoC2', 'tufA', 'ycf3', 'ycf4',
        'ycf12', 'rrn16', 'rrn23']
VALIDATED_PART_GENES = {  # C. reinhardtii orthologs used as validated expression parts (curated rows)
    'psbA': 'CE057/CE058/CE059/CE062', 'psbD': 'CE057/CE059/CE060', 'atpA': 'CE057/CE058/CE059/CE062',
    'rbcL': 'CE057/CE064', 'rrn16': 'CE058', 'psaA': 'CE061', 'atpB': 'CE004/CE050', 'petD': 'CE022-CE024'}


def load(P, tag, accs):
    out = {}
    for acc in accs:
        out[acc] = pickle.load(open(P/f'results/scan_cache/{tag}/{acc}.pkl', 'rb'))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..'); ap.add_argument('--tag', default='main')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t'); man = man[man.included_in_comparative_set]
    arch = pd.read_csv(T/'plastome_architecture.tsv', sep='\t')
    res = load(P, a.tag, man.insdc_accession)
    S = pd.concat([r['summary'] for r in res.values()]).merge(arch, on='genome')
    S = S.merge(man[['insdc_accession', 'role']], left_on='genome', right_on='insdc_accession').drop(columns='insdc_accession')
    S['group'] = np.where(S.clade == 'Sphaeropleales', 'Sphaeropleales',
                          np.where(S.clade == 'Chlamydomonadales', 'Chlamydomonadales', 'Outgroups (OCC/Trebouxiophyceae)'))
    for k in ('oct', 'hex', 'm35', 'pep', 'sd'):
        S[f'{k}_q_enriched_BH_across_genomes'] = v2.bh_adjust(S[f'{k}_P_enriched'].values)
    S.to_csv(T/'genome_level_signal_summary.tsv', sep='\t', index=False)
    A = pd.concat([r['anchors'] for r in res.values()])
    A = A.merge(S[['genome', 'species', 'family', 'clade', 'group']], on='genome')
    # start-region accessibility only for true translation starts (excludes internal trans-spliced exons)
    valid = A[(A.feature_type == 'CDS') & A.start_codon.isin(['ATG', 'GTG', 'TTG'])]
    A.loc[(A.feature_type == 'CDS') & ~A.start_codon.isin(['ATG', 'GTG', 'TTG']), ['start_region_z', 'start_region_MFE']] = np.nan
    sz = valid.groupby('genome').start_region_z.agg(['median', lambda x: (x >= 1).mean(), 'count'])
    sz.columns = ['start_z_median', 'start_frac_z_ge_1', 'n_start_contexts']
    S = S.drop(columns=['start_z_median', 'start_frac_z_ge_1']).merge(sz, left_on='genome', right_index=True)
    S.to_csv(T/'genome_level_signal_summary.tsv', sep='\t', index=False)
    A.to_csv(T/'all_boundary_contexts.tsv.gz', sep='\t', index=False)
    H = pd.concat([r['hits'].assign(genome=g) for g, r in res.items() if len(r['hits'])])
    H.to_csv(T/'all_motif_hits.tsv.gz', sep='\t', index=False)
    HP = pd.concat([r['hairpins'].assign(genome=g) for g, r in res.items() if len(r['hairpins'])])
    HP.to_csv(T/'all_best_hairpins.tsv.gz', sep='\t', index=False)
    PR = pd.concat([r['profile'] for r in res.values()]); PR.to_csv(T/'octamer_positional_profiles.tsv', sep='\t', index=False)
    BG = pd.concat([r['background'] for r in res.values()]); BG.to_csv(T/'distal_background_windows.tsv', sep='\t', index=False)
    KM = pd.concat([r['kmers'] for r in res.values()])
    # cross-genome k-mer meta-analysis (8-mers, -100..-1 non-coding)
    KM['enriched'] = (KM.OE >= 2) & (KM.observed >= 3)
    meta = KM.groupby('kmer').agg(n_genomes_enriched=('enriched', 'sum'), total_observed=('observed', 'sum'),
                                  total_expected=('expected', 'sum')).reset_index()
    meta['pooled_OE'] = (meta.total_observed + .5) / (meta.total_expected + .5)
    meta = meta.sort_values(['n_genomes_enriched', 'pooled_OE'], ascending=False)
    meta.head(500).to_csv(T/'kmer8_cross_genome_enrichment.tsv', sep='\t', index=False)
    KM[KM.genome == FOCAL].head(200).to_csv(T/'kmer8_Afalcatus_top.tsv', sep='\t', index=False)
    # Statistical comparisons
    st = []
    def add(test, x, y=None, **kw):
        st.append(dict(test=test, **kw))
    sph = S[S.group == 'Sphaeropleales']; chl = S[S.group == 'Chlamydomonadales']
    for k, lab in (('oct', 'TATAATAT O/E'), ('hex', 'TATAAT O/E'), ('pep', 'PEP pair O/E'), ('sd', 'SD-like O/E')):
        u = stats.mannwhitneyu(sph[f'{k}_OE'].dropna(), chl[f'{k}_OE'].dropna(), alternative='two-sided')
        st.append(dict(test=f'{lab}: Sphaeropleales vs Chlamydomonadales (Mann-Whitney U)', statistic=u.statistic, P=u.pvalue,
                       n1=len(sph), n2=len(chl), median1=sph[f'{k}_OE'].median(), median2=chl[f'{k}_OE'].median()))
        w = stats.wilcoxon(S[f'{k}_OE'].dropna() - 1)
        st.append(dict(test=f'{lab}: all genomes vs 1 (Wilcoxon signed-rank)', statistic=w.statistic, P=w.pvalue,
                       n1=S[f'{k}_OE'].notna().sum(), median1=S[f'{k}_OE'].median()))
        r = stats.spearmanr(S[f'{k}_OE'], S.upstream_nc_GC, nan_policy='omit')
        st.append(dict(test=f'{lab} vs upstream non-coding GC (Spearman)', statistic=r.statistic, P=r.pvalue, n1=len(S)))
    sel = S[S.family == 'Selenastraceae']; osph = S[(S.group == 'Sphaeropleales') & (S.family != 'Selenastraceae')]
    u = stats.mannwhitneyu(sel.oct_OE, osph.oct_OE, alternative='less')
    st.append(dict(test='TATAATAT O/E: Selenastraceae vs other Sphaeropleales (Mann-Whitney, one-sided less)', statistic=u.statistic,
                   P=u.pvalue, n1=len(sel), n2=len(osph), median1=sel.oct_OE.median(), median2=osph.oct_OE.median()))
    tab = [[int(sel.oct_observed.sum()), float(sel.oct_expected.sum())], [int(osph.oct_observed.sum()), float(osph.oct_expected.sum())]]
    # pooled observed vs expected per group: Poisson rate-ratio test via conditional binomial
    k1, e1, k2, e2 = tab[0][0], tab[0][1], tab[1][0], tab[1][1]
    bt = stats.binomtest(k1, k1 + k2, e1 / (e1 + e2), alternative='less')
    st.append(dict(test='Pooled TATAATAT O/E Selenastraceae vs other Sphaeropleales (conditional binomial rate-ratio test)',
                   statistic=(k1/e1)/(k2/e2), P=bt.pvalue, n1=k1, n2=k2, median1=k1/e1, median2=k2/e2))
    o2 = osph[osph.genome != 'KM462873']  # sensitivity: Neocystis brevis falls outside Chlorophyceae in the ML tree
    k2b, e2b = int(o2.oct_observed.sum()), float(o2.oct_expected.sum())
    bt2 = stats.binomtest(k1, k1 + k2b, e1 / (e1 + e2b), alternative='less')
    st.append(dict(test='Pooled TATAATAT O/E Selenastraceae vs other Sphaeropleales excluding Neocystis (conditional binomial)',
                   statistic=(k1/e1)/(k2b/e2b), P=bt2.pvalue, n1=k1, n2=k2b, median1=k1/e1, median2=k2b/e2b))
    w = stats.wilcoxon(S.down_z_median, S.bg_down_z_median)
    st.append(dict(test="3' window MFE z (median) vs distal non-coding background (paired Wilcoxon)", statistic=w.statistic, P=w.pvalue,
                   n1=len(S), median1=S.down_z_median.median(), median2=S.bg_down_z_median.median()))
    w = stats.wilcoxon(S.start_z_median)
    st.append(dict(test='Start-codon region MFE z median vs 0 (Wilcoxon; positive = less structured than shuffled)',
                   statistic=w.statistic, P=w.pvalue, n1=len(S), median1=S.start_z_median.median()))
    # pooled within-genome: start z for CDS (A. falcatus)
    af = A[(A.genome == FOCAL) & (A.feature_type == 'CDS') & A.start_codon.isin(['ATG', 'GTG', 'TTG'])]
    w = stats.wilcoxon(af.start_region_z.dropna())
    st.append(dict(test='A. falcatus CDS start-region z vs 0 (Wilcoxon)', statistic=w.statistic, P=w.pvalue,
                   n1=af.start_region_z.notna().sum(), median1=af.start_region_z.median()))
    afb = BG[BG.genome == FOCAL]
    u = stats.mannwhitneyu(A[A.genome == FOCAL].down_z.dropna(), afb.z.dropna(), alternative='less')
    st.append(dict(test="A. falcatus 3' windows vs distal background windows, MFE z (Mann-Whitney, one-sided)",
                   statistic=u.statistic, P=u.pvalue, n1=len(A[A.genome == FOCAL]), n2=len(afb),
                   median1=A[A.genome == FOCAL].down_z.median(), median2=afb.z.median()))
    STATS = pd.DataFrame(st); STATS['q_BH'] = v2.bh_adjust(STATS.P.values)
    STATS.to_csv(T/'comparative_statistics.tsv', sep='\t', index=False)
    # Gene x genome conservation matrices (first context per gene)
    A1 = A.sort_values(['genome', 'gene', 'context_kind', 'anchor_id']).drop_duplicates(['genome', 'gene'])
    order = list(man.insdc_accession)
    mats = {}
    for col, fn in (('octamer', lambda d: (d.n_TATAATAT > 0).astype(float)),
                    ('TATAAT', lambda d: (d.n_TATAAT > 0).astype(float)),
                    ('SD_like', lambda d: d.SD_like.astype(float).where(d.feature_type == 'CDS')),
                    ('down_z', lambda d: d.down_z), ('start_z', lambda d: d.start_region_z),
                    ('down_q05', lambda d: (d.down_q_BH < .05).astype(float))):
        m = A1.assign(v=fn(A1)).pivot(index='gene', columns='genome', values='v').reindex(columns=order)
        m = m.loc[[g for g in CORE if g in m.index]]
        mats[col] = m
        m.to_csv(T/f'matrix_{col}_gene_x_genome.tsv', sep='\t')
    # gene-level conservation summary for A. falcatus genes
    sph_accs = list(S[S.group == 'Sphaeropleales'].genome)
    sel_accs = list(S[S.family == 'Selenastraceae'].genome)
    fa = A[A.genome == FOCAL].copy()
    rows = []
    for r in fa.itertuples():
        def frac(mat, accs):
            if r.gene not in mat.index:
                return np.nan
            v = mat.loc[r.gene, accs]
            return float(np.nanmean(v.values)) if v.notna().any() else np.nan
        rows.append(dict(anchor_id=r.anchor_id, gene=r.gene,
                         SD_cons_Sphaeropleales=frac(mats['SD_like'], sph_accs),
                         SD_cons_Selenastraceae=frac(mats['SD_like'], sel_accs),
                         down_q05_cons_Sphaeropleales=frac(mats['down_q05'], sph_accs),
                         down_q05_cons_Selenastraceae=frac(mats['down_q05'], sel_accs),
                         octamer_cons_Sphaeropleales=frac(mats['octamer'], sph_accs),
                         TATAAT_cons_Sphaeropleales=frac(mats['TATAAT'], sph_accs)))
    cons = pd.DataFrame(rows)
    fa = fa.merge(cons, on=['anchor_id', 'gene'])
    # Biotechnology-oriented ranking of endogenous A. falcatus parts (hypotheses)
    fa['validated_Cr_ortholog_part'] = fa.gene.map(VALIDATED_PART_GENES).fillna('')
    p_score = ((fa.upstream_noncoding_len >= 150).astype(int) +
               ((fa.n_TATAATAT + fa.n_TATAAT + fa.n_PEP_pairs) > 0).astype(int) +
               fa.SD_like.astype(int) + (fa.SD_cons_Sphaeropleales.fillna(0) >= .5).astype(int) +
               (fa.start_region_z.fillna(-9) >= 0).astype(int) + (fa.validated_Cr_ortholog_part != '').astype(int))
    t_score = ((fa.down_q_BH < .05).astype(int) + (fa.best_hairpin_energy.fillna(0) <= -10).astype(int) +
               (fa.down_q05_cons_Sphaeropleales.fillna(0) >= .5).astype(int) + (fa.downstream_noncoding_len >= 60).astype(int) +
               (fa.validated_Cr_ortholog_part != '').astype(int))
    fa['promoter_5UTR_part_score_0to6'] = p_score; fa['terminator_3UTR_part_score_0to5'] = t_score
    fa = fa.sort_values(['promoter_5UTR_part_score_0to6', 'terminator_3UTR_part_score_0to5'], ascending=False)
    fa.to_csv(T/'Afalcatus_integrated_candidates.tsv', sep='\t', index=False)
    print(S[['species', 'group', 'oct_OE', 'oct_P_enriched', 'sd_OE', 'sd_P_enriched', 'down_z_median', 'bg_down_z_median',
             'start_z_median', 'SD_fraction_CDS']].to_string())
    print(STATS.to_string())
    print(fa[['gene', 'upstream_noncoding_len', 'n_TATAATAT', 'n_TATAAT', 'best_PEP', 'SD_words', 'start_region_z', 'down_z',
              'down_q_BH', 'promoter_5UTR_part_score_0to6', 'terminator_3UTR_part_score_0to5']].head(25).to_string())


if __name__ == '__main__':
    main()
