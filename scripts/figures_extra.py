#!/usr/bin/env python3
"""Fig 7 (footprints), Fig 8 (A. falcatus parts) and supplementary figures."""
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import logomaker
import numpy as np
import pandas as pd
import RNA
from Bio import Phylo

import figstyle as fs
from figstyle import plt
import matplotlib.patches as mpatches

FOCAL, CR = 'QB176266.1', 'FJ423446'


def meme_pwm(xml, consensus):
    root = ET.parse(xml).getroot()
    for m in root.iter('motif'):
        if m.get('name') == consensus:
            arr = m.find('probabilities').find('alphabet_matrix').findall('alphabet_array')
            return pd.DataFrame([{v.get('letter_id'): float(v.text) for v in a.findall('value')} for a in arr])[list('ACGT')]


def info_logo(ax, pwm, title):
    ic = logomaker.transform_matrix(pwm, from_type='probability', to_type='information')
    logomaker.Logo(ic, ax=ax, color_scheme='classic', font_name='Liberation Sans', vpad=0.02, width=0.9)
    ax.set_ylim(0, 2); ax.set_yticks([0, 1, 2]); ax.set_ylabel('bits', fontsize=5.5, labelpad=1)
    ax.set_xticks([]); ax.set_title(title, fontsize=5.8, loc='left', fontweight='normal')
    ax.spines['bottom'].set_visible(False)


def classify(cons):
    if 'GGAG' in cons or 'GAGG' in cons:
        return 'SD-like (GGAG family)', '#008300'
    if 'TATAAT' in cons or 'TTATAA' in cons:
        return 'TATAAT-containing', '#111111'
    if 'CGGA' in cons:
        return 'CGGA core', '#6250d6'
    return 'Other conserved block', '#2a78d6'


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; W = P/'results/motifs'
    FM, FS = P/'figures/main', P/'figures/supplementary'
    manifest = []
    Q = pd.read_csv(T/'footprint_motifs_QC.tsv', sep='\t')
    fp = Q[Q.positional_footprint].copy()
    fp['lo'] = fp.IQR_pos.str.split('..', regex=False).str[0].astype(float)
    # ---------------- Fig 7 ----------------
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.78))
    ax = fig.add_axes([0.07, 0.08, 0.40, 0.86])
    genes = sorted(fp.gene.unique(), key=lambda g: fp[fp.gene == g].median_pos.max())
    for i, g in enumerate(genes):
        ax.plot([-130, 0], [i, i], color=fs.GRID, lw=0.8, zorder=0)
        for r in fp[fp.gene == g].itertuples():
            lab, c = classify(r.consensus)
            x0 = r.median_pos
            ax.add_patch(mpatches.FancyBboxPatch((x0, i - 0.3), r.width, 0.6, boxstyle='round,pad=0,rounding_size=0.15',
                                                 color=c, lw=0, alpha=0.9))
            ax.text(x0 + r.width / 2, i + 0.42, f'{r.sites:.0f}/{r.n_sequences}', fontsize=3.8, ha='center', va='top', color=fs.MUTED)
        ax.text(-133, i, g, fontsize=5.5, ha='right', va='center', fontstyle='italic')
    ax.set_xlim(-130, 2); ax.set_ylim(len(genes) - 0.4, -0.8); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xlabel('Position of footprint relative to start codon (nt; median site)')
    handles = [mpatches.Patch(color=c, label=l) for l, c in (('SD-like (GGAG family)', '#008300'), ('TATAAT-containing', '#111111'),
                                                              ('CGGA core', '#6250d6'), ('Other conserved block', '#2a78d6'))]
    ax.legend(handles=handles, loc='lower left', fontsize=5, ncol=2, bbox_to_anchor=(0, -0.16))
    fs.panel(ax, 'a', -0.12, 1.0)
    ax.set_title('Positionally fixed 5′-leader footprints (MEME, ≥60% of genomes)', fontsize=6.3)
    picks = [('psbD', 'WTTWTTTTGGAGATCDMSBM', 'psbD −20: PRB1 GGAG region'),
             ('psbA', 'CGGAGAAWTWWAAATTTTW', 'psbA −33: CGGAG SD-like'),
             ('atpH', 'TTWHTAGGAGGAAAT', 'atpH −19: GGAGG'),
             ('rbcL', 'WWWYHCWCGGAYA', 'rbcL −55: CGGAC (cf. Cr element +38/+47)'),
             ('rbcL', 'GATTWTTATAGGATCGACAA', 'rbcL −42 block'),
             ('psbH', 'TATKACAGAAAGTAAAAAAA', 'psbH −40: AAGTAAA-like core (41/41)'),
             ('psbB', 'AAGHAATTAAGTAAAMAAAA', 'psbB −37: AAGTAAA-like core (40/40)'),
             ('psaB', 'GCAAAGTTCGCTTTGCAAAT', 'psaB −49: palindromic block')]
    for k, (g, cons, title) in enumerate(picks):
        r, c = divmod(k, 2)
        lx = fig.add_axes([0.56 + c * 0.235, 0.80 - r * 0.235, 0.195, 0.13])
        pw = meme_pwm(W/f'meme_{g}'/'meme.xml', cons)
        info_logo(lx, pw, title)
        if k == 0:
            fs.panel(lx, 'b', -0.25, 1.15)
    fs.save(fig, FM, 'Fig7_phylogenetic_footprints', manifest, 'Phylogenetic footprinting of 5′ leaders and logos of conserved elements')

    # ---------------- Fig 8: A. falcatus candidate parts ----------------
    C = pd.read_csv(T/'Afalcatus_integrated_candidates.tsv', sep='\t')
    C = C.drop_duplicates('gene')
    comp_p = {'nc ≥150 nt': C.upstream_noncoding_len >= 150, 'promoter motif': (C.n_TATAATAT + C.n_TATAAT + C.n_PEP_pairs) > 0,
              'SD-like': C.SD_like.astype(bool), 'SD conserved ≥50%': C.SD_cons_Sphaeropleales.fillna(0) >= .5,
              'start z ≥ 0': C.start_region_z.fillna(-9) >= 0, 'Cr part validated': C.validated_Cr_ortholog_part.fillna('') != ''}
    comp_t = {"3′ q<0.05": C.down_q_BH < .05, 'stem ΔG ≤ −10': C.best_hairpin_energy.fillna(0) <= -10,
              "3′ conserved ≥50%": C.down_q05_cons_Sphaeropleales.fillna(0) >= .5, '3′ nc ≥60 nt': C.downstream_noncoding_len >= 60,
              'Cr part validated ': C.validated_Cr_ortholog_part.fillna('') != ''}
    top = C.assign(tot=C['promoter_5UTR_part_score_0to6'] + C['terminator_3UTR_part_score_0to5']).sort_values(
        ['promoter_5UTR_part_score_0to6', 'tot'], ascending=False).head(18)
    M = pd.DataFrame({**{k: v for k, v in comp_p.items()}, **{k: v for k, v in comp_t.items()}}).loc[top.index].astype(float)
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.62))
    ax = fig.add_axes([0.09, 0.30, 0.42, 0.64])
    cm = np.where(np.arange(M.shape[1]) < len(comp_p), 0, 1)
    rgb = np.zeros((*M.shape, 3)) + 1
    for j in range(M.shape[1]):
        base = np.array([0.165, 0.47, 0.84]) if cm[j] == 0 else np.array([0.89, 0.29, 0.28])
        for i in range(M.shape[0]):
            rgb[i, j] = base if M.iloc[i, j] else np.array([0.94, 0.94, 0.93])
    ax.imshow(rgb, aspect='auto')
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene, fontstyle='italic', fontsize=5.5)
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=60, ha='right', fontsize=5.2)
    ax.axvline(len(comp_p) - 0.5, color='white', lw=2)
    for i, r in enumerate(top.itertuples()):
        ax.text(M.shape[1] - 0.3, i, f'  P {r.promoter_5UTR_part_score_0to6}/6   T {r.terminator_3UTR_part_score_0to5}/5',
                fontsize=5, va='center')
    ax.set_xticks(np.arange(-.5, M.shape[1]), minor=True); ax.set_yticks(np.arange(-.5, len(top)), minor=True)
    ax.grid(which='minor', color='white', lw=0.6); ax.tick_params(which='both', length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title('Evidence for endogenous A. falcatus parts (promoter/5′UTR | 3′UTR)', fontsize=6.3)
    fs.panel(ax, 'a', -0.2, 1.02)
    # cassette schematic using best promoter and best terminator
    bp = C.sort_values(['promoter_5UTR_part_score_0to6', 'terminator_3UTR_part_score_0to5'], ascending=False).iloc[0]
    bt = C.sort_values(['terminator_3UTR_part_score_0to5', 'down_z'], ascending=[False, True]).iloc[0]
    sx = fig.add_axes([0.63, 0.30, 0.36, 0.55]); sx.axis('off'); sx.set_xlim(0, 7.6); sx.set_ylim(0, 6)
    parts = [(0.2, 1.6, '#c9dcf3', f'P{bp.gene}\n5′UTR\n({bp.upstream_noncoding_len} nt)'),
             (1.9, 3.6, '#f6d28b', 'GOI\n(e.g. reporter)'), (5.6, 1.6, '#f4c1c0', f'{bt.gene}\n3′UTR\n(+1..+120)')]
    for x, w, c, lab in parts:
        sx.add_patch(mpatches.FancyBboxPatch((x, 2.4), w, 1.2, boxstyle='round,pad=0.02,rounding_size=0.2', color=c, lw=0))
        sx.text(x + w / 2, 3.0, lab, fontsize=5.2, ha='center', va='center', fontstyle='italic')
    sx.annotate('', xy=(1.4, 4.1), xytext=(0.4, 4.1), arrowprops=dict(arrowstyle='-|>', lw=0.8, color=fs.INK, mutation_scale=6))
    sx.plot([0.4, 0.4], [3.6, 4.1], color=fs.INK, lw=0.8)
    sx.text(0.2, 4.4, f'conserved leader footprints; SD {bp.SD_words or "–"}', fontsize=4.8)
    sx.text(5.6, 1.9, f'MFE z = {bt.down_z:.1f}; q = {bt.down_q_BH:.3f}; stem ΔG = {bt.best_hairpin_energy:.1f}', fontsize=4.8)
    sx.text(0.2, 0.8, 'Hypothesis for testing (not validated): endogenous parts\nranked by in silico evidence; insertion via flanking\nhomology arms at a neutral intergenic site.',
            fontsize=4.8, color=fs.MUTED)
    fs.panel(sx, 'b', -0.07, 0.855)
    sx.set_title('Candidate endogenous expression cassette', fontsize=6.3, loc='left', y=0.86)
    fs.save(fig, FM, 'Fig8_Afalcatus_candidate_parts', manifest, 'In silico ranking of endogenous A. falcatus regulatory parts and a candidate expression cassette')
    pd.DataFrame(manifest).to_csv((Path(fs.FIG_DIR) if fs.FIG_DIR else P/'figures')/'figure_manifest_main2.tsv', sep='\t', index=False)

    # ================= Supplementary =================
    sman = []
    S = pd.read_csv(T/'genome_level_signal_summary.tsv', sep='\t')
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t'); sp = dict(zip(man.insdc_accession, man.species))
    tree = Phylo.read(P/'results/phylogeny/plastome_ML.treefile', 'newick'); tree.ladderize(reverse=True)
    accmap = {x.split('.')[0]: x for x in man.insdc_accession}
    order = [accmap[t.name.rsplit('_', 1)[1]] for t in tree.get_terminals()]
    grp = dict(zip(S.genome, S.group))
    # S1 architecture
    fig, axs = plt.subplots(1, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.32))
    for gname, c in fs.GROUP_COLOR.items():
        d = S[S.group == gname]
        axs[0].scatter(d.length_bp / 1e3, d.intergenic_fraction * 100, s=9, color=c, lw=0, label=gname)
        axs[1].scatter(d.length_bp / 1e3, d.IR_length_bp / 1e3, s=9, color=c, lw=0)
        axs[2].scatter(d.GC, d.intergenic_GC, s=9, color=c, lw=0)
    for ax in axs:
        r = S[S.genome == FOCAL].iloc[0]
    af = S[S.genome == FOCAL].iloc[0]
    for ax, (x, y) in zip(axs, ((af.length_bp/1e3, af.intergenic_fraction*100), (af.length_bp/1e3, af.IR_length_bp/1e3), (af.GC, af.intergenic_GC))):
        ax.scatter([x], [y], s=22, facecolor='none', edgecolor=fs.FOCAL, lw=1); ax.annotate('Af', (x, y), xytext=(4, 2), textcoords='offset points', fontsize=5.5, fontweight='bold')
    axs[0].set_xlabel('Plastome size (kb)'); axs[0].set_ylabel('Intergenic DNA (%)'); axs[0].legend(fontsize=4.8)
    axs[1].set_xlabel('Plastome size (kb)'); axs[1].set_ylabel('Inverted repeat (kb)')
    axs[2].set_xlabel('Genome GC (%)'); axs[2].set_ylabel('Intergenic GC (%)')
    for i, ax in enumerate(axs):
        fs.panel(ax, 'abc'[i], -0.25)
    fig.tight_layout()
    fs.save(fig, FS, 'FigS1_plastome_architecture', sman, 'Plastome architecture of the 41 genomes')
    # S2 gene content
    G = pd.read_csv(T/'gene_content_matrix.tsv', sep='\t', index_col=0)
    keep = [c for c in G.columns if (G[c] > 0).mean() >= 0.15 and not c.startswith('trn')]
    G = (G.loc[order, keep] > 0).astype(int)
    G = G[G.mean().sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(fs.DOUBLE, fs.DOUBLE * 0.62))
    ax.imshow(G.values, aspect='auto', cmap='Greys', vmin=0, vmax=1.6, interpolation='nearest')
    ax.set_yticks(range(len(order))); ax.set_yticklabels([sp[g] for g in order], fontsize=4.2, fontstyle='italic')
    ax.set_xticks(range(len(G.columns))); ax.set_xticklabels(G.columns, rotation=90, fontsize=3.8, fontstyle='italic')
    ax.tick_params(length=0)
    fs.save(fig, FS, 'FigS2_gene_content', sman, 'Presence/absence of plastid protein and rRNA genes (genomes in phylogenetic order)')
    # S3 window sensitivity
    SE = pd.read_csv(T/'promoter_window_sensitivity.tsv', sep='\t')
    fig, axs = plt.subplots(1, 4, figsize=(fs.DOUBLE, fs.DOUBLE * 0.3), sharey=False)
    for ax, m in zip(axs, ('TATAATAT', 'TATAAT', 'TTGACA', 'PEP_pair')):
        for k, (L, mask) in enumerate([(w, mk) for w in (150, 300, 500) for mk in ('noncoding', 'unmasked')]):
            d = SE[(SE.motif == m) & (SE.window == L) & (SE.masking == mask)].merge(S[['genome', 'group']], on='genome')
            for gname, c in fs.GROUP_COLOR.items():
                v = d[d.group == gname].OE
                ax.scatter(np.full(len(v), k) + np.random.default_rng(k).uniform(-.15, .15, len(v)), v, s=4, color=c, lw=0)
            ax.plot([k - .25, k + .25], [d.OE.median()] * 2, color=fs.INK, lw=1)
            f = d[d.genome == FOCAL]
            if len(f):
                ax.scatter([k], f.OE, s=14, marker='D', facecolor='none', edgecolor=fs.FOCAL, lw=0.8)
        ax.axhline(1, ls='--', lw=0.5, color=fs.MUTED); ax.set_title(m, fontsize=6.3)
        ax.set_xticks(range(6)); ax.set_xticklabels([f'{w}\n{"nc" if mk=="noncoding" else "all"}' for w in (150, 300, 500) for mk in ('noncoding', 'unmasked')], fontsize=4.8)
        ax.set_ylabel('O/E' if m == 'TATAATAT' else '')
        if m == 'TATAATAT':
            ax.set_yscale('symlog', linthresh=1); ax.set_ylim(-0.2, 40); ax.set_yticks([0, 1, 10]); ax.set_yticklabels(['0', '1', '10'])
    fig.tight_layout()
    fs.save(fig, FS, 'FigS3_promoter_window_sensitivity', sman, 'Sensitivity of genome-level promoter-motif enrichment to window length and coding masking')
    # S4 octamer positional profiles small multiples
    PR = pd.read_csv(T/'octamer_positional_profiles.tsv', sep='\t')
    fig, axs = plt.subplots(7, 6, figsize=(fs.DOUBLE, fs.DOUBLE * 1.05), sharex=True)
    for ax, g in zip(axs.flat, order):
        d = PR[PR.genome == g]
        ax.bar(d.bin_start + 12.5, d.observed, width=22, color=fs.GROUP_COLOR[grp[g]] if g != FOCAL else fs.FOCAL, lw=0)
        ax.plot(d.bin_start + 12.5, d.expected, color=fs.INK, lw=0.6)
        ax.set_title(sp[g], fontsize=4.2, fontstyle='italic', loc='center', fontweight='normal')
        ax.tick_params(labelsize=3.8)
    for ax in list(axs.flat)[len(order):]:
        ax.axis('off')
    fig.supxlabel('Bin relative to first transcribed base (nt)', fontsize=6); fig.supylabel('TATAATAT count (bars) vs expected (line)', fontsize=6)
    fig.tight_layout(h_pad=0.4, w_pad=0.3)
    fs.save(fig, FS, 'FigS4_octamer_profiles_per_genome', sman, 'Per-genome TATAATAT positional profiles against dinucleotide-shuffled expectation')
    # S5 RNA structures of top A. falcatus 3' candidates
    A = pd.read_csv(T/'all_boundary_contexts.tsv.gz', sep='\t')
    HP = pd.read_csv(T/'all_best_hairpins.tsv.gz', sep='\t')
    af = A[A.genome == FOCAL].drop_duplicates('gene').sort_values('down_z').head(8)
    md = RNA.md(); md.temperature = 25
    fig, axs = plt.subplots(2, 4, figsize=(fs.DOUBLE, fs.DOUBLE * 0.55))
    for ax, r in zip(axs.flat, af.itertuples()):
        seq = r.down120.replace('T', 'U'); st, e = RNA.fold_compound(seq, md).mfe()
        coords = RNA.get_xy_coordinates(st)
        xy = np.array([(coords.get(i).X, coords.get(i).Y) for i in range(len(seq))])
        pt = RNA.ptable(st)
        for i in range(1, len(seq) + 1):
            if pt[i] > i:
                ax.plot(*xy[[i-1, pt[i]-1]].T, color='#e34948', lw=0.4)
        ax.plot(*xy.T, color=fs.INK, lw=0.5)
        ax.scatter(*xy.T, s=1.5, color=fs.INK, lw=0)
        ax.scatter(*xy[0], s=8, color='#2a78d6', lw=0)
        ax.set_aspect('equal'); ax.axis('off')
        ax.set_title(f'{r.gene} 3′ (+1..+120)\nMFE {e:.1f} kcal/mol, z = {r.down_z:.1f}', fontsize=5.2, fontstyle='italic', loc='center', fontweight='normal')
    fs.save(fig, FS, 'FigS5_Afalcatus_3prime_structures', sman, 'Predicted MFE structures (ViennaRNA, 25 °C) of the eight most structured A. falcatus 3′ flanks')
    # S6 gene-pair co-orientation conservation
    GP = pd.read_csv(T/'gene_pair_coorientation_conservation.tsv', sep='\t').sort_values('frac_conserved', ascending=True)
    fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 1.9))
    ax.barh(range(len(GP)), GP.frac_conserved * 100, color=fs.SPH, height=0.7, lw=0)
    ax.set_yticks(range(len(GP))); ax.set_yticklabels([f'{u}→{d} ({g} bp)' for u, d, g in zip(GP.upstream_gene, GP.downstream_gene, GP.focal_gap_bp)], fontsize=3.6, fontstyle='italic')
    ax.set_xlabel('Genomes conserving the co-oriented adjacency (%)'); ax.axvline(50, ls='--', lw=0.5, color=fs.MUTED)
    fs.save(fig, FS, 'FigS6_cooriented_gene_pairs', sman, 'Conservation of A. falcatus co-oriented adjacent gene pairs across 40 plastomes')
    # S7 intraspecific replicate concordance
    a1 = A[A.genome == 'DQ396875'].drop_duplicates('gene').set_index('gene'); a2 = A[A.genome == 'KX756229'].drop_duplicates('gene').set_index('gene')
    gs = a1.index.intersection(a2.index)
    fig, axs = plt.subplots(1, 2, figsize=(fs.SINGLE * 1.6, fs.SINGLE * 0.75))
    axs[0].scatter(a1.loc[gs, 'down_z'], a2.loc[gs, 'down_z'], s=6, color=fs.SPH, lw=0)
    rr = np.corrcoef(a1.loc[gs, 'down_z'], a2.loc[gs, 'down_z'])[0, 1]
    axs[0].plot([-12, 3], [-12, 3], ls='--', lw=0.5, color=fs.MUTED); axs[0].set_title(f"3′ MFE z (r = {rr:.2f})", fontsize=6)
    axs[0].set_xlabel('T. obliquus UTEX 393 (DQ396875)'); axs[0].set_ylabel('T. obliquus KX756229')
    axs[1].scatter(a1.loc[gs, 'start_region_z'], a2.loc[gs, 'start_region_z'], s=6, color=fs.SPH, lw=0)
    rr2 = np.corrcoef(a1.loc[gs, 'start_region_z'].fillna(0), a2.loc[gs, 'start_region_z'].fillna(0))[0, 1]
    axs[1].plot([-5, 4], [-5, 4], ls='--', lw=0.5, color=fs.MUTED); axs[1].set_title(f'Start-region z (r = {rr2:.2f})', fontsize=6)
    axs[1].set_xlabel('UTEX 393'); fig.tight_layout()
    fs.save(fig, FS, 'FigS7_intraspecific_concordance', sman, 'Concordance of signals between two Tetradesmus obliquus plastome assemblies')
    # S8 k-mer meta
    KM = pd.read_csv(T/'kmer8_cross_genome_enrichment.tsv', sep='\t').head(30)
    fig, ax = plt.subplots(figsize=(fs.SINGLE, fs.SINGLE * 1.2))
    ax.barh(range(len(KM))[::-1], KM.n_genomes_enriched, color=fs.SPH, height=0.7, lw=0)
    ax.set_yticks(range(len(KM))[::-1]); ax.set_yticklabels([f'{k} (pooled O/E {o:.1f})' for k, o in zip(KM.kmer, KM.pooled_OE)], fontsize=4.6, family='monospace')
    ax.set_xlabel('Genomes with O/E ≥ 2 and ≥3 occurrences')
    fs.save(fig, FS, 'FigS8_kmer_meta_enrichment', sman, 'Most consistently enriched 8-mers in proximal (−100..−1) non-coding upstream DNA')
    # S9 logos of all positional footprints
    n = len(fp); ncol = 4; nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(fs.DOUBLE, 0.62 * nrow))
    for ax, r in zip(axs.flat, fp.sort_values(['gene', 'median_pos']).itertuples()):
        info_logo(ax, meme_pwm(W/f'meme_{r.gene}'/'meme.xml', r.consensus),
                  f'{r.gene} {r.median_pos:.0f} ({r.sites:.0f}/{r.n_sequences}; E={r.E_value:.0e})')
    for ax in list(axs.flat)[n:]:
        ax.axis('off')
    fig.tight_layout(h_pad=0.5)
    fs.save(fig, FS, 'FigS9_all_footprint_logos', sman, 'Sequence logos of all positionally fixed 5′-leader footprints')
    pd.DataFrame(sman).to_csv((Path(fs.FIG_DIR) if fs.FIG_DIR else P/'figures')/'figure_manifest_supplementary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
