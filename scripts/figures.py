#!/usr/bin/env python3
"""Main and supplementary figures (600-dpi JPEG + vector SVG)."""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO

import figstyle as fs
from figstyle import plt

import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

FOCAL, CR = 'QB176266.1', 'FJ423446'
FUNC = [  # (regex on gene, class, colour)
    (r'^psb|^ycf12', 'Photosystem II', '#2a78d6'), (r'^psa|^ycf3|^ycf4', 'Photosystem I', '#1baf7a'),
    (r'^pet', 'Cytochrome b6f', '#eda100'), (r'^atp', 'ATP synthase', '#e87ba4'), (r'^rbcL', 'RuBisCO', '#008300'),
    (r'^rp[ls]|^infA|^tufA', 'Translation', '#6250d6'), (r'^rpo', 'RNA polymerase', '#eb6834'),
    (r'^rrn|^trn', 'rRNA / tRNA', '#e34948'), (r'.', 'Other / ORF', '#9a9893')]


def fclass(g):
    for pat, name, col in FUNC:
        if re.search(pat, g):
            return name, col


def tree_layout(tree):
    tree.ladderize(reverse=True)
    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}
    depth = tree.depths()
    def sety(c):
        if c.is_terminal():
            return y[c]
        y[c] = np.mean([sety(k) for k in c.clades]); return y[c]
    sety(tree.root)
    return tips, y, depth


def draw_tree(ax, tree, tips, y, depth, label_fn, color_fn, support=True, fontsize=5.5):
    for c in tree.find_clades():
        for k in c.clades:
            ax.plot([depth[c], depth[c]], [y[c], y[k]], color=fs.INK, lw=0.6)
            ax.plot([depth[c], depth[k]], [y[k], y[k]], color=fs.INK, lw=0.6)
        if support and not c.is_terminal() and c.confidence is not None and c != tree.root:
            pass
        if support and not c.is_terminal() and c.name:
            m = re.match(r'([\d.]+)/([\d.]+)', c.name)
            if m:
                sh, ab = float(m.group(1)), float(m.group(2))
                if sh >= 80 and ab >= 0.95:
                    ax.plot(depth[c], y[c], 'o', ms=2.0, color=fs.INK, mec='none')
                else:
                    ax.text(depth[c] - 0.004, y[c] - 0.15, f'{sh:.0f}/{ab:.2f}', fontsize=4.0, ha='right', va='bottom', color=fs.MUTED)
    xmax = max(depth.values())
    for t in tips:
        ax.text(depth[t] + 0.01 * xmax, y[t], label_fn(t.name), va='center', fontsize=fontsize, color=color_fn(t.name),
                fontstyle='italic', fontweight='bold' if t.name.endswith('QB176266') else 'normal')
    ax.set_ylim(len(tips) + 0.6, -0.5); ax.set_xlim(0, xmax * 1.75); ax.axis('off')
    ax.plot([0.01, 0.11], [len(tips)] * 2, color=fs.INK, lw=0.8)
    ax.text(0.06, len(tips) + 0.25, '0.1 substitutions/site', fontsize=5, ha='center', va='top')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'
    FM, FS = P/'figures/main', P/'figures/supplementary'
    manifest = []
    S = pd.read_csv(T/'genome_level_signal_summary.tsv', sep='\t')
    A = pd.read_csv(T/'all_boundary_contexts.tsv.gz', sep='\t')
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t')
    sp = dict(zip(man.insdc_accession, man.species))
    grp = dict(zip(S.genome, S.group))
    tree = Phylo.read(P/'results/phylogeny/plastome_ML.treefile', 'newick')
    acc_of = lambda name: name.rsplit('_', 1)[1] if name else name
    # map tip label suffix -> manifest accession
    accmap = {x.split('.')[0]: x for x in man.insdc_accession}
    tip_acc = lambda name: accmap[acc_of(name)]
    tips, ty, tdepth = tree_layout(tree)
    order = [tip_acc(t.name) for t in tips]
    lab = lambda name: sp[tip_acc(name)]
    col = lambda name: fs.GROUP_COLOR[grp[tip_acc(name)]] if tip_acc(name) != FOCAL else fs.FOCAL

    # ---------------- Fig 1: A. falcatus regulatory map ----------------
    rec = SeqIO.read(P/'data/standardized/QB176266.1.gb', 'genbank'); n = len(rec.seq)
    fa = A[A.genome == FOCAL]
    H = pd.read_csv(T/'all_motif_hits.tsv.gz', sep='\t'); Hf = H[H.genome == FOCAL]
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.78))
    ax = fig.add_axes([0.0, 0.02, 0.72, 0.96], projection='polar')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1); ax.axis('off'); ax.set_ylim(0, 10.6)
    th = lambda x: 2 * np.pi * x / n
    for f in rec.features:
        if f.type == 'misc_feature' or f.type == 'repeat_region':
            note = ' '.join(f.qualifiers.get('note', []) + f.qualifiers.get('rpt_type', []))
            if re.search(r'inverted|IR', note):
                ax.bar(th(int(f.location.start)), 10.6, width=th(len(f)), bottom=0, align='edge', color='#f1f0ec', lw=0, zorder=0)
    for f in rec.features:
        if f.type not in ('CDS', 'rRNA', 'tRNA'):
            continue
        g = f.qualifiers['gene'][0]; cname, c = fclass(g)
        for p in f.location.parts:
            r0 = 9.2 if p.strand == 1 else 8.3
            ax.bar(th(int(p.start)), 0.8, width=max(th(int(p.end) - int(p.start)), 0.003), bottom=r0, align='edge',
                   color=c, lw=0, zorder=2)
        if f.type != 'tRNA' and not g.startswith('orf'):
            mid = th(int(f.location.parts[0].start))
            r = 10.15 if f.location.parts[0].strand == 1 else 7.85
            ang = np.degrees(mid)
            rot = 90 - ang if ang <= 180 else 270 - ang
            ax.text(mid, r, g, fontsize=3.6, rotation=rot, rotation_mode='anchor',
                    ha='left' if (ang <= 180) == (r > 9) else 'right', va='center', fontstyle='italic', color=fs.INK)
    ax.plot(np.linspace(0, 2*np.pi, 400), [9.15]*400, color=fs.INK, lw=0.4)
    # motif tracks
    tracks = [('TATAAT', 6.4, '#2a78d6'), ('TATAATAT', 6.4, '#111111'), ('PEP_pair', 5.8, '#eb6834'), ('SD_like', 5.2, '#008300')]
    for motif, r, c in tracks:
        hh = Hf[Hf.motif == motif]
        if motif == 'SD_like':
            pos = fa.set_index('anchor_id').loc[hh.anchor_id].first_base_1b.values
        else:
            pos = hh.genomic_start_1b.values
        ax.vlines(th(np.asarray(pos, float)), r, r + (0.5 if motif != 'TATAATAT' else 0.75), color=c, lw=0.6 if motif != 'TATAATAT' else 1.1)
    for r in (6.4, 5.8, 5.2):
        ax.plot(np.linspace(0, 2*np.pi, 400), [r]*400, color=fs.GRID, lw=0.3)
    # 3' structure z bars
    zz = fa.dropna(subset=['down_z'])
    for rr in zz.itertuples():
        L = min(-rr.down_z, 12) * 0.17 if rr.down_z < 0 else 0
        ax.bar(th(rr.last_base_1b), L, width=0.009, bottom=4.6 - L, color='#e34948' if rr.down_q_BH < .05 else '#bdbbb5', lw=0, align='center')
    ax.plot(np.linspace(0, 2*np.pi, 400), [4.6]*400, color=fs.GRID, lw=0.3)
    ax.text(0, 0, f'Ankistrodesmus falcatus\nUNAP_012026\nQB176266.1\n{n:,} bp', ha='center', va='center', fontsize=6.5, fontstyle='italic')
    for r, tl in ((6.7, 'promoter words'), (5.95, 'PEP pairs'), (5.35, 'SD-like'), (4.1, "3′ structure")):
        ax.text(np.radians(-4), r, tl, fontsize=3.8, color=fs.MUTED, ha='right', va='center')
    for k in range(0, n, 20000):
        ax.text(th(k), 10.9, f'{k//1000} kb', fontsize=4, ha='center', va='center', color=fs.MUTED)
    lax = fig.add_axes([0.72, 0.08, 0.27, 0.85]); lax.axis('off')
    handles = [mpatches.Patch(color=c, label=nm) for _, nm, c in FUNC]
    handles += [mpatches.Patch(color='#f1f0ec', label='Inverted repeats (IRa/IRb)')]
    l1 = lax.legend(handles=handles, loc='upper left', title='Genes (outer: + strand; inner: − strand)', title_fontsize=6, fontsize=5.8)
    lax.add_artist(l1)
    h2 = [Line2D([], [], color='#111111', lw=1.3, label='TATAATAT (Chlamydomonas-type −10 octamer)'),
          Line2D([], [], color='#2a78d6', lw=1, label='TATAAT hexamer (upstream, non-coding)'),
          Line2D([], [], color='#eb6834', lw=1, label='PEP-like −35/−10 pair (≤1 mismatch, 16–18 nt)'),
          Line2D([], [], color='#008300', lw=1, label='SD-like word (3–15 nt before start)'),
          mpatches.Patch(color='#e34948', label="3′ window MFE z (bar length = −z), BH q < 0.05"),
          mpatches.Patch(color='#bdbbb5', label="3′ window MFE z, q ≥ 0.05")]
    lax.legend(handles=h2, loc='lower left', title='Regulatory-signal tracks (inner rings)', title_fontsize=6, fontsize=5.8)
    fs.save(fig, FM, 'Fig1_Afalcatus_regulatory_map', manifest,
            'Plastome map of A. falcatus with uniform regulatory-signal tracks')

    # ---------------- Fig 2: phylogeny + architecture + signal summary ----------------
    arch = S.set_index('genome')
    cols = [('length_bp', 'Plastome size (kb)', 1e-3, None), ('intergenic_fraction', 'Intergenic (%)', 100, None),
            ('IR_length_bp', 'IR length (kb)', 1e-3, None), ('oct_OE', 'TATAATAT O/E', 1, 1),
            ('sd_OE', 'SD-like O/E', 1, 1), ('down_z_median', "3′ MFE z (median)", 1, 0)]
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.95))
    tx = fig.add_axes([0.0, 0.04, 0.40, 0.93])
    draw_tree(tx, tree, tips, ty, tdepth, lab, col)
    w = 0.088
    for i, (c, title, mult, ref) in enumerate(cols):
        bx = fig.add_axes([0.41 + i * (w + 0.012), 0.04, w, 0.93])
        vals = np.array([arch.loc[g, c] * mult for g in order], float)
        colors = [fs.FOCAL if g == FOCAL else fs.GROUP_COLOR[grp[g]] for g in order]
        if c in ('oct_OE', 'sd_OE'):
            pcol = 'oct_P_enriched' if c == 'oct_OE' else 'sd_P_enriched'
            for j, (v, g) in enumerate(zip(vals, order)):
                sig = arch.loc[g, pcol] < .05
                bx.plot(v, j, 'o', ms=3.2, color=colors[j], mfc=colors[j] if sig else 'white', mew=0.7)
            bx.axvline(1, color=fs.MUTED, lw=0.5, ls='--')
        elif c == 'down_z_median':
            bz = [arch.loc[g, 'bg_down_z_median'] for g in order]
            for j, (v, b) in enumerate(zip(vals, bz)):
                bx.plot([b, v], [j, j], color=fs.GRID, lw=0.8)
                bx.plot(b, j, '|', ms=4, color=fs.MUTED, mew=0.8)
                bx.plot(v, j, 'o', ms=3, color=colors[j], mec='none')
            bx.axvline(0, color=fs.MUTED, lw=0.5, ls='--')
        else:
            bx.barh(range(len(order)), vals, color=colors, height=0.62, lw=0)
        bx.set_ylim(len(order) + 0.6, -0.5); bx.set_yticks([])
        bx.set_title(title, fontsize=5.8, loc='center', fontweight='bold')
        bx.spines['left'].set_visible(False); bx.tick_params(axis='x', labelsize=5)
        bx.grid(axis='x', color=fs.GRID, lw=0.3); bx.set_axisbelow(True)
    leg = [mpatches.Patch(color=v, label=k) for k, v in fs.GROUP_COLOR.items()] + [mpatches.Patch(color=fs.FOCAL, label='A. falcatus (focal)'),
           Line2D([], [], marker='o', ls='', color=fs.INK, mfc='white', label='O/E, P ≥ 0.05'),
           Line2D([], [], marker='o', ls='', color=fs.INK, label='O/E, empirical P < 0.05'),
           Line2D([], [], marker='|', ls='', color=fs.MUTED, label='distal non-coding background z')]
    fig.legend(handles=leg, loc='lower center', ncol=4, bbox_to_anchor=(0.6, -0.04), fontsize=5.5)
    fs.save(fig, FM, 'Fig2_phylogeny_architecture_signals', manifest,
            'ML phylogeny of 41 plastomes with architecture and genome-level regulatory-signal statistics')

    # ---------------- Fig 3: benchmark ----------------
    B = pd.read_csv(T/'benchmark_Creinhardtii.tsv', sep='\t')
    PR = pd.read_csv(T/'octamer_positional_profiles.tsv', sep='\t')
    fig, axs = plt.subplots(1, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.36), gridspec_kw=dict(width_ratios=[1.35, 1, 0.9]))
    ax = axs[0]
    lanes = [('rbcL', -92, [-104], [(-55, 10, 'stability element +38/+47')]),
             ('atpB', -340, [-352], [(-310, 12, 'stability element +31/+42')]),
             ('psbD', -74, [-224, -88], [(-13, 4, 'PRB1 GGAG')])]
    for i, (g, tss, octs, els) in enumerate(lanes):
        yv = i
        ax.plot([tss - 40, 0], [yv, yv], color=fs.GRID, lw=2.5, solid_capstyle='butt')
        ax.plot([tss, 0], [yv, yv], color='#9fc2ee', lw=2.5, solid_capstyle='butt')
        ax.annotate('', xy=(tss + 22, yv - 0.42), xytext=(tss, yv - 0.42),
                    arrowprops=dict(arrowstyle='-|>', lw=0.6, color=fs.INK, mutation_scale=5))
        ax.plot([tss, tss], [yv - 0.42, yv], color=fs.INK, lw=0.6)
        ax.text(tss + 24, yv - 0.42, 'TSS', fontsize=4.4, va='center')
        for o in octs:
            ok = tss - 15 <= o <= tss - 5
            ax.add_patch(mpatches.Rectangle((o, yv - 0.16), 8, 0.32, color='#111111' if ok else '#bdbbb5', lw=0))
            ax.text(o + 2, yv + 0.22, f'TATAATAT ({o - tss:+d})', fontsize=4.6, ha='right', va='top')
        for p, wdt, nm in els:
            ax.add_patch(mpatches.Rectangle((p, yv - 0.16), wdt, 0.32, color='#008300', lw=0))
            ax.text(p + wdt / 2, yv + 0.22, nm, fontsize=4.4, ha='center', va='top', color='#006400')
        ax.text(tss - 46, yv, g, fontsize=6.5, ha='right', va='center', fontstyle='italic')
        ax.text(2, yv, 'AUG', fontsize=4.6, ha='left', va='center')
    ax.set_xlim(-420, 30); ax.set_ylim(2.6, -0.75); ax.set_yticks([])
    ax.set_xlabel('Position relative to start codon (nt)'); ax.spines['left'].set_visible(False)
    fs.panel(ax, 'a', -0.05)
    ax.set_title('TSS-anchored −10 octamers (C. reinhardtii)', fontsize=6.5)
    ax = axs[1]
    for label, accs, c in (('Chlamydomonadales', list(S[S.group == 'Chlamydomonadales'].genome), fs.CHL),
                           ('Sphaeropleales', list(S[S.group == 'Sphaeropleales'].genome), fs.SPH),
                           ('A. falcatus', [FOCAL], fs.FOCAL)):
        g = PR[PR.genome.isin(accs)].groupby('bin_start')[['observed', 'expected']].sum()
        ratio = (g.observed + 0.5) / (g.expected + 0.5)
        ax.plot(g.index + 12.5, ratio, '-o', ms=2.2, lw=1, color=c, label=label)
    ax.axhline(1, color=fs.MUTED, lw=0.5, ls='--'); ax.set_xlabel('Bin centre relative to start codon / 5′ end (nt)')
    ax.set_ylabel('TATAATAT observed/expected'); ax.legend(fontsize=5.3, loc='upper left'); fs.panel(ax, 'b')
    ax.set_title('Positional octamer enrichment', fontsize=6.5)
    ax = axs[2]
    BS = pd.read_csv(T/'benchmark_summary.tsv', sep='\t')
    names = {'3prime_IR': "3′ IR (q<0.05)", 'promoter_-10_octamer': '−10 octamer', 'SD_like_required': 'Required SD-like',
             'promoter_-35/-10': '16S −35/−10'}
    BS = BS[BS.element_class.isin(names)]
    yy = np.arange(len(BS))
    ax.barh(yy, BS.n_detected_lenient / BS.n_validated, color='#c9dcf3', height=0.6, label='lenient criteria')
    ax.barh(yy, BS.recall, color=fs.SPH, height=0.6, label='strict scan criteria')
    for j, r in enumerate(BS.itertuples()):
        ax.text(1.02, j, f'{r.n_detected}/{r.n_validated} ({r.n_detected_lenient}/{r.n_validated})', fontsize=5, va='center')
    ax.set_yticks(yy); ax.set_yticklabels([names[x] for x in BS.element_class]); ax.set_xlim(0, 1.45)
    ax.set_xlabel('Recall of validated elements'); ax.legend(fontsize=5, loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=2); fs.panel(ax, 'c', -0.45)
    ax.set_title('Benchmark recall', fontsize=6.5)
    fig.tight_layout(w_pad=1.2)
    fs.save(fig, FM, 'Fig3_benchmark_validated_elements', manifest, 'Benchmark of the uniform scan against validated C. reinhardtii elements')

    # ---------------- Fig 4: genome-level motif enrichment ----------------
    fig, axs = plt.subplots(1, 4, figsize=(fs.DOUBLE, fs.DOUBLE * 0.33))
    groups = ['Sphaeropleales', 'Chlamydomonadales', 'Outgroups (OCC/Trebouxiophyceae)']
    short = {'Sphaeropleales': 'Sphaero-\npleales', 'Chlamydomonadales': 'Chlamydo-\nmonadales', 'Outgroups (OCC/Trebouxiophyceae)': 'Out-\ngroups'}
    rng = np.random.default_rng(1)
    for ax, (k, title) in zip(axs[:3], (('oct', 'TATAATAT'), ('pep', 'PEP −35/−10 pair'), ('sd', 'SD-like (CDS)'))):
        for i, gname in enumerate(groups):
            d = S[S.group == gname]
            x = i + rng.uniform(-0.18, 0.18, len(d))
            sig = (d[f'{k}_P_enriched'] < .05).values
            sel = (d.family == 'Selenastraceae').values
            for msk, mk, sz in ((~sel, 'o', 9), (sel, '^', 16)):
                ax.scatter(x[~sig & msk], d[f'{k}_OE'][~sig & msk], s=sz, marker=mk, facecolor='white', edgecolor=fs.GROUP_COLOR[gname], lw=0.7, zorder=3)
                ax.scatter(x[sig & msk], d[f'{k}_OE'][sig & msk], s=sz, marker=mk, color=fs.GROUP_COLOR[gname], lw=0, zorder=3)
            ax.plot([i - 0.28, i + 0.28], [d[f'{k}_OE'].median()] * 2, color=fs.INK, lw=1)
            for gg, lbl in ((FOCAL, 'Af'), (CR, 'Cr')):
                if gg in set(d.genome):
                    j = list(d.genome).index(gg)
                    ax.annotate(lbl, (x[j], d[f'{k}_OE'].iloc[j]), xytext=(4, 2), textcoords='offset points', fontsize=5.5, fontweight='bold')
        ax.axhline(1, color=fs.MUTED, lw=0.5, ls='--'); ax.set_xticks(range(3)); ax.set_xticklabels([short[g] for g in groups], fontsize=5.3)
        ax.set_ylabel('Observed / expected'); ax.set_title(title, fontsize=6.5)
    ax = axs[3]
    for gname in groups:
        d = S[S.group == gname]
        ax.scatter(d.upstream_nc_GC, d.oct_OE, s=9, color=fs.GROUP_COLOR[gname], lw=0)
    for gg, lbl in ((FOCAL, 'Af'), (CR, 'Cr')):
        r = S[S.genome == gg].iloc[0]; ax.annotate(lbl, (r.upstream_nc_GC, r.oct_OE), xytext=(3, 2), textcoords='offset points', fontsize=5.5, fontweight='bold')
    ax.set_xlabel('Upstream non-coding GC (%)'); ax.set_ylabel('TATAATAT O/E'); ax.axhline(1, color=fs.MUTED, lw=0.5, ls='--')
    ax.set_title('Composition control', fontsize=6.5)
    axs[0].legend(handles=[Line2D([], [], marker='^', ls='', mfc='white', mec=fs.SPH, label='Selenastraceae'),
                           Line2D([], [], marker='o', ls='', mfc=fs.INK, mec=fs.INK, label='P < 0.05'),
                           Line2D([], [], marker='o', ls='', mfc='white', mec=fs.INK, label='P ≥ 0.05')], fontsize=5, loc='upper right')
    for i, ax in enumerate(axs):
        fs.panel(ax, 'abcd'[i], -0.3)
    fig.tight_layout(w_pad=0.8)
    fs.save(fig, FM, 'Fig4_genome_level_motif_enrichment', manifest, 'Genome-level enrichment of promoter and SD-like words against dinucleotide-preserving nulls')

    # ---------------- Fig 5: RNA-structure signals ----------------
    BG = pd.read_csv(T/'distal_background_windows.tsv', sep='\t')
    fig, axs = plt.subplots(1, 3, figsize=(fs.DOUBLE, fs.DOUBLE * 0.34))
    ax = axs[0]
    for gname in groups:
        d = S[S.group == gname]
        for r in d.itertuples():
            ax.plot([0, 1], [r.bg_down_z_median, r.down_z_median], color=fs.GROUP_COLOR[gname], lw=0.5, alpha=0.7)
    for gg, lbl in ((FOCAL, 'Af'), (CR, 'Cr')):
        r = S[S.genome == gg].iloc[0]; ax.plot([0, 1], [r.bg_down_z_median, r.down_z_median], color=fs.FOCAL, lw=1.2)
        ax.text(1.04, r.down_z_median, lbl, fontsize=5.5, va='center', fontweight='bold')
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Distal non-coding\nbackground', "3′ flanks\n(+1..+120)"]); ax.set_xlim(-0.2, 1.25)
    ax.set_ylabel('Median MFE z-score (per genome)'); ax.set_title("3′ flank structure", fontsize=6.5)
    ax = axs[1]
    af = A[A.genome == FOCAL]; bgf = BG[BG.genome == FOCAL]
    bins = np.linspace(-12, 4, 33)
    ax.hist(bgf.z, bins=bins, density=True, color='#bdbbb5', label=f'distal background (n={len(bgf)})')
    ax.hist(af.down_z.dropna(), bins=bins, density=True, histtype='step', lw=1.2, color='#e34948', label=f"3′ flanks (n={af.down_z.notna().sum()})")
    ax.set_xlabel('MFE z-score'); ax.set_ylabel('Density'); ax.legend(fontsize=5.2, loc='upper left')
    ax.set_title('A. falcatus', fontsize=6.5, fontstyle='italic')
    ax = axs[2]
    for i, gname in enumerate(groups):
        d = S[S.group == gname]
        x = i + rng.uniform(-0.18, 0.18, len(d))
        ax.scatter(x, d.start_z_median, s=9, color=fs.GROUP_COLOR[gname], lw=0)
        ax.plot([i - 0.28, i + 0.28], [d.start_z_median.median()] * 2, color=fs.INK, lw=1)
    ax.axhline(0, color=fs.MUTED, lw=0.5, ls='--'); ax.set_xticks(range(3)); ax.set_xticklabels([short[g] for g in groups], fontsize=5.3)
    ax.set_ylabel('Median start-region (−30..+30) MFE z'); ax.set_title('Start-codon accessibility', fontsize=6.5)
    for i, ax in enumerate(axs):
        fs.panel(ax, 'abc'[i], -0.25)
    fig.tight_layout(w_pad=1)
    fs.save(fig, FM, 'Fig5_RNA_structure_signals', manifest, "3′-flank structure and start-codon accessibility relative to shuffled and genomic backgrounds")

    # ---------------- Fig 6: gene x genome conservation heatmaps ----------------
    msd = pd.read_csv(T/'matrix_SD_like_gene_x_genome.tsv', sep='\t', index_col=0)[order]
    mdz = pd.read_csv(T/'matrix_down_z_gene_x_genome.tsv', sep='\t', index_col=0)[order]
    keep = [g for g in mdz.index if mdz.loc[g].notna().mean() > .7]
    msd, mdz = msd.reindex(keep), mdz.loc[keep]
    fig, axs = plt.subplots(1, 2, figsize=(fs.DOUBLE, fs.DOUBLE * 0.95), gridspec_kw=dict(wspace=0.05))
    names = [sp[g] for g in order]
    for ax, M, cmap, vmin, vmax, title in ((axs[0], msd, 'Greens', 0, 1, 'SD-like word at 3–15 nt before start (presence)'),
                                           (axs[1], mdz, fs.DIV, -8, 8, "3′ flank MFE z-score")):
        im = ax.imshow(M.T.values, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_xticks(range(len(keep))); ax.set_xticklabels(keep, rotation=90, fontsize=4.3, fontstyle='italic')
        ax.set_title(title, fontsize=6.3, pad=8)
        if ax is axs[1]:
            cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.025, pad=0.07); cb.ax.tick_params(labelsize=5)
            cb.set_label('MFE z-score (negative = more structured than shuffles)', fontsize=5)
        else:
            ax.legend(handles=[mpatches.Patch(color=plt.get_cmap('Greens')(1.0), label='SD-like word present'),
                               mpatches.Patch(color=plt.get_cmap('Greens')(0.0), label='absent'),
                               mpatches.Patch(facecolor='white', edgecolor=fs.GRID, label='gene absent / not CDS')],
                      loc='upper center', bbox_to_anchor=(0.5, -0.07), ncol=3, fontsize=5)
        ax.set_yticks(range(len(order)))
        if ax is axs[0]:
            ax.set_yticklabels(names, fontsize=4.3, fontstyle='italic')
            for t, g in zip(ax.get_yticklabels(), order):
                t.set_color(fs.FOCAL if g == FOCAL else fs.GROUP_COLOR[grp[g]]); t.set_fontweight('bold' if g in (FOCAL, CR) else 'normal')
        else:
            ax.set_yticklabels([])
        ax.tick_params(length=0)
        ax.set_xticks(np.arange(-.5, len(keep)), minor=True); ax.set_yticks(np.arange(-.5, len(order)), minor=True)
        ax.grid(which='minor', color='white', lw=0.25); ax.tick_params(which='minor', length=0)
    fs.panel(axs[0], 'a', -0.32, 1.01); fs.panel(axs[1], 'b', -0.05, 1.03)
    fs.save(fig, FM, 'Fig6_gene_by_genome_conservation', manifest, 'Gene-by-genome conservation of SD-like words and 3′-flank structure, genomes in phylogenetic order')

    pd.DataFrame(manifest).to_csv((Path(fs.FIG_DIR) if fs.FIG_DIR else P/'figures')/'figure_manifest_main.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
