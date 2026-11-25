#!/usr/bin/env python3
"""Fig 1.1: improved A. falcatus plastome regulatory map (serif labels, de-overlapped
gene names, coloured LSC/SSC/IR regions, minimal centre text, no on-ring track labels)."""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

import figstyle as fs
from figstyle import plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from figures import FOCAL, FUNC, fclass

# region palette (soft, print-friendly; IRa and IRb share a colour)
REGION_RING = {'LSC': '#5B9BD5', 'SSC': '#70AD47', 'IRa': '#C0589E', 'IRb': '#C0589E'}
IR_WEDGE = '#F4ECF2'


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; FM = P/'figures/main'
    import matplotlib.font_manager as _fm
    for _p in ['/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf']:
        try: _fm.fontManager.addfont(_p)
        except Exception: pass
    SERIF = {'font.family': 'serif', 'font.serif': ['Liberation Serif', 'Times New Roman', 'DejaVu Serif'],
             'svg.fonttype': 'none'}
    rec = SeqIO.read(P/'data/standardized/QB176266.1.gb', 'genbank'); n = len(rec.seq)
    A = pd.read_csv(T/'all_boundary_contexts.tsv.gz', sep='\t'); fa = A[A.genome == FOCAL]
    H = pd.read_csv(T/'all_motif_hits.tsv.gz', sep='\t'); Hf = H[H.genome == FOCAL]
    ar = pd.read_csv(T/'plastome_architecture.tsv', sep='\t').query('genome == @FOCAL').iloc[0]
    irs = [tuple(int(x) for x in seg.strip().split('-')) for seg in str(ar.IR_coords).split('/')]
    (ira0, ira1), (irb0, irb1) = sorted(irs)
    regions = [('LSC', 1, ira0 - 1), ('IRa', ira0, ira1), ('SSC', ira1 + 1, irb0 - 1), ('IRb', irb0, irb1)]
    th = lambda x: 2 * np.pi * x / n
    manifest = []

    with plt.rc_context(SERIF):
        fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.95))
        ax = fig.add_axes([0.145, 0.275, 0.71, 0.71], projection='polar')
        ax.set_theta_zero_location('N'); ax.set_theta_direction(-1); ax.axis('off'); ax.set_ylim(0, 12.2)

        # translucent region bands spanning the full radius (centre -> rim), coloured by region
        for nm, a0, a1 in regions:
            ax.bar(th(a0 - 1), 11.45, width=th(a1 - a0 + 1), bottom=0, align='edge',
                   color=REGION_RING[nm], alpha=0.20, lw=0, zorder=0)

        # genes
        labels = []  # (angle, strand, gene, start, end)
        for f_ in rec.features:
            if f_.type not in ('CDS', 'rRNA', 'tRNA'):
                continue
            g = f_.qualifiers['gene'][0]; cname, c = fclass(g)
            for p in f_.location.parts:
                r0 = 9.2 if p.strand == 1 else 8.3
                ax.bar(th(int(p.start)), 0.8, width=max(th(int(p.end) - int(p.start)), 0.003), bottom=r0, align='edge',
                       color=c, lw=0, zorder=2)
            if f_.type != 'tRNA' and not g.startswith('orf'):
                st = f_.location.parts[0]
                labels.append((th(int(st.start)), st.strand, g, int(st.start), int(st.end)))
        ax.plot(np.linspace(0, 2*np.pi, 400), [9.15]*400, color=fs.INK, lw=0.4)

        # ---- de-overlap gene labels: two radial tiers per strand ----
        gap = np.radians(3.0)                       # minimum angular separation within a tier
        tiers = {1: (10.12, 10.78), -1: (7.88, 7.22)}  # (near, far) radius per strand
        for strand in (1, -1):
            grp = sorted([x for x in labels if x[1] == strand], key=lambda x: x[0])
            last = [-9, -9]                         # last placed angle per tier
            for ang, _s, g, s0, s1 in grp:
                t = 0 if (ang - last[0]) >= gap else (1 if (ang - last[1]) >= gap else
                        (0 if (ang - last[0]) >= (ang - last[1]) else 1))
                last[t] = ang
                r = tiers[strand][t]
                deg = np.degrees(ang)
                rot = 90 - deg if deg <= 180 else 270 - deg
                ha = 'left' if (deg <= 180) == (strand == 1) else 'right'
                if t == 1:                          # connector from gene ring to pushed-out label
                    r_in = 10.0 if strand == 1 else 8.0
                    ax.plot([ang, ang], [r_in, r - (0.05 if strand == 1 else -0.05)], color='#c8c7c2', lw=0.3, zorder=1)
                ax.text(ang, r, g, fontsize=4.0, rotation=rot, rotation_mode='anchor',
                        ha=ha, va='center', fontstyle='italic', color=fs.INK, zorder=4)

        # regulatory-signal tracks (no on-ring text labels)
        tracks = [('TATAAT', 6.4, '#2a78d6'), ('TATAATAT', 6.4, '#111111'), ('PEP_pair', 5.8, '#eb6834'), ('SD_like', 5.2, '#008300')]
        for motif, r, c in tracks:
            hh = Hf[Hf.motif == motif]
            pos = (fa.set_index('anchor_id').loc[hh.anchor_id].first_base_1b.values if motif == 'SD_like'
                   else hh.genomic_start_1b.values)
            ax.vlines(th(np.asarray(pos, float)), r, r + (0.5 if motif != 'TATAATAT' else 0.75),
                      color=c, lw=0.6 if motif != 'TATAATAT' else 1.1)
        for r in (6.4, 5.8, 5.2):
            ax.plot(np.linspace(0, 2*np.pi, 400), [r]*400, color=fs.GRID, lw=0.3)
        zz = fa.dropna(subset=['down_z'])
        for rr in zz.itertuples():
            L = min(-rr.down_z, 12) * 0.17 if rr.down_z < 0 else 0
            ax.bar(th(rr.last_base_1b), L, width=0.009, bottom=4.6 - L,
                   color='#e34948' if rr.down_q_BH < .05 else '#bdbbb5', lw=0, align='center')
        ax.plot(np.linspace(0, 2*np.pi, 400), [4.6]*400, color=fs.GRID, lw=0.3)

        # minimal centre: scientific name + size
        ax.text(0, 0.75, 'Ankistrodesmus falcatus', ha='center', va='center', fontsize=9.5, fontstyle='italic', color=fs.INK)
        ax.text(0, -0.55, f'{n:,} bp', ha='center', va='center', fontsize=7.5, color=fs.INK)

        # kb ruler
        for k in range(0, n, 10000):
            ax.vlines(th(k), 10.95, 11.15, color=fs.MUTED, lw=0.4)
            if k % 20000 == 0:
                ax.text(th(k), 10.66, f'{k//1000}', fontsize=4.6, ha='center', va='center', color=fs.MUTED)

        # coloured structural region rim (regions identified by colour; see key), 50% opacity, eccentric
        for nm, a0, a1 in regions:
            ax.bar(th(a0 - 1), 0.42, width=th(a1 - a0 + 1), bottom=11.58, align='edge',
                   color=REGION_RING[nm], alpha=0.50, lw=0.4, edgecolor='white', zorder=3)

        # bottom legends: Genes (left) and Regulatory-signal tracks (right)
        gax = fig.add_axes([0.015, 0.005, 0.47, 0.235]); gax.axis('off')
        rax = fig.add_axes([0.515, 0.005, 0.47, 0.235]); rax.axis('off')
        gh = [mpatches.Patch(color=c, label=nm) for _, nm, c in FUNC]
        gh += [mpatches.Patch(color=REGION_RING['LSC'], label='LSC (large single-copy)'),
               mpatches.Patch(color=REGION_RING['SSC'], label='SSC (small single-copy)'),
               mpatches.Patch(color=REGION_RING['IRa'], label='IRa / IRb (inverted repeats)')]
        gax.legend(handles=gh, loc='center', ncol=2, title='Genes (outer ring: + strand; inner ring: − strand) and regions',
                   title_fontsize=6.6, fontsize=6.0, handlelength=1.3, columnspacing=1.2, labelspacing=0.5, frameon=False)
        h2 = [Line2D([], [], color='#111111', lw=1.3, label='TATAATAT (Chlamydomonas-type −10 octamer)'),
              Line2D([], [], color='#2a78d6', lw=1, label='TATAAT hexamer (upstream, non-coding)'),
              Line2D([], [], color='#eb6834', lw=1, label='PEP-like −35/−10 pair (≤1 mismatch, 16–18 nt)'),
              Line2D([], [], color='#008300', lw=1, label='SD-like sequence (3–15 nt before start)'),
              mpatches.Patch(color='#e34948', label="3′ window MFE z (bar = −z), BH q < 0.05"),
              mpatches.Patch(color='#bdbbb5', label="3′ window MFE z, q ≥ 0.05")]
        rax.legend(handles=h2, loc='center', ncol=1,
                   title='Regulatory-signal tracks (inner rings, outer→inner:\npromoter motifs · PEP pairs · SD-like · 3′ structure)',
                   title_fontsize=6.2, fontsize=6.0, handlelength=1.6, labelspacing=0.5, frameon=False)
        fs.save(fig, FM, 'Fig1.1_Afalcatus_regulatory_map', manifest,
                'Improved plastome map of A. falcatus (Fig. 1.1)')
    print('Fig1.1 done')


if __name__ == '__main__':
    main()
