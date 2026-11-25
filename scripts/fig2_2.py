#!/usr/bin/env python3
"""Fig 2.2: ML phylogeny + architecture + genome-level regulatory signatures, with the
GenBank accession in parentheses after each species, and serif typography matching Fig 1.1."""
import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
from Bio import Phylo

import figstyle as fs
from figstyle import plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.offsetbox import TextArea, HPacker, AnnotationBbox
from figures import FOCAL, tree_layout


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; FM = P/'figures/main'
    import matplotlib.font_manager as _fm
    for _p in ['/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf',
               '/usr/share/fonts/truetype/liberation/LiberationSerif-BoldItalic.ttf']:
        try: _fm.fontManager.addfont(_p)
        except Exception: pass
    SERIF = {'font.family': 'serif', 'font.serif': ['Liberation Serif', 'Times New Roman', 'DejaVu Serif'], 'svg.fonttype': 'none'}
    S = pd.read_csv(T/'genome_level_signal_summary.tsv', sep='\t')
    man = pd.read_csv(T/'genome_manifest.tsv', sep='\t')
    sp = dict(zip(man.insdc_accession, man.species))
    grp = dict(zip(S.genome, S.group))
    tree = Phylo.read(P/'results/phylogeny/plastome_ML.treefile', 'newick')
    accmap = {x.split('.')[0]: x for x in man.insdc_accession}
    acc_of = lambda name: name.rsplit('_', 1)[1] if name else name
    tip_acc = lambda name: accmap[acc_of(name)]
    ACC_VER = {'QB176266': 'QB176266.1', 'KT199254': 'KT199254.1', 'KT625421': 'KT625421.1', 'KT199253': 'KT199253.1', 'DQ630521': 'DQ630521.1', 'KT199251': 'KT199251.1', 'KY114065': 'KY114065.1', 'MF276976': 'MF276976.1', 'MF276977': 'MF276977.1', 'KY114064': 'KY114064.1', 'MF276982': 'MF276982.1', 'MF276984': 'MF276984.1', 'MF276985': 'MF276985.1', 'MF276986': 'MF276986.1', 'KT199249': 'KT199249.1', 'KT625411': 'KT625411.1', 'KT199252': 'KT199252.1', 'KT199248': 'KT199248.1', 'KT199256': 'KT199256.1', 'KM462873': 'KM462873.1', 'MF423120': 'MF423120.1', 'CM031383': 'CM031383.1', 'ON619891': 'ON619891.1', 'OP124049': 'OP124049.1', 'PV295633': 'PV295633.1', 'KU847995': 'KU847995.1', 'OR350844': 'OR350844.1', 'DQ396875': 'DQ396875.1', 'CM007919': 'CM007919.1', 'KX756229': 'KX756229.1', 'OR502664': 'OR502664.1', 'KT199250': 'KT199250.1', 'OM329943': 'OM329943.1', 'CM002678': 'CM002678.1', 'KT199255': 'KT199255.1', 'MG786420': 'MG786420.1', 'FJ423446': 'FJ423446.1', 'GQ250046': 'GQ250046.1', 'AP012494': 'AP012494.1', 'EU677193': 'EU677193.1', 'HQ700713': 'HQ700713.1', 'AB001684': 'AB001684.1'}
    disp_acc = lambda a: ACC_VER.get(a.split('.')[0], a)
    tips, ty, tdepth = tree_layout(tree)
    order = [tip_acc(t.name) for t in tips]
    col = lambda acc: fs.FOCAL if acc == FOCAL else fs.GROUP_COLOR[grp[acc]]

    arch = S.set_index('genome')
    cols = [('length_bp', 'Size (kb)', 1e-3, None), ('intergenic_fraction', 'Intergenic (%)', 100, None),
            ('IR_length_bp', 'IR (kb)', 1e-3, None), ('oct_OE', 'TATAATAT O/E', 1, 1),
            ('sd_OE', 'SD-like O/E', 1, 1), ('down_z_median', "3′ MFE z", 1, 0)]

    with plt.rc_context(SERIF):
        fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 0.98))
        tx = fig.add_axes([0.0, 0.115, 0.46, 0.855])
        # ---- tree ----
        for c in tree.find_clades():
            for k in c.clades:
                tx.plot([tdepth[c], tdepth[c]], [ty[c], ty[k]], color=fs.INK, lw=0.6)
                tx.plot([tdepth[c], tdepth[k]], [ty[k], ty[k]], color=fs.INK, lw=0.6)
            pass
        # support values (aBayes) at every internal node, placed on the incoming branch
        xmax = max(tdepth.values())
        parent = {}
        for c in tree.find_clades():
            for k in c.clades:
                parent[k] = c
        for c in tree.find_clades():
            if c.is_terminal() or c == tree.root or not c.name:
                continue
            m = re.match(r'([\d.]+)/([\d.]+)', c.name)
            if not m:
                continue
            sh = float(m.group(1))
            tx.text(tdepth[c] - 0.0015 * xmax, ty[c] - 0.10, f'{sh:.0f}', fontsize=3.4, ha='right', va='bottom',
                    color=('#444444' if sh >= 80 else '#c0392b'),
                    bbox=dict(boxstyle='round,pad=0.05', fc='white', ec='none', alpha=0.7), zorder=5)
        xmax = max(tdepth.values())
        for t in tips:
            acc = tip_acc(t.name); focal = (acc == FOCAL)
            spname = TextArea(sp[acc], textprops=dict(fontstyle='italic', fontsize=5.0, color=col(acc),
                              fontweight='bold' if focal else 'normal', family='serif'))
            accname = TextArea(f'  ({disp_acc(acc)})', textprops=dict(fontsize=4.5, color=fs.MUTED, family='serif'))
            box = HPacker(children=[spname, accname], align='baseline', pad=0, sep=0)
            tx.add_artist(AnnotationBbox(box, (tdepth[t] + 0.01 * xmax, ty[t]), frameon=False,
                                         box_alignment=(0, 0.5), xycoords='data', pad=0))
        tx.set_ylim(len(tips) + 0.6, -0.9); tx.set_xlim(0, xmax * 2.25); tx.axis('off')
        tx.plot([0.01, 0.11], [len(tips)] * 2, color=fs.INK, lw=0.8)
        tx.text(0.06, len(tips) + 0.55, '0.1 substitutions/site', fontsize=5.2, ha='center', va='top')
        # ---- data columns ----
        w = 0.078
        for i, (c, title, mult, ref) in enumerate(cols):
            bx = fig.add_axes([0.49 + i * (w + 0.0095), 0.115, w, 0.855])
            vals = np.array([arch.loc[g, c] * mult for g in order], float)
            colors = [col(g) for g in order]
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
            bx.set_ylim(len(order) + 0.6, -0.9); bx.set_yticks([])
            bx.set_title(title, fontsize=5.8, loc='center', fontweight='bold')
            bx.spines['left'].set_visible(False); bx.tick_params(axis='x', labelsize=5.0)
            bx.grid(axis='x', color=fs.GRID, lw=0.3); bx.set_axisbelow(True)
        leg = [mpatches.Patch(color=v, label=k) for k, v in fs.GROUP_COLOR.items()] + \
              [mpatches.Patch(color=fs.FOCAL, label='A. falcatus (focal)'),
               Line2D([], [], marker='o', ls='', color=fs.INK, mfc='white', label='O/E, P ≥ 0.05'),
               Line2D([], [], marker='o', ls='', color=fs.INK, label='O/E, empirical P < 0.05'),
               Line2D([], [], marker='|', ls='', color=fs.MUTED, label='distal non-coding background z')]
        fig.legend(handles=leg, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.012), fontsize=5.8,
                   frameon=False, columnspacing=1.6, handletextpad=0.5)
        fs.save(fig, FM, 'Fig2.2_phylogeny_architecture_signals', [],
                'ML phylogeny of 41 plastomes with architecture and genome-level regulatory-signal statistics (Fig. 2.2)')
    print('Fig2.2 done')


if __name__ == '__main__':
    main()
