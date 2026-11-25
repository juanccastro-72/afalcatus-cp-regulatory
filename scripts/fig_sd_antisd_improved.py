#!/usr/bin/env python3
"""Figure summarising the improved SD / anti-SD inventory."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import logomaker
from pathlib import Path

for p in ['/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf',
          '/usr/share/fonts/truetype/liberation/LiberationSerif-Italic.ttf',
          '/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf']:
    try: fm.fontManager.addfont(p)
    except Exception: pass
plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif','Times New Roman','DejaVu Serif'],
                     'svg.fonttype':'none','font.size':8})

T = Path('results/tables')
AD = pd.read_csv(T/'antiSD_16S_3prime_improved.tsv', sep='\t')
SDs = pd.read_csv(T/'SD_summary_by_genome_improved.tsv', sep='\t')
SDg = pd.read_csv(T/'SD_per_gene_improved.tsv', sep='\t')

fig = plt.figure(figsize=(174/25.4, 185/25.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.15], hspace=0.62, wspace=0.32)

# (a) anti-SD detection: current vs improved
ax = fig.add_subplot(gs[0, 0])
cur = 22; reann = int((AD['method'] == 'reannotated_homology').sum()); ann = int((AD['method'] == 'annotated').sum())
absent = int((AD['method'] == 'absent_in_assembly').sum())
ax.bar(['Current\n(annotated\ntail only)'], [cur], color='#9ab7d3', width=0.6)
ax.bar(['Improved\n(+ homology\nreannotation)'], [ann], color='#2a78d6', width=0.6, label='annotated')
ax.bar(['Improved\n(+ homology\nreannotation)'], [reann], bottom=[ann], color='#7fc97f', width=0.6, label='reannotated')
ax.set_ylabel('plastomes with detectable anti-SD'); ax.set_ylim(0, 41)
ax.axhline(41, ls=':', c='#999', lw=0.8); ax.text(1.5, 41.4, 'n = 41', fontsize=6, ha='right')
ax.set_title('a  Anti-SD recovery', loc='left', fontsize=9, fontweight='bold')
ax.legend(fontsize=6, loc='upper left', frameon=False)
ax.text(0, cur+0.6, str(cur), ha='center', fontsize=7)
ax.text(1, ann+reann+0.6, f'{ann+reann}', ha='center', fontsize=7)

# (b) homology-identity gap
ax = fig.add_subplot(gs[0, 1])
pres = AD[AD['present'] == 'yes']['hom_identity']; ab = AD[AD['present'] == 'absent']['hom_identity']
ax.plot(np.random.RandomState(1).normal(0, 0.05, len(pres)), pres, 'o', ms=3, color='#2a78d6', label='present')
ax.plot(np.random.RandomState(2).normal(1, 0.05, len(ab)), ab, 'o', ms=3, color='#d9534f', label='absent')
ax.axhline(0.73, ls='--', c='#444', lw=0.8); ax.text(1.45, 0.735, 'threshold 0.73', fontsize=6, ha='right')
ax.set_xticks([0, 1]); ax.set_xticklabels([f'present\n(n={len(pres)})', f'absent\n(n={len(ab)})'])
ax.set_ylabel('16S 3′ homology identity'); ax.set_ylim(0.45, 1.03)
ax.set_title('b  Clean present/absent gap', loc='left', fontsize=9, fontweight='bold')

# (c) anti-SD logo from recovered tails (aligned on CCTCC core)
ax = fig.add_subplot(gs[1, :])
LW, RW = 8, 3
wins = []
for t in AD[AD['present'] == 'yes']['tail']:
    t = str(t).upper().replace('U', 'T'); p = t.rfind('CCTCC')
    if p-LW >= 0 and p+5+RW <= len(t): wins.append(t[p-LW:p+5+RW])
L = LW+5+RW
counts = pd.DataFrame(0, index=range(L), columns=list('ACGT'))
for w in wins:
    for i, b in enumerate(w):
        if b in 'ACGT': counts.loc[i, b] += 1
info = logomaker.transform_matrix(counts, from_type='counts', to_type='information')
logomaker.Logo(info, ax=ax, color_scheme='classic')
ax.set_ylabel('bits'); ax.set_ylim(0, 2.1)
ax.set_xticks(range(L)); ax.set_xticklabels(list(('·'*LW+'CCTCC'+'·'*RW)[:L]), fontsize=6)
ax.set_title(f'c  16S 3′ anti-SD region, aligned on the CCTCC core (n = {len(wins)} plastomes, improved set)',
             loc='left', fontsize=9, fontweight='bold')

# (d) SD detection: three methods
ax = fig.add_subplot(gs[2, 0])
tot = len(SDg)
vals = [SDg.strict.sum(), SDg.expanded.sum(), SDg.thermo.sum()]
pct = [100*v/tot for v in vals]
bars = ax.bar(['strict\n(exact GGAG)', 'expanded\ncores', 'thermodynamic\n(ΔG ≤ −4)'], vals,
              color=['#9ab7d3', '#2a78d6', '#7b5fb0'], width=0.65)
ax.set_ylabel(f'CDS with SD-like element\n(of {tot})')
for b, v, p in zip(bars, vals, pct): ax.text(b.get_x()+b.get_width()/2, v+4, f'{v}\n{p:.1f}%', ha='center', fontsize=6)
ax.set_ylim(0, max(vals)*1.25)
ax.set_title('d  SD detection, three methods', loc='left', fontsize=9, fontweight='bold')

# (e) per-genome %CDS with SD, distributions
ax = fig.add_subplot(gs[2, 1])
data = [SDs.pct_strict, SDs.pct_expanded, SDs.pct_thermo]
bp = ax.boxplot(data, tick_labels=['strict', 'expanded', 'thermo'], widths=0.55, patch_artist=True,
                medianprops=dict(color='black'))
for patch, c in zip(bp['boxes'], ['#9ab7d3', '#2a78d6', '#7b5fb0']): patch.set_facecolor(c)
ax.set_ylabel('% CDS with SD per genome')
ax.set_title('e  Per-genome SD rate', loc='left', fontsize=9, fontweight='bold')

fig.savefig('results/figures/supplementary/FigS12_SD_antiSD_improved.jpg', dpi=600, bbox_inches='tight')
fig.savefig('results/figures/supplementary/FigS12_SD_antiSD_improved.pdf', bbox_inches='tight')
fig.savefig('results/figures/supplementary/FigS12_SD_antiSD_improved.svg', bbox_inches='tight')
print('wrote FigS12_SD_antiSD_improved (n present =', len(wins), ')')
