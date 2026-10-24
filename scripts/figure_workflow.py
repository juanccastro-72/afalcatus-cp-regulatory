#!/usr/bin/env python3
"""Supplementary Fig. S10: framework of the v3 regulatory-element workflow (JPEG + SVG)."""
import argparse
from pathlib import Path

import figstyle as fs
from figstyle import plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

C = {'in': '#e8eef8', 'prep': '#eef3e6', 'scan': '#fdf0e2', 'stat': '#f3eaf6', 'cmp': '#e6f4f1', 'val': '#fbe9e9', 'out': '#f1f0ec'}
E = {'in': '#2a78d6', 'prep': '#008300', 'scan': '#eb6834', 'stat': '#6250d6', 'cmp': '#1baf7a', 'val': '#e34948', 'out': '#6b6a66'}


def box(ax, x, y, w, h, kind, title, lines, script=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.004,rounding_size=0.012',
                                fc=C[kind], ec=E[kind], lw=0.9))
    ax.text(x + 0.012, y + h - 0.012, title, fontsize=6.6, fontweight='bold', va='top', color=fs.INK)
    yy = y + h - 0.036
    for ln in lines:
        ax.text(x + 0.014, yy, '• ' + ln, fontsize=5.3, va='top', color=fs.INK)
        yy -= 0.0185
    if script:
        ax.text(x + w - 0.01, y + 0.008, script, fontsize=4.8, ha='right', va='bottom', color=E[kind], family='monospace')


def arrow(ax, a, b, style='-|>'):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=7, lw=0.8, color=fs.MUTED, shrinkA=0, shrinkB=0))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    P = Path(ap.parse_args().project)
    fig = plt.figure(figsize=(fs.DOUBLE, fs.DOUBLE * 1.07))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0.095, 1); ax.axis('off')

    # 1 inputs
    box(ax, 0.03, 0.865, 0.30, 0.11, 'in', '1  Inputs', [
        'A. falcatus QB176266.1 (GenBank)', '40 plastomes: Table S1 panel + references', '(INSDC records via ENA; SHA-256 checksums)'])
    box(ax, 0.35, 0.865, 0.30, 0.11, 'in', '1  Prior knowledge', [
        '67 validated cis-elements (57 studies)', 'C. reinhardtii and other microalgae', 'PMIDs/DOIs verified (Europe PMC)'])
    box(ax, 0.67, 0.865, 0.30, 0.11, 'in', '1  Parameters', [
        'upstream ≤300 nt; 3′ flank +1..+120', 'ViennaRNA 2.7.2, Turner 2004, 25 °C', 'seed 20261023 (fully reproducible)'])
    # 2 preparation
    box(ax, 0.03, 0.70, 0.45, 0.14, 'prep', '2  Genome preparation and QC', [
        'parse GenBank/EMBL (ENA CON records fixed)', 'gene-name normalisation (rplB→rpl2, 16S→rrn16 …)',
        'homology annotation of WGS scaffolds: tblastn/blastn + ARAGORN', 'translation QC; Table S1 audit; duplicate records removed'],
        'prepare_genomes.py · cpreg_genomes.py')
    box(ax, 0.52, 0.70, 0.45, 0.14, 'prep', '3  Boundary contexts (v2 core, corrected)', [
        'CDS/rRNA products → oriented 5′/3′ boundary contexts', 'mixed-strand trans-splicing split into exon contexts',
        'circular wrap-around; exact 1-based coordinate paths', 'contiguous non-coding masking (no exon on either strand)'],
        'chloroplast_regulatory_inventory_improved.py')
    # 4 scans
    ax.text(0.03, 0.675, '4  Uniform signal scans in every genome (cpreg_scan.py)', fontsize=7, fontweight='bold')
    sw = 0.228
    scans = [('Promoter words', ['TATAATAT −10 octamer (exact)', 'TATAAT, TTGACA (exact)', 'PEP pairs: ≤1 mismatch,', '16–18 nt spacer']),
             ('Shine–Dalgarno-like', ['GGAG-family words', 'spacer 3–15 nt before', 'start codon (1st exon)', 'nested words collapsed']),
             ('Start-codon region', ['−30..+30 nt folded', '(ViennaRNA MFE)', 'z vs 50 shuffles', '(accessibility)']),
             ("3′ flank structure", ['+1..+120 MFE, z vs 100', 'shuffles; perfect', 'antiparallel stems', '(4–12 bp, 3–10 nt loop)'])]
    for k, (t, ls) in enumerate(scans):
        box(ax, 0.03 + k * (sw + 0.013), 0.535, sw, 0.125, 'scan', t, ls)
    # 5 statistics
    box(ax, 0.03, 0.365, 0.45, 0.15, 'stat', '5  Null models and statistics', [
        'dinucleotide-preserving shuffles of each window (Euler path)', 'per-context empirical P; Benjamini–Hochberg q',
        'genome-level O/E from summed shuffle replicates', 'positional profiles (25-nt bins); distal non-coding background',
        'window-length (150/300/500 nt) and masking sensitivity'],
        'cpreg_scan.py · aggregate.py · promoter_sensitivity.py')
    box(ax, 0.52, 0.365, 0.45, 0.15, 'cmp', '6  Comparative and evolutionary layer', [
        'phylogenetic footprinting: MEME (ZOOPS) on orthologous leaders', 'positional QC, ORF-extension risk, TOMTOM vs validated motifs',
        'STREME by clade; 8-mer census', '64-protein ML phylogeny (MAFFT, trimAl, IQ-TREE)', 'IR, gene content, co-oriented gene-order conservation'],
        'footprinting.py · motif_qc.py · phylogenomics.py · architecture.py')
    # 7 validation
    box(ax, 0.03, 0.205, 0.45, 0.14, 'val', '7  Validation', [
        'benchmark on C. reinhardtii: TSS-anchored −10 octamers,', '  required SD-like elements, 16S promoter, 3′ IRs (recall, 95% CI)',
        'projection of validated sequences to all 41 plastomes', '36 regression tests; independent T. obliquus assemblies'],
        'benchmark_validated.py · test_*.py')
    box(ax, 0.52, 0.205, 0.45, 0.14, 'out', '8  Integration and outputs', [
        'per-gene evidence + conservation across Sphaeropleales', 'additive part scores: promoter/5′UTR (0–6), 3′UTR (0–5)',
        'Excel workbooks, BED/TSV tables, 600-dpi JPEG, PDF and SVG figures', 'manuscript tables and Supplementary Materials'],
        'aggregate.py · build_excel.py · figures*.py')
    # note box
    ax.add_patch(FancyBboxPatch((0.03, 0.105), 0.94, 0.08, boxstyle='round,pad=0.004,rounding_size=0.012', fc='white', ec=fs.INK, lw=0.8, ls='--'))
    ax.text(0.045, 0.178, 'Design principle', fontsize=6.8, fontweight='bold', va='top')
    note = ('Knowledge-based and statistically controlled; no machine-learning classifier (e.g. a convolutional neural network) is trained. '
            'Candidate elements are defined by explicit sequence rules derived from experimentally validated elements, thermodynamic RNA folding '
            '(ViennaRNA) and composition-preserving null models. The only model fitted to data is the MEME motif model (expectation maximisation), '
            'used for unsupervised discovery. A supervised classifier was not used because fewer than 15 experimentally mapped chloroplast elements '
            'can be positioned precisely, almost all from C. reinhardtii, which is too few to train or evaluate one without overfitting.')
    import textwrap
    ax.text(0.045, 0.158, '\n'.join(textwrap.wrap(note, 175)), fontsize=5.4, va='top', color=fs.INK, linespacing=1.35)
    # arrows
    for x in (0.18, 0.50, 0.82):
        arrow(ax, (x, 0.865), (x if x < 0.5 else 0.745 if x > 0.6 else 0.30, 0.84))
    arrow(ax, (0.48, 0.77), (0.52, 0.77))
    arrow(ax, (0.745, 0.70), (0.745, 0.665))
    for k in range(4):
        xc = 0.03 + k * (sw + 0.013) + sw / 2
        arrow(ax, (xc, 0.535), (0.255 if k < 2 else 0.745, 0.515))
    arrow(ax, (0.255, 0.365), (0.255, 0.345))
    arrow(ax, (0.745, 0.365), (0.745, 0.345))
    arrow(ax, (0.48, 0.275), (0.52, 0.275))
    arrow(ax, (0.48, 0.44), (0.52, 0.44))
    fs.save(fig, P/'figures/supplementary', 'FigS10_workflow_framework')


if __name__ == '__main__':
    main()
