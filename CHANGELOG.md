# Changelog

## v1.1.0 — 2026-11-11 (PMB revision)
- Added internal-node ancestral reconstructions and a three-state sensitivity recoding for
  the TATAATAT −10 octamer (ancestral_octamer.py).
- Added homology re-annotation of the 16S 3′ terminus and a thermodynamic SD-detection method
  (sd_antisd_improved.py, fig_sd_antisd_improved.py).
- Added the RNA-seq transcriptional-support pipeline (rnaseq_pipeline/) using a merged BAM of
  four A. falcatus libraries; coverage reported as read depth (libraries are non-strand-specific).
- Figures now exported as editable-text SVG (svg.fonttype=none) in addition to JPEG.
- Supplementary workbooks extended to Tables S27–S34 (build_excel_PMB.py); SM1 author list fixed.
- CITATION.cff: added sixth author (M. A. Grandez).

## v1.0.0 — 2026-10
- Initial reproducible workflow for in silico identification of chloroplast regulatory elements.
