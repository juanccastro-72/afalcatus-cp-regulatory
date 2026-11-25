# afalcatus-cp-regulatory

[![regression-tests](https://github.com/juanccastro-72/afalcatus-cp-regulatory/actions/workflows/tests.yml/badge.svg)](https://github.com/juanccastro-72/afalcatus-cp-regulatory/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- After the first Zenodo release, add its DOI badge here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Reproducible in silico identification of chloroplast regulatory elements in microalgal plastomes (v3 workflow; Supplementary Material 5 of the article).

This is the code for the article "Computational characterization of transcriptional and translational regulatory elements in the *Ankistrodesmus falcatus* plastome reveals a lineage-specific promoter signature, relevant to chloroplast engineering" (*Plant Molecular Biology*). It reproduces every table and figure from public sequence data with a fixed random seed (20261023). The workflow is summarised in Fig. S10.

## Quick start

    git clone https://github.com/juanccastro-72/afalcatus-cp-regulatory.git
    cd afalcatus-cp-regulatory
    conda env create -f environment.yml
    conda activate cpreg
    bash run_all.sh            # about 3 h on 2 CPUs; THREADS=8 bash run_all.sh to use more

External tools are found through the `CPREG_BIN` environment variable. If it is not set, the scripts use the directory that contains `mafft` on your `PATH`.

## Workflow and scripts

| Step | Script(s) | Output (Table/Fig.) |
|---|---|---|
| Regression tests | `test_inventory.py` (20 tests), `test_comparative.py` (16 tests) | – |
| Retrieve plastomes (ENA) | `fetch_ena.sh`, `ena_lookup.py` | `data/genbank/` |
| Standardise, QC and homology annotation; audit of the original panel table | `prepare_genomes.py`, `cpreg_genomes.py` | Table S1 |
| Boundary contexts (corrected v2 core) | `chloroplast_regulatory_inventory_improved.py`, `inventory_reporting.py` | – |
| Uniform signal scans and null models | `cpreg_scan.py`, `run_scans.py` | Tables S21–S29 |
| Benchmark against validated elements | `benchmark_validated.py` | Tables S4–S6 |
| Comparative statistics and candidate parts | `aggregate.py` | Tables S7, S8, S10–S12, S21 |
| Sensitivity analysis | `promoter_sensitivity.py` | Table S9, Fig. S3 |
| Phylogenetic footprinting and motifs | `footprinting.py`, `motif_qc.py` | Tables S13–S17, Figs. 7, S9 |
| Architecture and gene order | `architecture.py` | Tables S2, S18, S19, Figs. S1, S2, S6 |
| Phylogenomics | `phylogenomics.py` | Table S20, Fig. 2 |
| Figures (PMB format: `CPREG_FIG_STYLE=PMB`) | `figstyle.py`, `figures.py`, `figures_extra.py`, `figure_workflow.py` | Figs. 1–8, S1–S10 |
| Supplementary workbooks | `build_excel_PMB.py` (`build_excel.py` = earlier layout) | Supplementary Materials 2–4 |
| Manuscript tools: verified bibliography, Zotero-linked Word file, relinking | `manuscript_tools/` | – |
| Original v1 script comparison | `compare_original.py`, `rebuild_reports.py`, `verify_deliverables.py` | – |

## Data included

- `data/genbank/QB176266.1.gb`: the *A. falcatus* UNAP_012026 plastome.
- `data/genome_manifest_input.tsv`: the 42 INSDC accessions (the original panel plus references).
- `data/validated_elements_curated.tsv`: 67 experimentally validated chloroplast elements, with PMIDs and DOIs.
- `data/Table_S1.xlsx`: the original comparative panel table, before it was audited.

## Method summary

The workflow is knowledge-based and statistically controlled. It uses no machine-learning classifier.

- **Scans.** Promoter words, SD-like words, start-region folding and 3′-flank folding (ViennaRNA 2.7.2, Turner 2004 parameters, 25 °C).
- **Statistics.** Each window is compared with dinucleotide-preserving shuffles of itself; P values are BH-corrected and observed/expected ratios are computed per genome.
- **Comparative layer.** Phylogenetic footprinting (MEME/TOMTOM) and a phylogeny from 64 concatenated plastid proteins (IQ-TREE).

## Tests

The regression tests (`scripts/test_inventory.py`, `scripts/test_comparative.py`) run automatically on GitHub Actions for every change (`.github/workflows/tests.yml`, environment `ci-environment.yml`).

## How to cite

Please cite the article and the archived software release (see `CITATION.cff`; the Zenodo DOI is added at the first release):

Castro JC, Castro CG, Diaz KF, Acosta O, Cobos M (2026) Computational characterization of transcriptional and translational regulatory elements in the *Ankistrodesmus falcatus* plastome reveals a lineage-specific promoter signature, relevant to chloroplast engineering. Plant Molecular Biology (submitted).

## Licence

MIT (see `LICENSE`). External programs called by the workflow (ViennaRNA, MEME Suite, IQ-TREE, MAFFT, trimAl, BLAST+, ARAGORN) are not redistributed here and remain under their own licences.

## Contact

Juan C. Castro, Natural Resources Research Center (CIRNA), Universidad Nacional de la Amazonía Peruana (UNAP), Iquitos, Peru (juan.castro@unapiquitos.edu.pe); Marianela Cobos (marianela.cobos@unapiquitos.edu.pe)

## Revision update (2026-11)

This release accompanies the revised *Plant Molecular Biology* submission and adds:

- **Octamer ancestral-state analysis** now emits **internal-node reconstructions**
  (`results/tables/octamer_ancestral_internal_nodes.tsv`) and a **three-state sensitivity
  recoding** (`octamer_detection_states_sensitivity.tsv`), so non-significant enrichment is
  not equated with motif absence (`scripts/ancestral_octamer.py`).
- **Anti-SD / SD detection** with homology-based re-annotation of the 16S 3′ terminus and a
  thermodynamic (anti-SD:mRNA duplex) SD method (`scripts/sd_antisd_improved.py`,
  `scripts/fig_sd_antisd_improved.py`).
- **RNA-seq transcriptional support** pipeline (`rnaseq_pipeline/`): alignment of public
  non-strand-specific libraries to QB176266.1 and read-depth / 5′-read-start profiling.
  The data document transcription, not promoter activity.
- **Figures** are emitted as 600-dpi JPEG and **editable-text SVG** (`svg.fonttype=none`),
  compatible with CorelDRAW (`scripts/figstyle.py`).
- Supplementary workbooks updated to Tables S27–S34 (`scripts/build_excel_PMB.py`).
- Author list harmonised across the manuscript, Supplementary Material 1 and `CITATION.cff`.

Chloroplast transformation has been reported in another member of the Sphaeropleales
(*Asterarcys* sp., Scenedesmaceae; Kumari et al. 2022, Algal Research 66:102770); the
candidate parts here target the Selenastraceae, for which no transformation has yet been
reported. All *A. falcatus* elements are computational predictions awaiting experimental
validation.
