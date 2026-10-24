#!/usr/bin/env bash
# Reproduce the complete analysis (Supplementary Material 5).
# 1) conda env create -f environment.yml && conda activate cpreg
# 2) bash run_all.sh            (about 3 h on 2 CPUs; set THREADS to use more)
set -euo pipefail
cd "$(dirname "$0")/scripts"
THREADS=${THREADS:-2}
export CPREG_BIN=${CPREG_BIN:-$(dirname "$(command -v mafft)")}
python test_inventory.py && python test_comparative.py          # regression tests
bash fetch_ena.sh                                                # 40 INSDC plastomes from ENA
python prepare_genomes.py --project ..                           # standardise + QC + homology annotation
python architecture.py --project ..
python run_scans.py --project .. --threads "$THREADS" --nrep 200 --nrep-fold 100 --nrep-start 50 --temperature 25 --seed 20261023
python benchmark_validated.py --project ..
python aggregate.py --project ..
python promoter_sensitivity.py --project .. --threads "$THREADS"
python footprinting.py --project .. && python motif_qc.py --project ..
python phylogenomics.py --project .. --threads "$THREADS"
# Figures in Plant Molecular Biology format (Fig1..Fig8, FigS1..FigS10; JPEG 600 dpi, PDF, SVG; <=174 x 234 mm)
export CPREG_FIG_STYLE=PMB CPREG_FIG_DIR=../figures_PMB
python figures.py --project .. && python figures_extra.py --project .. && python figure_workflow.py --project ..
unset CPREG_FIG_STYLE CPREG_FIG_DIR
# Supplementary Materials 2-4 (sheets 'Table S1' ... 'Table S29')
python build_excel_PMB.py --project .. --out ../supplementary_materials
