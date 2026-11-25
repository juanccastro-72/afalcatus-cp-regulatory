# RNA-seq transcriptional support pipeline (optional validation)

Standard (not strand-specific TSS) RNA-seq cannot pinpoint transcription start sites, but it
shows that candidate promoter/5′-UTR regions are transcribed and where coverage begins/ends.

Datasets (`datasets.tsv`): six public paired-end RNA-seq runs — two from the authors'
laboratory (same strain) and four from independent studies, incl. a Monoraphidium griffithi
cross-species control.

Steps (on any machine with internet; not runnable inside the restricted analysis sandbox):
1. `bash download_rnaseq.sh`                      # fetch FASTQs from ENA
2. `python rnaseq_plastome_coverage.py --genome ../data/genbank/QB176266.1.gb --reads reads/ --out cov/ [--max-reads 8000000]`
3. `python plot_rnaseq_validation.py --cov cov/ --genome ../data/genbank/QB176266.1.gb --candidates ../results/tables/Afalcatus_integrated_candidates.tsv`

Output: per-run plastome mapping rate, strand-specific coverage (.npz), and Fig. 10.
Send the small `cov/*.npz` files back to regenerate Fig. 10 and the Results text.
