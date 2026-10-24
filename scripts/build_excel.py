#!/usr/bin/env python3
"""Formatted Excel workbooks (Online Resources) from the TSV result tables."""
import argparse
import re
from pathlib import Path

import pandas as pd

FONT = 'Arial'


def write_book(path, sheets, readme):
    with pd.ExcelWriter(path, engine='xlsxwriter') as xw:
        wb = xw.book
        base = wb.add_format({'font_name': FONT, 'font_size': 9})
        head = wb.add_format({'font_name': FONT, 'font_size': 9, 'bold': True, 'bg_color': '#DCE6F2', 'border': 1,
                              'text_wrap': True, 'valign': 'top'})
        title = wb.add_format({'font_name': FONT, 'font_size': 12, 'bold': True})
        wrap = wb.add_format({'font_name': FONT, 'font_size': 9, 'text_wrap': True, 'valign': 'top'})
        ws = wb.add_worksheet('README'); ws.set_column(0, 0, 28, base); ws.set_column(1, 1, 120, wrap)
        ws.write(0, 0, readme['title'], title)
        r = 2
        for k, v in readme['rows']:
            ws.write(r, 0, k, head if k.isupper() else base); ws.write(r, 1, v, wrap); r += 1
        r += 1; ws.write(r, 0, 'SHEET', head); ws.write(r, 1, 'CONTENT', head); r += 1
        for name, (df, desc) in sheets.items():
            ws.write_url(r, 0, f"internal:'{name}'!A1", string=name); ws.write(r, 1, f'{desc} ({len(df):,} rows)', wrap); r += 1
        for name, (df, desc) in sheets.items():
            df = df.copy()
            for c in df.columns:
                if df[c].dtype == object:
                    df[c] = df[c].astype(str).str.slice(0, 32000).replace('nan', '')
            df.to_excel(xw, sheet_name=name, index=False, startrow=1)
            sh = xw.sheets[name]
            sh.write(0, 0, desc, wb.add_format({'font_name': FONT, 'italic': True, 'font_size': 9}))
            for j, c in enumerate(df.columns):
                sh.write(1, j, str(c), head)
                w = min(max(len(str(c)), *(df[c].astype(str).str.len().head(200).tolist() or [8])) + 2, 45)
                sh.set_column(j, j, w, base)
            sh.freeze_panes(2, 1)
            if len(df):
                sh.autofilter(1, 0, len(df) + 1, len(df.columns) - 1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..')
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'; O = P/'results/excel'; O.mkdir(parents=True, exist_ok=True)
    rd = lambda f, **k: pd.read_csv(T/f, sep='\t', **k)
    man = rd('genome_manifest.tsv')
    common = [('Workflow', 'Comparative chloroplast regulatory-element workflow v3 (scripts/run_all.sh). Seed 20261023; ViennaRNA 2.7.2 Turner 2004 parameters at 25 °C.'),
              ('Coordinates', 'Positions relative to the first transcribed base of the annotated product (start codon for CDS); negative = upstream. Genomic coordinates 1-based inclusive.'),
              ('Evidence level', 'All candidates are computational predictions. Only rows in Validated_elements are experimentally supported (primary literature, C. reinhardtii and other microalgae).'),
              ('Statistics', 'Null model: randomised Euler-path dinucleotide-preserving shuffles (approximate). Empirical P = (1+#null≥obs)/(B+1). BH = Benjamini–Hochberg within the stated family.')]
    sheets = {
        'Genome_manifest': (man, 'Genomes analysed, INSDC accessions retrieved from ENA, checksums, annotation source and Table S1 audit'),
        'Architecture': (rd('plastome_architecture.tsv'), 'Plastome size, GC, inverted repeat, intergenic fraction, introns'),
        'Genome_signal_summary': (rd('genome_level_signal_summary.tsv'), 'Genome-level observed/expected motif counts, structure z-scores and background controls'),
        'Statistics': (rd('comparative_statistics.tsv'), 'Comparative statistical tests (BH across this table)'),
        'Benchmark': (rd('benchmark_Creinhardtii.tsv'), 'Benchmark against experimentally validated C. reinhardtii elements'),
        'Benchmark_summary': (rd('benchmark_summary.tsv'), 'Recall by element class (exact binomial 95% CI)'),
        'Validated_elements': (pd.read_csv(P/'data/validated_elements_curated.tsv', sep='\t'), 'Curated experimentally validated chloroplast cis-elements (primary literature, verified PMIDs/DOIs)'),
        'Validated_projection': (rd('validated_element_projection.tsv'), 'Search for validated C. reinhardtii element sequences in orthologous upstream regions of all genomes'),
        'Window_sensitivity': (rd('promoter_window_sensitivity.tsv'), 'Promoter-word enrichment for 150/300/500-nt windows, with and without coding masking'),
        'Footprints_QC': (rd('footprint_motifs_QC.tsv'), 'Per-gene MEME footprints with positional QC, ORF-extension risk and TOMTOM best match'),
        'Footprint_sites': (rd('footprint_sites.tsv.gz'), 'Individual MEME sites per genome'),
        'Motif_discovery': (rd('motif_discovery.tsv'), 'All MEME (per gene) and STREME (pooled, by clade) motifs'),
        'TOMTOM': (rd('motif_tomtom_vs_validated.tsv'), 'TOMTOM comparisons of discovered motifs with validated reference motifs'),
        'Kmer8_meta': (rd('kmer8_cross_genome_enrichment.tsv'), 'Cross-genome 8-mer enrichment in −100..−1 non-coding upstream DNA (top 500)'),
        'Gene_pairs': (rd('gene_pair_coorientation_conservation.tsv'), 'Conservation of A. falcatus co-oriented adjacent gene pairs'),
        'Clusters': (rd('conserved_cooriented_clusters.tsv'), 'Conserved co-oriented gene clusters (adjacencies conserved in ≥50% of genomes)'),
        'Octamer_profiles': (rd('octamer_positional_profiles.tsv'), 'TATAATAT positional profiles (25-nt bins) vs shuffled expectation'),
        'Background_windows': (rd('distal_background_windows.tsv'), 'Distal non-coding 120-nt background windows (MFE z)'),
    }
    for m in ('octamer', 'TATAAT', 'SD_like', 'down_z', 'down_q05', 'start_z'):
        sheets[f'Matrix_{m}'] = (rd(f'matrix_{m}_gene_x_genome.tsv'), f'Gene × genome matrix: {m}')
    ph = P/'results/phylogeny'
    if (ph/'gene_alignments.tsv').exists():
        sheets['Phylogeny_genes'] = (pd.read_csv(ph/'gene_alignments.tsv', sep='\t'), 'Protein alignments used for the concatenated ML phylogeny')
    fm = pd.concat([pd.read_csv(f, sep='\t') for f in sorted((P/'figures').glob('figure_manifest_*.tsv'))]) if list((P/'figures').glob('figure_manifest_*.tsv')) else pd.DataFrame()
    sheets['Figure_manifest'] = (fm, 'Figures (JPEG 600 dpi + vector SVG) and captions')
    write_book(O/'ESM_2_Comparative_regulatory_genomics_results.xlsx', sheets,
               dict(title='Online Resource 2 — Comparative chloroplast regulatory genomics of 41 green-algal plastomes', rows=common))
    A = rd('all_boundary_contexts.tsv.gz'); H = rd('all_motif_hits.tsv.gz'); HP = rd('all_best_hairpins.tsv.gz')
    F = 'QB176266.1'
    sheets2 = {
        'Integrated_candidates': (rd('Afalcatus_integrated_candidates.tsv').drop(columns=['upstream300', 'down120'], errors='ignore'),
                                  'A. falcatus boundary contexts with all evidence, conservation and part scores (hypotheses)'),
        'Boundary_contexts': (A[A.genome == F], 'All A. falcatus CDS/rRNA boundary contexts with sequences and metrics'),
        'Motif_hits': (H[H.genome == F], 'All promoter-word, PEP-pair and SD-like hits with exact coordinates'),
        'Hairpins_all': (HP[HP.genome == F], 'All enumerated perfect-stem 3′ hairpins (ΔG ≤ −2 kcal/mol)'),
        'Footprint_sites_Af': (rd('footprint_sites.tsv.gz').query('genome == @F'), 'A. falcatus sites of conserved MEME footprints'),
        'Kmer8_Afalcatus': (rd('kmer8_Afalcatus_top.tsv'), 'Top enriched 8-mers in A. falcatus proximal upstream DNA'),
    }
    write_book(O/'ESM_3_Afalcatus_regulatory_inventory.xlsx', sheets2,
               dict(title='Online Resource 3 — Ankistrodesmus falcatus UNAP_012026 (QB176266.1) chloroplast regulatory inventory', rows=common))
    sheets3 = {'All_contexts': (A.drop(columns=['upstream300', 'down120', 'down_MFE_structure'], errors='ignore'), 'All boundary contexts of 41 genomes (sequences in TSV)'),
               'All_motif_hits': (H, 'All motif hits, 41 genomes'), 'Best_hairpins': (HP, 'Best 3′ hairpin per context (all hairpins for focal/benchmark genomes)')}
    write_book(O/'ESM_4_All_genomes_contexts_and_hits.xlsx', sheets3,
               dict(title='Online Resource 4 — Per-context results for all genomes', rows=common))


if __name__ == '__main__':
    main()
