#!/usr/bin/env python3
"""Supplementary Materials 2-4 for Plant Molecular Biology.

Each workbook opens with a 'Contents' sheet; every table is on its own sheet
named 'Table S<n>' with the PMB-style caption in row 1 (bold label, no final period).
"""
import argparse
from pathlib import Path

import pandas as pd

FONT = 'Arial'
TITLE = ('Computational characterization of transcriptional and translational regulatory elements in the '
         'Ankistrodesmus falcatus plastome reveals a lineage-specific promoter signature, relevant to chloroplast engineering')
NOTES = [
    ('Journal', 'Plant Molecular Biology'),
    ('Coordinates', 'Positions are relative to the start codon (CDS) or the first base of the mature product (rRNA); '
                    'negative = upstream. Genomic coordinates are 1-based, inclusive'),
    ('Evidence', 'All A. falcatus elements are computational predictions; only Table S3 lists experimentally validated elements'),
    ('Statistics', 'Null model: dinucleotide-preserving shuffles; empirical P = (1 + #null >= obs)/(B + 1); '
                   'BH = Benjamini-Hochberg correction within the stated family'),
    ('Software', 'v3 workflow, Supplementary Material 5 (seed 20261023; ViennaRNA 2.7.2, Turner 2004 parameters, 25 °C)'),
]


def write_book(path, sm_label, sm_desc, tables):
    with pd.ExcelWriter(path, engine='xlsxwriter') as xw:
        wb = xw.book
        base = wb.add_format({'font_name': FONT, 'font_size': 9})
        bold = wb.add_format({'font_name': FONT, 'font_size': 9, 'bold': True})
        head = wb.add_format({'font_name': FONT, 'font_size': 9, 'bold': True, 'bottom': 1, 'top': 1,
                              'text_wrap': True, 'valign': 'top'})
        big = wb.add_format({'font_name': FONT, 'font_size': 12, 'bold': True})
        wrap = wb.add_format({'font_name': FONT, 'font_size': 9, 'text_wrap': True, 'valign': 'top'})
        link = wb.add_format({'font_name': FONT, 'font_size': 9, 'bold': True, 'font_color': '#1F4E79', 'underline': 1})
        ws = wb.add_worksheet('Contents')
        ws.set_column(0, 0, 16, base); ws.set_column(1, 1, 130, wrap)
        ws.write(0, 0, sm_label, big)
        ws.write(1, 0, TITLE, wrap)
        ws.write(2, 0, sm_desc, base)
        r = 4
        for k, v in NOTES:
            ws.write(r, 0, k, bold); ws.write(r, 1, v, wrap); r += 1
        r += 1
        ws.write(r, 0, 'Table', head); ws.write(r, 1, 'Caption', head); r += 1
        for name, (df, cap) in tables.items():
            ws.write_url(r, 0, f"internal:'{name}'!A1", link, string=name)
            ws.write(r, 1, f'{cap} ({len(df):,} rows)', wrap); r += 1
        for name, (df, cap) in tables.items():
            df = df.copy()
            for c in df.columns:
                if df[c].dtype == object:
                    df[c] = df[c].astype(str).str.slice(0, 32000).replace('nan', '')
            df.to_excel(xw, sheet_name=name, index=False, startrow=2)
            sh = xw.sheets[name]
            sh.write_rich_string(0, 0, bold, name + ' ', base, cap)
            for j, c in enumerate(df.columns):
                sh.write(2, j, str(c), head)
                w = min(max(len(str(c)), *(df[c].astype(str).str.len().head(200).tolist() or [8])) + 2, 45)
                sh.set_column(j, j, w, base)
            sh.freeze_panes(3, 1)
            if len(df):
                sh.autofilter(2, 0, len(df) + 2, len(df.columns) - 1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--project', default='..'); ap.add_argument('--out', default=None)
    a = ap.parse_args(); P = Path(a.project); T = P/'results/tables'
    O = Path(a.out) if a.out else P/'results/excel_PMB'; O.mkdir(parents=True, exist_ok=True)
    rd = lambda f, **k: pd.read_csv(T/f, sep='\t', **k)
    mats = []
    for m in ('octamer', 'TATAAT', 'SD_like', 'down_z', 'down_q05', 'start_z'):
        x = rd(f'matrix_{m}_gene_x_genome.tsv').melt(id_vars='gene', var_name='genome', value_name='value')
        x.insert(2, 'metric', m); mats.append(x)
    F = 'QB176266.1'
    ph = P/'results/phylogeny/gene_alignments.tsv'
    sm2 = [
        (rd('genome_manifest.tsv'), 'Plastomes analysed: INSDC accessions retrieved from ENA, SHA-256 checksums, annotation source and audit of the original panel table'),
        (rd('plastome_architecture.tsv'), 'Plastome architecture: size, GC content, inverted repeat, intergenic fraction and introns'),
        (pd.read_csv(P/'data/validated_elements_curated.tsv', sep='\t'), 'Experimentally validated chloroplast cis-elements curated from primary literature, with verified PMIDs and DOIs'),
        (rd('benchmark_Creinhardtii.tsv'), 'Benchmark of the workflow against experimentally validated C. reinhardtii elements'),
        (rd('benchmark_summary.tsv'), 'Benchmark recall by element class, with exact binomial 95% confidence intervals'),
        (rd('validated_element_projection.tsv'), 'Search for validated C. reinhardtii element sequences in orthologous upstream regions of all plastomes'),
        (rd('genome_level_signal_summary.tsv'), 'Genome-level observed/expected motif counts, RNA-structure z-scores and background controls'),
        (rd('comparative_statistics.tsv'), 'Comparative statistical tests (BH correction across this table)'),
        (rd('promoter_window_sensitivity.tsv'), 'Promoter-word enrichment for 150-, 300- and 500-nt windows, with and without coding masking'),
        (rd('octamer_positional_profiles.tsv'), 'TATAATAT positional profiles (25-nt bins) versus shuffled expectation'),
        (rd('distal_background_windows.tsv'), 'Distal non-coding 120-nt background windows used for RNA-structure comparisons (MFE z-scores)'),
        (pd.concat(mats, ignore_index=True), 'Gene × genome matrices (long format) of octamer, TATAAT and SD-like counts and 3′-flank and start-region structure metrics'),
        (rd('footprint_motifs_QC.tsv'), 'Per-gene MEME phylogenetic footprints with positional QC, ORF-extension risk and best TOMTOM match'),
        (rd('footprint_sites.tsv.gz'), 'Individual MEME footprint sites per plastome'),
        (rd('motif_discovery.tsv'), 'All MEME (per gene) and STREME (pooled, by clade) motifs'),
        (rd('motif_tomtom_vs_validated.tsv'), 'TOMTOM comparisons of discovered motifs with validated reference motifs'),
        (rd('kmer8_cross_genome_enrichment.tsv'), 'Cross-genome 8-mer enrichment in non-coding upstream DNA (−100 to −1; top 500)'),
        (rd('gene_pair_coorientation_conservation.tsv'), 'Conservation of A. falcatus co-oriented adjacent gene pairs across 40 plastomes'),
        (rd('conserved_cooriented_clusters.tsv'), 'Conserved co-oriented gene clusters (adjacencies conserved in ≥50% of plastomes)'),
        (pd.read_csv(ph, sep='\t') if ph.exists() else pd.DataFrame(), 'Protein alignments (64 genes) used for the concatenated maximum-likelihood phylogeny'),
    ]
    A = rd('all_boundary_contexts.tsv.gz'); H = rd('all_motif_hits.tsv.gz'); HP = rd('all_best_hairpins.tsv.gz')
    sm3 = [
        (rd('Afalcatus_integrated_candidates.tsv').drop(columns=['upstream300', 'down120'], errors='ignore'),
         'A. falcatus candidate regulatory regions with integrated evidence, conservation, part scores and exact coordinates'),
        (A[A.genome == F], 'All A. falcatus CDS and rRNA boundary contexts with sequences and metrics'),
        (H[H.genome == F], 'All A. falcatus promoter-word, PEP-like pair and SD-like hits with exact coordinates'),
        (HP[HP.genome == F], 'All enumerated perfect-stem 3′ hairpins in A. falcatus (ΔG ≤ −2 kcal mol⁻¹)'),
        (rd('footprint_sites.tsv.gz').query('genome == @F'), 'A. falcatus sites of the conserved MEME footprints'),
        (rd('kmer8_Afalcatus_top.tsv'), 'Most enriched 8-mers in A. falcatus proximal upstream DNA'),
    ]
    sm4 = [
        (A.drop(columns=['upstream300', 'down120', 'down_MFE_structure'], errors='ignore'), 'Boundary contexts of all 41 plastomes'),
        (H, 'Motif hits in all 41 plastomes'),
        (HP, 'Best 3′ hairpin per context (all hairpins for the focal and benchmark plastomes)'),
    ]
    n = 1; books = []
    for label, desc, fname, lst in (
            ('Supplementary Material 2', 'Tables S1–S20. Comparative results for 41 green-algal plastomes', 'Supplementary_Material_2_Tables_S1-S20.xlsx', sm2),
            ('Supplementary Material 3', 'Tables S21–S26. Ankistrodesmus falcatus UNAP_012026 (QB176266.1) regulatory inventory', 'Supplementary_Material_3_Tables_S21-S26.xlsx', sm3),
            ('Supplementary Material 4', 'Tables S27–S29. Per-context results for all plastomes', 'Supplementary_Material_4_Tables_S27-S29.xlsx', sm4)):
        tabs = {}
        for df, cap in lst:
            tabs[f'Table S{n}'] = (df, cap); n += 1
        write_book(O/fname, label, desc, tabs); books.append((fname, tabs))
    with open(O/'supplementary_tables_index.tsv', 'w') as fh:
        fh.write('table\tsupplementary_material\tcaption\trows\n')
        for fname, tabs in books:
            for k, (df, cap) in tabs.items():
                fh.write(f'{k}\t{fname}\t{cap}\t{len(df)}\n')


if __name__ == '__main__':
    main()
