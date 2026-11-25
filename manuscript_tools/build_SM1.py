"""Supplementary Material 1 (Figs. S1-S10) as a PDF in Plant Molecular Biology style."""
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image, Spacer, PageBreak, KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PI
import glob

def font(name, pats):
    for p in pats:
        f = glob.glob(p, recursive=True)
        if f: pdfmetrics.registerFont(TTFont(name, f[0])); return name
    raise SystemExit('font missing ' + name)
R = font('TNR', ['/usr/share/fonts/**/LiberationSerif-Regular.ttf'])
B = font('TNR-B', ['/usr/share/fonts/**/LiberationSerif-Bold.ttf'])
I = font('TNR-I', ['/usr/share/fonts/**/LiberationSerif-Italic.ttf'])
BI = font('TNR-BI', ['/usr/share/fonts/**/LiberationSerif-BoldItalic.ttf'])
from reportlab.lib.fonts import addMapping
addMapping('TNR', 0, 0, R); addMapping('TNR', 1, 0, B); addMapping('TNR', 0, 1, I); addMapping('TNR', 1, 1, BI)

FIG = Path(__import__('sys').argv[1] if len(__import__('sys').argv) > 1 else '../figures_PMB/supplementary')
OUT = Path(__import__('sys').argv[2] if len(__import__('sys').argv) > 2 else 'Supplementary_Material_1_Figs_S1-S10.pdf')
TITLE = ('Computational characterization of transcriptional and translational regulatory elements in the '
         '<i>Ankistrodesmus falcatus</i> plastome reveals a lineage-specific promoter signature, relevant to chloroplast engineering')
CAP = {
 1: 'Plastome architecture of the 41 genomes. <b>a</b> Size versus intergenic fraction. <b>b</b> Size versus inverted-repeat length. <b>c</b> Genome versus intergenic GC content. Af, <i>A. falcatus</i>',
 2: 'Presence/absence of plastid protein-coding and rRNA genes (genomes in phylogenetic order)',
 3: 'Sensitivity of genome-level promoter-motif enrichment (observed/expected against 200 dinucleotide-preserving shuffles) to upstream window length (150, 300 and 500 nt) and to non-coding masking (nc) versus unmasked windows (all). Diamonds, <i>A. falcatus</i>; bars, medians',
 4: 'Per-genome TATAATAT positional profiles in 25-nt bins upstream of the start codon. Bars, observed counts; lines, expectation from dinucleotide-preserving shuffles',
 5: 'Predicted minimum-free-energy structures (ViennaRNA 2.7.2, 25 °C) of the eight most structured <i>A. falcatus</i> 3′ flanks (+1 to +120). The blue dot marks the first nucleotide after the stop codon',
 6: 'Conservation of <i>A. falcatus</i> co-oriented adjacent gene pairs (transcription order; intergenic gap in <i>A. falcatus</i> in parentheses) across the other 40 plastomes',
 7: 'Concordance of per-gene 3′-flank MFE z-scores and start-region z-scores between two independent <i>Tetradesmus obliquus</i> plastome assemblies',
 8: 'The 30 most consistently enriched 8-mers in proximal (−100 to −1) non-coding upstream DNA across genomes',
 9: 'Information-content logos of all 28 positionally fixed 5′-leader footprints (gene, median position, genomes with a site/genomes analysed and MEME <i>E</i> value)',
 11: 'Improved detection of the 16S rRNA 3′ anti-SD and of SD-like elements. <b>a</b> Anti-SD recovery: homology-based reannotation of the 16S 3′ terminus raises detection from 22/41 (annotated tail only) to 28/41 (annotated + reannotated). <b>b</b> 16S 3′ homology identity separates present (0.77–1.00) from absent (0.53–0.67) plastomes with a clean gap at 0.73; the 13 non-recovered genomes lack the conserved 3′ terminus in the deposited sequence (assembly-limited, scattered across lineages). <b>c</b> Information-content logo of the recovered 16S 3′ anti-SD region aligned on the CCUCC core (n = 26 plastomes). <b>d</b> SD-like elements detected across 3,288 CDS by three methods: strict exact GGAG (251; 7.6%), expanded cores (443; 13.5%) and a thermodynamic anti-SD:mRNA duplex (ΔG ≤ −4 kcal mol⁻¹; 338; 10.3%). <b>e</b> Per-genome percentage of CDS with an SD-like element under each method. The strict exact-word scan undercounts SD-like elements by 40–75%',
 10: 'Framework of the reproducible workflow (v3) for identifying chloroplast regulatory elements. Steps 1–8 show the inputs, genome preparation, boundary contexts, the four uniform signal scans, null models and statistics, the comparative layer, validation against experimentally characterised elements, and integration into ranked candidate parts; the scripts responsible for each step are given in each box (Supplementary Material 5). The workflow is knowledge-based and statistically controlled and does not train a machine-learning classifier such as a convolutional neural network',
}
st = dict(
 h=ParagraphStyle('h', fontName='TNR', fontSize=14, leading=18, spaceAfter=8),
 t=ParagraphStyle('t', fontName='TNR', fontSize=11, leading=15, spaceAfter=10),
 n=ParagraphStyle('n', fontName='TNR', fontSize=10.5, leading=14, spaceAfter=6),
 c=ParagraphStyle('c', fontName='TNR', fontSize=10, leading=13, spaceBefore=6))
W, H = 174 * mm, 230 * mm
def header(c, d):
    c.saveState(); c.setFont(R, 8.5)
    c.drawString(18 * mm, A4[1] - 12 * mm, 'Plant Molecular Biology — Supplementary Material 1')
    c.drawRightString(A4[0] - 18 * mm, 12 * mm, str(d.page)); c.restoreState()
doc = SimpleDocTemplate(str(OUT), pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=18 * mm, bottomMargin=18 * mm,
                        title='Supplementary Material 1', author='Castro JC, Castro CG, Diaz KF, Acosta O, Grandez MA, Cobos M')
s = [Paragraph('<b>Supplementary Material 1</b>', st['h']),
     Paragraph(TITLE, st['t']),
     Paragraph('Juan C. Castro, Carlos G. Castro, Katia F. Diaz, Oscar Acosta, Miguel A. Grandez, Marianela Cobos*', st['n']),
     Paragraph('*Corresponding author: juan.castro@unapiquitos.edu.pe; Universidad Nacional de la Amazonía Peruana, Iquitos, Peru', st['n']),
     Paragraph('<i>Plant Molecular Biology</i>', st['n']), Spacer(1, 6 * mm),
     Paragraph('<b>Contents</b>', st['n'])]
SHORT = {11: 'Improved SD / anti-SD detection', 1: 'Plastome architecture', 2: 'Gene content', 3: 'Sensitivity of promoter-motif enrichment',
         4: 'TATAATAT positional profiles per genome', 5: '<i>A. falcatus</i> 3′-flank structures',
         6: 'Conservation of co-oriented gene pairs', 7: 'Concordance between <i>T. obliquus</i> assemblies',
         8: 'Cross-genome 8-mer enrichment', 9: 'Logos of all 5′-leader footprints', 10: 'Workflow framework'}
for k in list(range(1, 11))+[11]:
    s.append(Paragraph(f'<b>Fig. S{k}</b> ' + SHORT[k], st['n']))
for k in list(range(1, 11))+[11]:
    s.append(PageBreak())
    im = PI.open(FIG / f'FigS{k}.jpg'); w, h = im.size; dpi = im.info.get('dpi', (600, 600))[0]
    pw, ph = w / dpi * 25.4 * mm, h / dpi * 25.4 * mm
    capH = 45 * mm if k == 10 else 28 * mm
    f = min(1, W / pw, (H - capH) / ph)
    s.append(KeepTogether([Image(str(FIG / f'FigS{k}.jpg'), width=pw * f, height=ph * f),
                           Paragraph(f'<b>Fig. S{k}</b> ' + CAP[k], st['c'])]))
doc.build(s, onFirstPage=header, onLaterPages=header)
print(OUT, OUT.stat().st_size)
