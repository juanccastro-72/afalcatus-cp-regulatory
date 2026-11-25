"""Shared publication figure style: Arial-metric font, validated categorical palette,
600-dpi JPEG + fully vector SVG (text converted to paths)."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.dpi': 600, 'figure.dpi': 150,
    'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'legend.frameon': False,
})

# Plant Molecular Biology (Springer) artwork settings: CPREG_FIG_STYLE=PMB
#   * final width <= 174 mm (double column) or 84 mm (single), height <= 234 mm
#   * Arial-metric sans-serif lettering, bold lowercase panel labels
#   * files named Fig1, Fig2 ... FigS1 ...: 600-dpi JPEG, vector PDF (fonts embedded) and vector SVG
import os as _os
PMB = _os.environ.get('CPREG_FIG_STYLE', '').upper() == 'PMB'
FIG_DIR = _os.environ.get('CPREG_FIG_DIR')
if PMB:
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8, 'xtick.labelsize': 7,
                         'ytick.labelsize': 7, 'legend.fontsize': 7})
MAX_W_MM, MAX_H_MM = 174, 234

# Categorical slots (validated reference palette, fixed order)
SPH, CHL, OUT = '#2a78d6', '#eb6834', '#1baf7a'
GROUP_COLOR = {'Sphaeropleales': SPH, 'Chlamydomonadales': CHL, 'Outgroups (OCC/Trebouxiophyceae)': OUT}
FOCAL = '#111111'
INK, MUTED, GRID = '#222222', '#6b6a66', '#d9d8d4'
DIV = 'RdBu_r'          # diverging (z-scores): blue = negative (more structured), red = positive
SEQ = 'Blues'

MM = 1 / 25.4
SINGLE, DOUBLE = 84 * MM, 174 * MM


def _save_pmb(fig, outdir, name, manifest, caption):
    """Save Fig<n> in PMB format, scaled (never enlarged) to fit 174 x 234 mm."""
    import io
    import re as _re
    from pypdf import PdfReader, PdfWriter, Transformation
    from PIL import Image
    short = name.split('_')[0]
    sub = 'supplementary' if short.startswith('FigS') else 'main'
    out = Path(FIG_DIR) / sub if FIG_DIR else Path(outdir); out.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO(); fig.savefig(buf, format='pdf', bbox_inches='tight'); buf.seek(0)
    page = PdfReader(buf).pages[0]
    w_pt, h_pt = float(page.mediabox.width), float(page.mediabox.height)
    s = min(1.0, MAX_W_MM / (w_pt / 72 * 25.4), MAX_H_MM / (h_pt / 72 * 25.4))
    # vector PDF
    if s < 1:
        page.add_transformation(Transformation().scale(s, s))
        page.mediabox.upper_right = (w_pt * s, h_pt * s); page.mediabox.lower_left = (0, 0)
        page.cropbox = page.mediabox
    wr = PdfWriter(); wr.add_page(page); wr.add_metadata({'/Title': short})
    with open(out / f'{short}.pdf', 'wb') as fh:
        wr.write(fh)
    # vector SVG: same drawing, physical size scaled through width/height attributes
    sb = io.StringIO(); fig.savefig(sb, format='svg', bbox_inches='tight'); svg = sb.getvalue()
    def _sc(m):
        return f'{m.group(1)}="{float(m.group(2)) * s * 25.4 / 72:.2f}mm"'
    head, rest = svg.split('<svg', 1)
    tag, body = rest.split('>', 1)
    tag = _re.sub(r'(width|height)="([0-9.]+)pt"', _sc, tag)
    (out / f'{short}.svg').write_text(head + '<svg' + tag + '>' + body)
    # 600-dpi JPEG at the final physical size
    jb = io.BytesIO(); fig.savefig(jb, format='png', dpi=600 * s, bbox_inches='tight'); jb.seek(0)
    Image.open(jb).convert('RGB').save(out / f'{short}.jpg', quality=95, dpi=(600, 600))
    plt.close(fig)
    if manifest is not None:
        manifest.append(dict(figure=short, source_name=name, jpeg=f'{short}.jpg', pdf=f'{short}.pdf', svg=f'{short}.svg',
                             width_mm=round(w_pt * s / 72 * 25.4, 1), height_mm=round(h_pt * s / 72 * 25.4, 1),
                             scale=round(s, 3), caption=caption))


def save(fig, outdir, name, manifest=None, caption=''):
    if PMB:
        return _save_pmb(fig, outdir, name, manifest, caption)
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f'{name}.svg', format='svg', bbox_inches='tight')
    fig.savefig(out / f'{name}.jpg', format='jpg', dpi=600, bbox_inches='tight', pil_kwargs={'quality': 95})
    plt.close(fig)
    if manifest is not None:
        manifest.append(dict(figure=name, jpeg=f'{name}.jpg', svg=f'{name}.svg', caption=caption))


def panel(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10 if PMB else 9, fontweight='bold', va='bottom', ha='left')
