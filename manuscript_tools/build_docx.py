"""Build the Planta-formatted manuscript .docx from manuscript.md (pandoc + python-docx post-processing)."""
import subprocess, sys
from pathlib import Path
from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt
M = Path(sys.argv[1]); out = Path(sys.argv[2]); ref = M.parent/'build/ref.docx'
raw = M.parent/'build/manuscript_raw.docx'
subprocess.run(['pandoc', str(M), '-o', str(raw), f'--reference-doc={ref}'], check=True)
d = Document(raw); sect = d.sections[0]
ln = OxmlElement('w:lnNumType'); ln.set(qn('w:countBy'), '1'); ln.set(qn('w:restart'), 'continuous'); ln.set(qn('w:distance'), '284'); sect._sectPr.append(ln)
p = sect.footer.paragraphs[0]; p.alignment = 1
for typ, txt in (('begin', None), (None, 'PAGE'), ('end', None)):
    r = p.add_run()
    if typ: e = OxmlElement('w:fldChar'); e.set(qn('w:fldCharType'), typ)
    else: e = OxmlElement('w:instrText'); e.set(qn('xml:space'), 'preserve'); e.text = txt
    r._r.append(e)
for t in d.tables:
    b = OxmlElement('w:tblBorders')
    for e in ('top', 'bottom', 'insideH'):
        el = OxmlElement(f'w:{e}'); el.set(qn('w:val'), 'single'); el.set(qn('w:sz'), '4'); el.set(qn('w:color'), '808080'); b.append(el)
    t._tbl.tblPr.append(b)
    for row in t.rows:
        for c in row.cells:
            for para in c.paragraphs:
                para.paragraph_format.line_spacing = 1.0
                for run in para.runs: run.font.size = Pt(8.5)
d.save(out)
