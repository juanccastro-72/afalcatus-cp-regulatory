#!/usr/bin/env python3
"""Build the manuscript .docx with live Zotero citation and bibliography fields.

Every author-year citation group in manuscript.md is converted into a Word field
  ADDIN ZOTERO_ITEM CSL_CITATION {...}
carrying full CSL-JSON item data (from references/references.json, built from
Europe PMC / Crossref records), and the reference list becomes an
  ADDIN ZOTERO_BIBL ... CSL_BIBLIOGRAPHY
field. Zotero document preferences (style: Springer Basic author-date, as used by
Plant Molecular Biology and Planta; field mode) are stored as ZOTERO_PREF_n custom document properties, so the
Zotero Word plugin recognises the document and can Refresh it. Item URIs use the
placeholder library 'users/local/PLANTAAF' until relink_zotero_citations.py
rewrites them to the items imported into your own Zotero library.
"""
import json
import random
import re
import subprocess
import sys
import unicodedata
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

HERE = Path(__file__).resolve().parent.parent
REFS = json.loads((HERE / 'references/references.json').read_text())
KEYCHARS = '23456789ABCDEFGHIJKLMNPQRSTUVWXYZ'
rng = random.Random(20261024)
for it in REFS:
    it['zkey'] = ''.join(rng.choice(KEYCHARS) for _ in range(8))
LIB = 'http://zotero.org/users/local/PLANTAAF'
STYLE = 'http://www.zotero.org/styles/springer-basic-author-date'


def norm(s):
    return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode().lower().replace(' ', '').replace('-', '')


INDEX = {}
for it in REFS:
    INDEX[(norm(it['author'][0].get('family', '')), str(it['issued']['date-parts'][0][0]))] = it

CIT = re.compile(r"^(?P<prefix>.*?)(?P<name>(?:de )?[A-Z][\wÀ-ſ\-]+)(?P<etal> et al\.| and [A-Z][\wÀ-ſ\-]+)? (?P<years>\d{4}(?:, \d{4})*)$")


def parse_group(g):
    """'Klein et al. 1992; TSS from Anthonisen et al. 2001' -> list of [prefix, item, suffix] or None.

    Non-citation parts before a citation become its prefix; parts after the last citation its suffix."""
    out, pending = [], []
    for part in [p.strip() for p in g.split(';')]:
        m = CIT.match(part)
        its = None
        if m:
            its = [INDEX.get((norm(m.group('name')), y)) for y in m.group('years').split(', ')]
            if None in its:
                its = None
        if its is None:
            pending.append(part)
            continue
        pre = ('; '.join(pending) + '; ' if pending else '') + m.group('prefix')
        pre = pre.rstrip()
        pending = []
        for k, it in enumerate(its):
            out.append([pre if k == 0 else '', it, ''])
    if not out:
        return None
    if pending:
        out[-1][2] = '; ' + '; '.join(pending)
    return out


def short_author(it):
    a = it['author']
    fam = lambda x: x.get('family', x.get('literal'))
    if len(a) == 1:
        return fam(a[0])
    if len(a) == 2:
        return f'{fam(a[0])} and {fam(a[1])}'
    return f'{fam(a[0])} et al.'


def formatted(cits):
    parts, last = [], None
    for pre, it, suf in cits:
        y = str(it['issued']['date-parts'][0][0])
        if last == short_author(it) and parts and not pre:
            parts[-1] += f', {y}{suf}'
        else:
            parts.append(f"{(pre + ' ') if pre else ''}{short_author(it)} {y}{suf}")
        last = short_author(it)
    return re.sub(r'\*([^*]+)\*', r'\1', '(' + '; '.join(parts) + ')')


def item_data(it):
    d = {k: v for k, v in it.items() if k not in ('zkey',)}
    d['id'] = it['zkey']
    return d


def citation_json(cits, cid):
    items = []
    for pre, it, suf in cits:
        ci = {'id': it['zkey'], 'uris': [f"{LIB}/items/{it['zkey']}"], 'itemData': item_data(it)}
        ital = lambda t: re.sub(r'\*([^*]+)\*', r'<i>\1</i>', t)
        if pre:
            ci['prefix'] = ital(pre)
        if suf:
            ci['suffix'] = ital(suf)
        items.append(ci)
    f = formatted(cits)
    return json.dumps({'citationID': cid, 'properties': {'formattedCitation': f, 'plainCitation': f, 'noteIndex': 0},
                       'citationItems': items, 'schema': 'https://github.com/citation-style-language/schema/raw/master/csl-citation.json'},
                      ensure_ascii=False), f


def bib_entry(it):
    fam = lambda a: a.get('family', a.get('literal', ''))
    ini = lambda a: ''.join(p[0] for p in re.split(r'[\s\.\-]+', a.get('given', '')) if p)
    au = ', '.join(f'{fam(a)} {ini(a)}'.strip() for a in it['author'])
    s = f"{au} ({it['issued']['date-parts'][0][0]}) {it['title']}. {it.get('container-title-short') or it['container-title']}"
    if it.get('volume'):
        s += f" {it['volume']}"
    if it.get('page'):
        s += f":{it['page'].replace('-', '–')}"
    if it.get('DOI'):
        s += f". https://doi.org/{it['DOI']}"
    return s


def main(md_in, docx_out, ref_docx):
    text = Path(md_in).read_text()
    body, refs_and_after = text.split('## References', 1)
    after = refs_and_after[refs_and_after.index('\n## '):]
    fields, cited = {}, []

    def repl(m):
        g = m.group(1)
        cits = parse_group(g)
        if not cits:
            return m.group(0)
        tok = f'ZCITE{len(fields):04d}'
        fields[tok] = cits
        cited.extend(it['id'] for _, it, _s in cits)
        return tok

    body2 = re.sub(r'\(([^()]*\d{4}[^()]*)\)', repl, body)
    after2 = re.sub(r'\(([^()]*\d{4}[^()]*)\)', repl, after)
    uncited = [it['id'] for it in REFS if it['id'] not in cited]
    if uncited:
        print('WARNING uncited references:', uncited)
    for m in re.finditer(r'[A-Z][\w\-]+ (?:et al\.|and [A-Z]\w+) \d{4}', body2 + after2):
        print('WARNING unconverted citation:', (body2 + after2)[max(0, m.start() - 60):m.end() + 10].replace('\n', ' '))
    refs_sorted = sorted([it for it in REFS if it['id'] in cited], key=lambda i: (norm(re.sub(r'^de ', '', i['author'][0].get('family', ''))), i['issued']['date-parts'][0][0]))
    md = body2 + '## References\n\nZBIBLSTART\n\n' + '\n\n'.join(f'ZBIBITEM{i:03d}' for i in range(len(refs_sorted))) + '\n\n' + after2
    tmp_md = Path(docx_out).with_suffix('.tmp.md'); tmp_md.write_text(md)
    raw = Path(docx_out).with_suffix('.raw.docx')
    subprocess.run(['pandoc', str(tmp_md), '-o', str(raw), f'--reference-doc={ref_docx}'], check=True)
    tmp_md.unlink()
    d = Document(raw)

    def mk_run(rpr, kind=None, text=None):
        r = OxmlElement('w:r')
        if rpr is not None:
            r.append(__import__('copy').deepcopy(rpr))
        if kind in ('begin', 'separate', 'end'):
            e = OxmlElement('w:fldChar'); e.set(qn('w:fldCharType'), kind); r.append(e)
        elif kind == 'instr':
            e = OxmlElement('w:instrText'); e.set(qn('xml:space'), 'preserve'); e.text = text; r.append(e)
        else:
            e = OxmlElement('w:t'); e.set(qn('xml:space'), 'preserve'); e.text = text; r.append(e)
        return r

    def all_paragraphs():
        for p in d.paragraphs:
            yield p
        for t in d.tables:
            for row in t.rows:
                for c in row.cells:
                    yield from c.paragraphs

    n_done = 0
    for p in list(all_paragraphs()):
        for r in list(p._p.iter(qn('w:r'))):
            t = r.find(qn('w:t'))
            if t is None or 'ZCITE' not in (t.text or ''):
                continue
            rpr = r.find(qn('w:rPr'))
            pieces = re.split(r'(ZCITE\d{4})', t.text)
            new = []
            for pc in pieces:
                if not pc:
                    continue
                if pc in fields:
                    cid = ''.join(rng.choice('abcdefghijklmnopqrstuvwxyz0123456789') for _ in range(8))
                    js, f = citation_json(fields[pc], cid)
                    new += [mk_run(rpr, 'begin'), mk_run(rpr, 'instr', f' ADDIN ZOTERO_ITEM CSL_CITATION {js} '),
                            mk_run(rpr, 'separate'), mk_run(rpr, None, f), mk_run(rpr, 'end')]
                    n_done += 1
                else:
                    new.append(mk_run(rpr, None, pc))
            parent = r.getparent(); idx = parent.index(r)
            parent.remove(r)
            for k, nr in enumerate(new):
                parent.insert(idx + k, nr)
    # bibliography field spanning the reference paragraphs
    paras = d.paragraphs
    start = next(i for i, p in enumerate(paras) if p.text.strip() == 'ZBIBLSTART')
    bib_paras = paras[start + 1:start + 1 + len(refs_sorted)]
    for p, it in zip(bib_paras, refs_sorted):
        for r in list(p._p.iter(qn('w:r'))):
            r.getparent().remove(r)
        p._p.append(mk_run(None, None, bib_entry(it)))
        p.paragraph_format.left_indent = Pt(18); p.paragraph_format.first_line_indent = Pt(-18)
    first, last = bib_paras[0]._p, bib_paras[-1]._p
    bibjson = json.dumps({'uncited': [], 'omitted': [], 'custom': []})
    for k, nr in enumerate([mk_run(None, 'begin'), mk_run(None, 'instr', f' ADDIN ZOTERO_BIBL {bibjson} CSL_BIBLIOGRAPHY '), mk_run(None, 'separate')]):
        ppr = first.find(qn('w:pPr'))
        first.insert((first.index(ppr) + 1 if ppr is not None else 0) + k, nr)
    last.append(mk_run(None, 'end'))
    paras[start]._p.getparent().remove(paras[start]._p)
    # line numbers, page numbers, table borders (Springer layout)
    sect = d.sections[0]
    ln = OxmlElement('w:lnNumType'); ln.set(qn('w:countBy'), '1'); ln.set(qn('w:restart'), 'continuous'); ln.set(qn('w:distance'), '284'); sect._sectPr.append(ln)
    fp = sect.footer.paragraphs[0]; fp.alignment = 1
    for kind, txt in (('begin', None), ('instr', 'PAGE'), ('end', None)):
        fp._p.append(mk_run(None, kind, txt))
    for t in d.tables:
        b = OxmlElement('w:tblBorders')
        for e in ('top', 'bottom', 'insideH'):
            el = OxmlElement(f'w:{e}'); el.set(qn('w:val'), 'single'); el.set(qn('w:sz'), '4'); el.set(qn('w:color'), '808080'); b.append(el)
        t._tbl.tblPr.append(b)
        for row in t.rows:
            for c in row.cells:
                for para in c.paragraphs:
                    para.paragraph_format.line_spacing = 1.0
                    for run in para.runs:
                        run.font.size = Pt(8.5)
    d.save(raw)
    # Zotero document preferences as custom properties
    pref = (f'<data data-version="3" zotero-version="7.0"><session id="{"".join(rng.choice("ABCDEFGHJKLMNPQRSTUVWXYZ23456789") for _ in range(8))}"/>'
            f'<style id="{STYLE}" hasBibliography="1" bibliographyStyleHasBeenSet="1"/>'
            '<prefs><pref name="fieldType" value="Field"/><pref name="automaticJournalAbbreviations" value="true"/></prefs></data>')
    chunks = [pref[i:i + 255] for i in range(0, len(pref), 255)]
    props = ''.join(f'<property fmtid="{{D5CDD505-2E9C-101B-9397-08002B2CF9AE}}" pid="{i + 2}" name="ZOTERO_PREF_{i + 1}">'
                    f'<vt:lpwstr>{escape(c)}</vt:lpwstr></property>' for i, c in enumerate(chunks))
    custom = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/custom-properties" '
              f'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">{props}</Properties>')
    zin = zipfile.ZipFile(raw)
    with zipfile.ZipFile(docx_out, 'w', zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            data = zin.read(info.filename)
            if info.filename == '[Content_Types].xml' and b'custom-properties' not in data:
                data = data.replace(b'</Types>', b'<Override PartName="/docProps/custom.xml" ContentType="application/vnd.openxmlformats-officedocument.custom-properties+xml"/></Types>')
            if info.filename == '_rels/.rels' and b'custom-properties' not in data:
                data = data.replace(b'</Relationships>', b'<Relationship Id="rIdZoteroCustom" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/custom-properties" Target="docProps/custom.xml"/></Relationships>')
            if info.filename == 'docProps/custom.xml':
                continue
            zout.writestr(info, data)
        zout.writestr('docProps/custom.xml', custom)
    raw.unlink()
    # item key table for the relink script / user
    (HERE / 'references/zotero_placeholder_keys.json').write_text(json.dumps({it['zkey']: it['id'] for it in REFS}, indent=1))
    print(f'{n_done} citation fields; {len(refs_sorted)} bibliography entries')


if __name__ == '__main__':
    main(*sys.argv[1:4])
