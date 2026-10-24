#!/usr/bin/env python3
"""Link the Zotero citations of the manuscript to the items in YOUR Zotero library.

Usage (Python 3.8+, standard library only; close Word first, Zotero may stay open):

    python3 relink_zotero_citations.py Manuscript_Afalcatus_..._Zotero.docx
    python3 relink_zotero_citations.py MANUSCRIPT.docx --zotero-db "~/Zotero/zotero.sqlite" --out linked.docx

Before running, import references_zotero.rdf (or references.ris) into Zotero
(File > Import...). The script reads a *copy* of zotero.sqlite (never modifies it),
finds each cited item in your personal library by DOI, then PMID, then title, and
rewrites the 'uris' of every ZOTERO_ITEM field so that they point to your items
(http://zotero.org/users/<userID>/items/<KEY> for synced accounts and
 http://zotero.org/users/local/<localUserKey>/items/<KEY>). Embedded item data are
kept as a fallback. Afterwards open the document in Word and click Zotero > Refresh:
the citations and the bibliography are then live links to your library.
"""
import argparse
import html
import json
import os
import re
import shutil
import sqlite3
import sys
import tempfile
import unicodedata
import zipfile


def norm_title(t):
    t = unicodedata.normalize('NFKD', re.sub(r'<[^>]+>', '', t or '')).encode('ascii', 'ignore').decode().lower()
    return re.sub(r'[^a-z0-9]', '', t)


def default_db():
    home = os.path.expanduser('~')
    for p in (os.path.join(home, 'Zotero', 'zotero.sqlite'),
              os.path.join(home, 'Documents', 'Zotero', 'zotero.sqlite')):
        if os.path.exists(p):
            return p
    return os.path.join(home, 'Zotero', 'zotero.sqlite')


def library_index(db):
    tmp = tempfile.mkdtemp()
    copy = os.path.join(tmp, 'zotero.sqlite')
    shutil.copy2(db, copy)                       # Zotero locks its DB; work on a copy
    con = sqlite3.connect(copy)
    cur = con.cursor()
    settings = dict(((k, v) for k, v in cur.execute("SELECT key, value FROM settings WHERE setting='account'")))
    user_lib = cur.execute("SELECT libraryID FROM libraries WHERE type='user'").fetchone()[0]
    rows = cur.execute("""
        SELECT i.key, f.fieldName, v.value
        FROM items i
        JOIN itemData d ON d.itemID = i.itemID
        JOIN fields f ON f.fieldID = d.fieldID
        JOIN itemDataValues v ON v.valueID = d.valueID
        WHERE i.libraryID = ? AND i.itemID NOT IN (SELECT itemID FROM deletedItems)
          AND f.fieldName IN ('DOI', 'title', 'extra')""", (user_lib,)).fetchall()
    con.close(); shutil.rmtree(tmp, ignore_errors=True)
    by_doi, by_pmid, by_title = {}, {}, {}
    for key, field, value in rows:
        if field == 'DOI':
            by_doi[value.lower().strip()] = key
        elif field == 'title':
            by_title[norm_title(value)] = key
        elif field == 'extra':
            for m in re.finditer(r'PMID:\s*(\d+)', value):
                by_pmid[m.group(1)] = key
            for m in re.finditer(r'DOI:\s*(\S+)', value):
                by_doi.setdefault(m.group(1).lower(), key)
    bases = []
    if settings.get('userID'):
        bases.append(f"http://zotero.org/users/{settings['userID']}/items/")
    if settings.get('localUserKey'):
        bases.append(f"http://zotero.org/users/local/{settings['localUserKey']}/items/")
    if not bases:
        sys.exit('Could not read the library identifiers from zotero.sqlite')
    return bases, by_doi, by_pmid, by_title


def match(item, by_doi, by_pmid, by_title):
    d = item.get('itemData', {})
    if d.get('DOI') and d['DOI'].lower() in by_doi:
        return by_doi[d['DOI'].lower()]
    if d.get('PMID') and str(d['PMID']) in by_pmid:
        return by_pmid[str(d['PMID'])]
    return by_title.get(norm_title(d.get('title', '')))


TOKEN = re.compile(r'<w:fldChar\b[^>]*w:fldCharType="(begin|separate|end)"[^>]*/>|(<w:instrText\b[^>]*>)(.*?)(</w:instrText>)', re.S)


def relink_xml(xml, index):
    bases, by_doi, by_pmid, by_title = index
    stack, edits = [], []
    stats = {'fields': 0, 'linked': 0, 'unmatched': set()}

    def process(field):
        pieces = field['pieces']
        if not pieces:
            return
        instr = ''.join(html.unescape(p[3]) for p in pieces)
        if 'ADDIN ZOTERO_ITEM CSL_CITATION' not in instr:
            return
        head, js = instr.split('CSL_CITATION', 1)
        data = json.loads(js.strip())
        stats['fields'] += 1
        for ci in data.get('citationItems', []):
            key = match(ci, by_doi, by_pmid, by_title)
            if key:
                ci['uris'] = [b + key for b in bases]
                stats['linked'] += 1
            else:
                stats['unmatched'].add(ci.get('itemData', {}).get('title', '?')[:70])
        new = f"{head}CSL_CITATION {json.dumps(data, ensure_ascii=False)} "
        s0, e0, o0, _, c0 = pieces[0]
        edits.append((s0, e0, o0 + html.escape(new, quote=False) + c0))
        for s1, e1, o1, _, c1 in pieces[1:]:
            edits.append((s1, e1, o1 + c1))

    for m in TOKEN.finditer(xml):
        kind = m.group(1)
        if kind == 'begin':
            stack.append({'pieces': [], 'closed': False})
        elif kind == 'separate':
            if stack and not stack[-1]['closed']:
                process(stack[-1]); stack[-1]['closed'] = True
        elif kind == 'end':
            if stack:
                f = stack.pop()
                if not f['closed']:
                    process(f)
        elif stack and not stack[-1]['closed']:
            stack[-1]['pieces'].append((m.start(), m.end(), m.group(2), m.group(3), m.group(4)))
    for s, e, rep in sorted(edits, reverse=True):
        xml = xml[:s] + rep + xml[e:]
    return xml, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('docx')
    ap.add_argument('--zotero-db', default=default_db())
    ap.add_argument('--out')
    a = ap.parse_args()
    db = os.path.expanduser(a.zotero_db)
    if not os.path.exists(db):
        sys.exit(f'zotero.sqlite not found at {db}; pass --zotero-db (Zotero > Settings > Advanced > Files and Folders shows the data directory).')
    out = a.out or re.sub(r'\.docx$', '', a.docx) + '_linked.docx'
    index = library_index(db)
    zin = zipfile.ZipFile(a.docx)
    total = {'fields': 0, 'linked': 0, 'unmatched': set()}
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            data = zin.read(info.filename)
            if re.match(r'word/(document|footnotes|endnotes)\.xml$', info.filename):
                xml, st = relink_xml(data.decode('utf-8'), index)
                data = xml.encode('utf-8')
                total['fields'] += st['fields']; total['linked'] += st['linked']; total['unmatched'] |= st['unmatched']
            zout.writestr(info, data)
    print(f"Citation fields: {total['fields']}; cited items linked to your library: {total['linked']}")
    if total['unmatched']:
        print('Not found in your library (import references_zotero.rdf first):')
        for t in sorted(total['unmatched']):
            print('  -', t)
    print(f'Written: {out}\nOpen it in Word and click Zotero > Refresh.')


if __name__ == '__main__':
    main()
