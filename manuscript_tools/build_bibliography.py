#!/usr/bin/env python3
"""Build a verified bibliography (CSL-JSON, RIS, BibTeX, Zotero RDF) from database records.

Every PubMed-indexed reference is fetched from Europe PMC by PMID (authors, title,
journal, volume, issue, pages, year, DOI, PMID, PMCID). Two non-PubMed items are
taken from their Crossref records (verified in this session). Nothing is typed by hand
except those two Crossref transcriptions.
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

PMIDS = """3870875 11453063 7584402 19458158 33760053 16231149 2152108 8388079 2897716 16638149 20003500 19505945
28985404 15840812 33521263 19304878 9520261 26903036 17996019 39688456 17324271 10567573 41184475 28481363 23329690
1565638 16410618 14704338 8756597 24335082 22115189 15123812 7798310 29722901 39871102 32011700 8039511 10330479
20230484 21535358 9484442 33989339 23656324 8455609 41484447 10758498 20525638 21540409""".split()

CROSSREF_ITEMS = [  # transcribed from api.crossref.org/works/<DOI> (retrieved 24 Oct 2026)
    dict(type='article-journal', title='Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing',
         **{'container-title': 'Journal of the Royal Statistical Society Series B: Statistical Methodology',
            'container-title-short': 'J R Stat Soc Ser B'},
         author=[dict(family='Benjamini', given='Yoav'), dict(family='Hochberg', given='Yosef')],
         issued={'date-parts': [[1995]]}, volume='57', issue='1', page='289-300', DOI='10.1111/j.2517-6161.1995.tb02031.x'),
    dict(type='article-journal', title='Isolation and Characterization of Native Microalgae from the Peruvian Amazon with Potential for Biodiesel Production',
         **{'container-title': 'Energies', 'container-title-short': 'Energies'},
         author=[dict(family='Cobos', given='Marianela'), dict(family='Paredes', given='Jae D.'), dict(family='Maddox', given='J. Dylan'),
                 dict(family='Vargas-Arana', given='Gabriel'), dict(family='Flores', given='Leenin'), dict(family='Aguilar', given='Carla'),
                 dict(family='Marapara', given='Jorge L.'), dict(family='Castro', given='Juan C.')],
         issued={'date-parts': [[2017]]}, volume='10', issue='2', page='224', DOI='10.3390/en10020224'),
]


def epmc(pmid):
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode(
        dict(query=f'EXT_ID:{pmid} AND SRC:MED', format='json', resultType='core'))
    for k in range(6):
        try:
            return json.load(urllib.request.urlopen(url, timeout=60))['resultList']['result'][0]
        except Exception:
            time.sleep(5 * (k + 1))
    raise RuntimeError(pmid)


def clean_title(t):
    t = re.sub(r'<[^>]+>', '', t).strip().rstrip('.')
    return t


def to_csl(r):
    ji = r.get('journalInfo', {})
    j = ji.get('journal', {})
    authors = []
    for a in r.get('authorList', {}).get('author', []):
        if 'lastName' in a:
            authors.append(dict(family=a['lastName'], given=a.get('firstName') or a.get('initials', '')))
        elif 'collectiveName' in a:
            authors.append(dict(literal=a['collectiveName']))
    item = dict(type='article-journal', title=clean_title(r['title']), author=authors,
                **{'container-title': j.get('title', '').rstrip('.'), 'container-title-short': j.get('isoabbreviation', '').replace('.', '')},
                issued={'date-parts': [[int(r['pubYear'])]]}, volume=ji.get('volume', ''), issue=ji.get('issue', ''),
                page=r.get('pageInfo', ''), DOI=r.get('doi', ''), PMID=r['pmid'])
    if r.get('pmcid'):
        item['PMCID'] = r['pmcid']
    if r.get('journalInfo', {}).get('journal', {}).get('issn'):
        item['ISSN'] = j['issn']
    return {k: v for k, v in item.items() if v not in ('', None, [])}


def citekey(item, used):
    fam = item['author'][0].get('family', item['author'][0].get('literal', 'Anon'))
    import unicodedata
    fam = unicodedata.normalize('NFKD', fam).encode('ascii', 'ignore').decode()
    base = re.sub(r'[^A-Za-z]', '', fam) + str(item['issued']['date-parts'][0][0])
    k, n = base, 0
    while k in used:
        n += 1; k = base + 'abcdefgh'[n - 1]
    used.add(k)
    return k


def ris(items):
    out = []
    for it in items:
        L = ['TY  - JOUR']
        for a in it['author']:
            L.append(f"AU  - {a['family']}, {a['given']}" if 'family' in a else f"AU  - {a['literal']}")
        L += [f"TI  - {it['title']}", f"T2  - {it['container-title']}", f"J2  - {it.get('container-title-short', '')}",
              f"PY  - {it['issued']['date-parts'][0][0]}", f"VL  - {it.get('volume', '')}", f"IS  - {it.get('issue', '')}"]
        if it.get('page'):
            p = it['page'].split('-')
            L.append(f'SP  - {p[0]}')
            if len(p) > 1:
                L.append(f'EP  - {p[1]}')
        if it.get('DOI'):
            L += [f"DO  - {it['DOI']}", f"UR  - https://doi.org/{it['DOI']}"]
        if it.get('PMID'):
            L.append(f"AN  - PMID:{it['PMID']}")
        if it.get('ISSN'):
            L.append(f"SN  - {it['ISSN']}")
        L.append(f"N1  - Citation key: {it['id']}")
        L.append('ER  - ')
        out.append('\n'.join(x for x in L if not x.endswith('-  ') and not re.match(r'^\w\w  - $', x)))
    return '\n\n'.join(out) + '\n'


def bibtex(items):
    out = []
    esc = lambda s: s.replace('&', r'\&')
    for it in items:
        au = ' and '.join(f"{a['family']}, {a['given']}" if 'family' in a else '{' + a['literal'] + '}' for a in it['author'])
        f = [f"  author = {{{au}}}", f"  title = {{{{{esc(it['title'])}}}}}", f"  journal = {{{esc(it['container-title'])}}}",
             f"  year = {{{it['issued']['date-parts'][0][0]}}}"]
        for k in ('volume', 'issue', 'page'):
            if it.get(k):
                f.append(f"  {'number' if k == 'issue' else 'pages' if k == 'page' else k} = {{{it[k].replace('-', '--')}}}")
        if it.get('DOI'):
            f.append(f"  doi = {{{it['DOI']}}}")
        if it.get('PMID'):
            f.append(f"  pmid = {{{it['PMID']}}}")
        out.append(f"@article{{{it['id']},\n" + ',\n'.join(f) + '\n}')
    return '\n\n'.join(out) + '\n'


def rdf(items):
    from xml.sax.saxutils import escape as x
    head = ('<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:z="http://www.zotero.org/namespaces/export#" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:bib="http://purl.org/net/biblio#" xmlns:foaf="http://xmlns.com/foaf/0.1/" '
            'xmlns:prism="http://prismstandard.org/namespaces/1.2/basic/" xmlns:link="http://purl.org/rss/1.0/modules/link/">\n')
    body = []
    for it in items:
        about = f"https://doi.org/{it['DOI']}" if it.get('DOI') else f"#item_{it['id']}"
        authors = ''.join(
            f"<rdf:li><foaf:Person><foaf:surname>{x(a['family'])}</foaf:surname><foaf:givenName>{x(a['given'])}</foaf:givenName></foaf:Person></rdf:li>"
            if 'family' in a else f"<rdf:li><foaf:Organization><foaf:name>{x(a['literal'])}</foaf:name></foaf:Organization></rdf:li>"
            for a in it['author'])
        jr = f"<dcterms:isPartOf><bib:Journal><dc:title>{x(it['container-title'])}</dc:title>" \
             f"{('<dcterms:alternative>' + x(it['container-title-short']) + '</dcterms:alternative>') if it.get('container-title-short') else ''}" \
             f"{('<prism:volume>' + x(it['volume']) + '</prism:volume>') if it.get('volume') else ''}" \
             f"{('<prism:number>' + x(it['issue']) + '</prism:number>') if it.get('issue') else ''}" \
             f"{('<dc:identifier>ISSN ' + x(it['ISSN']) + '</dc:identifier>') if it.get('ISSN') else ''}</bib:Journal></dcterms:isPartOf>"
        ids = (f"<dc:identifier>DOI {x(it['DOI'])}</dc:identifier>" if it.get('DOI') else '')
        extra = f"Citation Key: {it['id']}" + (f"\nPMID: {it['PMID']}" if it.get('PMID') else '') + (f"\nPMCID: {it['PMCID']}" if it.get('PMCID') else '')
        body.append(f'<bib:Article rdf:about="{x(about)}"><z:itemType>journalArticle</z:itemType>{jr}'
                    f'<bib:authors><rdf:Seq>{authors}</rdf:Seq></bib:authors><dc:title>{x(it["title"])}</dc:title>'
                    f'{("<bib:pages>" + x(it["page"]) + "</bib:pages>") if it.get("page") else ""}'
                    f'<dc:date>{it["issued"]["date-parts"][0][0]}</dc:date>{ids}'
                    f'<dc:description>{x(extra)}</dc:description></bib:Article>')
    return head + '\n'.join(body) + '\n</rdf:RDF>\n'


def main():
    out = Path(__file__).resolve().parent.parent / 'manuscript' / 'references'
    out.mkdir(parents=True, exist_ok=True)
    items = [to_csl(epmc(p)) for p in PMIDS] + CROSSREF_ITEMS
    for it in items:
        it['title'] = it['title'].replace('mRNas', 'mRNAs')
        if it.get('container-title', '').startswith('Science (New York'):
            it['container-title'] = 'Science'
    for it in items:  # corrections from Crossref where Europe PMC author lists are incomplete
        if it.get('DOI') == '10.1126/science.2897716' and not any(a.get('family') == 'Sanford' for a in it['author']):
            it['author'].append(dict(family='Sanford', given='J. C.'))
        if it.get('page') == 'W202-8':
            it['page'] = 'W202-W208'
    used = set()
    for it in sorted(items, key=lambda i: (i['author'][0].get('family', ''), i['issued']['date-parts'][0][0])):
        it['id'] = citekey(it, used)
    items.sort(key=lambda i: i['id'])
    (out / 'references.json').write_text(json.dumps(items, ensure_ascii=False, indent=1))
    (out / 'references.ris').write_text(ris(items))
    (out / 'references.bib').write_text(bibtex(items))
    (out / 'references_zotero.rdf').write_text(rdf(items))
    print(len(items), 'items')


if __name__ == '__main__':
    main()
