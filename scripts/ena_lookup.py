import urllib.request, urllib.parse, pandas as pd, io, sys, json
def q(query, fields="accession,description,base_count,organelle,tax_id,scientific_name,topology,first_public"):
    url="https://www.ebi.ac.uk/ena/portal/api/search?"+urllib.parse.urlencode(dict(result="sequence",query=query,fields=fields,limit=200,format="tsv"))
    with urllib.request.urlopen(url, timeout=60) as r: t=r.read().decode()
    return pd.read_csv(io.StringIO(t),sep="\t") if t.strip() else pd.DataFrame()
names = sys.argv[1:]
out=[]
for n in names:
    df=q(f'tax_name("{n}") AND organelle="plastid:chloroplast"') if not n.isdigit() else q(f'tax_eq({n}) AND organelle="plastid:chloroplast"')
    if len(df): df=df[df.base_count>50000]
    print("##",n); print(df.to_string() if len(df) else "none")
