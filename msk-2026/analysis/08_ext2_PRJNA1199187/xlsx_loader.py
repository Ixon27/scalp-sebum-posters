import zipfile,re,xml.etree.ElementTree as ET
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
def col2i(ref):
    n=0
    for ch in ref: n=n*26+(ord(ch)-64)
    return n-1
def load(path):
    z=zipfile.ZipFile(path)
    sst=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        r=ET.fromstring(z.read('xl/sharedStrings.xml'))
        for si in r.findall(NS+'si'): sst.append(''.join(t.text or '' for t in si.iter(NS+'t')))
    wb=z.read('xl/workbook.xml').decode()
    names=re.findall(r'<sheet name="([^"]+)"[^>]*r:id="(rId\d+)"',wb)
    rels=dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', z.read('xl/_rels/workbook.xml.rels').decode()))
    out={}
    for nm,rid in names:
        tgt=rels[rid]
        if not tgt.startswith('xl/'): tgt='xl/'+tgt.lstrip('/')
        rt=ET.fromstring(z.read(tgt)); grid=[]
        for row in rt.iter(NS+'row'):
            d={}
            for c in row.findall(NS+'c'):
                v=c.find(NS+'v')
                if v is None: continue
                d[col2i(re.match(r'[A-Z]+',c.get('r')).group(0))]= sst[int(v.text)] if c.get('t')=='s' else v.text
            if d: grid.append([d.get(i,'') for i in range(max(d)+1)])
        out[nm]=grid
    return out
