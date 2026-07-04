import csv, json, subprocess, time, urllib.parse
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
lfc=rdMolStandardize.LargestFragmentChooser(); unch=rdMolStandardize.Uncharger()
def key(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None: return None,None
    m=rdMolStandardize.Cleanup(m); m=lfc.choose(m); m=unch.uncharge(m)
    try: ik=inchi.InchiToInchiKey(inchi.MolToInchi(m))
    except Exception: return None,None
    return (ik, ik.split('-')[0]) if ik else (None,None)
def pubchem(name):
    url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/SMILES/JSON"%urllib.parse.quote(name,safe='')
    try:
        out=subprocess.run(["curl","-sS","-L","--max-time","25",url],capture_output=True,text=True).stdout
        return json.loads(out)['PropertyTable']['Properties'][0].get('SMILES')
    except Exception: return None
names=[l.strip() for l in open('dm_names.txt') if l.strip()]
rows=[]; nf=0
for i,name in enumerate(names):
    q=name.replace("''","'")           # normalize prime encoding
    smi=pubchem(q)
    if not smi and "'" in q: smi=pubchem(q.replace("'",""))
    ik=conn=None
    if smi: ik,conn=key(smi)
    if not ik: nf+=1
    rows.append({'name':name,'smiles':smi or '','inchikey':ik or '','connectivity':conn or ''})
    time.sleep(0.1)
    if (i+1)%100==0: print(f"  ...{i+1}/{len(names)} done, {nf} unresolved so far", flush=True)
with open('dm_keyed.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['name','smiles','inchikey','connectivity']); w.writeheader(); w.writerows(rows)
print('resolved',sum(1 for r in rows if r['inchikey']),'/',len(rows),' unresolved',nf)
