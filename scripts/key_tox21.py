import csv, sys
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

src="../tox21.csv"      # MoleculeNet Tox21 (SMILES + 12 assay columns)
out="tox21_keyed.csv"   # run from data/

lfc = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()

def standardize(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    m = rdMolStandardize.Cleanup(m)
    m = lfc.choose(m)          # strip salts/counterions
    m = uncharger.uncharge(m)  # neutralize
    return m

n=0; ok=0; fail=0
rows=[]
with open(src) as f:
    r=csv.DictReader(f)
    assays=[c for c in r.fieldnames if c!="SMILES"]
    for row in r:
        n+=1
        m=standardize(row["SMILES"])
        if m is None:
            fail+=1; continue
        try:
            ik=inchi.InchiToInchiKey(inchi.MolToInchi(m))
        except Exception:
            fail+=1; continue
        if not ik:
            fail+=1; continue
        ok+=1
        rows.append((ik, ik.split("-")[0], Chem.MolToSmiles(m), row))

# dedupe by full InChIKey, merging assay calls
from collections import defaultdict
byik={}
for ik, conn, csmi, row in rows:
    byik.setdefault(ik, (conn, csmi, []))[2].append(row)

with open(out,"w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["inchikey","connectivity","std_smiles"]+assays)
    for ik,(conn,csmi,rws) in byik.items():
        merged=[]
        for a in assays:
            vals=[rw[a] for rw in rws if rw[a] not in ("", None)]
            merged.append(vals[0] if vals else "")
        w.writerow([ik,conn,csmi]+merged)

print(f"input_rows={n} keyed_ok={ok} failed={fail}")
print(f"unique_full_inchikey={len(byik)}")
print(f"unique_connectivity={len(set(conn for _,(conn,_,_) in byik.items()))}")
