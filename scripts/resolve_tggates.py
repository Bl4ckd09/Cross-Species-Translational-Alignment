"""Resolve Open TG-GATEs compound names -> standardized InChIKeys.

Input : open_tggates_main.csv  (compound list from the LSDB Archive)
Output: tggates_keyed.csv

Names are resolved to SMILES via PubChem PUG-REST (through `curl`, because the
sandbox proxy uses a cert Python's ssl does not trust), then run through the
SAME RDKit standardization pipeline used for Tox21 so InChIKeys are comparable.
Matching downstream is done on both the full InChIKey and the first 14-char
connectivity block (salt/stereo/tautomer-insensitive).
"""
import csv, json, subprocess, time, urllib.parse
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

lfc = rdMolStandardize.LargestFragmentChooser()   # strip salts/counterions
unch = rdMolStandardize.Uncharger()               # neutralize


def key(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    m = rdMolStandardize.Cleanup(m); m = lfc.choose(m); m = unch.uncharge(m)
    try:
        ik = inchi.InchiToInchiKey(inchi.MolToInchi(m))
    except Exception:
        return None, None
    return (ik, ik.split('-')[0]) if ik else (None, None)


def pubchem_smiles(name):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/SMILES/JSON"
           % urllib.parse.quote(name, safe=''))
    try:
        out = subprocess.run(["curl", "-sS", "-L", "--max-time", "30", url],
                             capture_output=True, text=True).stdout
        return json.loads(out)['PropertyTable']['Properties'][0].get('SMILES')
    except Exception:
        return None


# a few TG-GATEs names need better synonyms for PubChem
SYNONYMS = {'acetamidofluorene': '2-acetamidofluorene',
            'naphthyl isothiocyanate': '1-naphthyl isothiocyanate'}


def liver(r):  return bool(r['Rat - in vivo - Liver - Single'] or r['Rat - in vivo - Liver - Repeat'])
def kidney(r): return bool(r['Rat - in vivo - Kidney - Single'] or r['Rat - in vivo - Kidney - Repeat'])


def main():
    rows = list(csv.DictReader(open('open_tggates_main.csv', encoding='latin-1')))
    out, fails = [], []
    for r in rows:
        name = r['COMPOUND_NAME']
        smi = pubchem_smiles(name) or pubchem_smiles(SYNONYMS.get(name, name))
        ik = conn = None
        if smi:
            ik, conn = key(smi)
        if not ik:
            fails.append(name)
        out.append({'name': name, 'abbr': r['COMPOUND_ABBREVIATION'],
                    'liver': liver(r), 'kidney': kidney(r),
                    'smiles': smi or '', 'inchikey': ik or '', 'connectivity': conn or ''})
        time.sleep(0.12)
    with open('tggates_keyed.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['name', 'abbr', 'liver', 'kidney',
                                          'smiles', 'inchikey', 'connectivity'])
        w.writeheader(); w.writerows(out)
    print('resolved', sum(1 for o in out if o['inchikey']), '/', len(out))
    print('unresolved (biologics / mixtures expected):', fails)


if __name__ == '__main__':
    main()
