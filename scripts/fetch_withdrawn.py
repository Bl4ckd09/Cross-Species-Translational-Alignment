#!/usr/bin/env python3
"""
fetch_withdrawn.py — secondary DILI target: market-withdrawal status (ChEMBL drug_warning).

Broader, noisier supplementary signal — kept SEPARATE from DILIrank (different label, own file).
Source: ChEMBL API `drug_warning?warning_type=Withdrawn`. Matched to the rat cohort via InChIKey
connectivity (same pipeline as the main build) — not fuzzy names.

Output: data/dili/withdrawn_chembl.csv  (chembl_id, name, inchikey, connectivity, withdrawn_year, reason)
"""
import os, time
import requests
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://www.ebi.ac.uk/chembl/api/data"

def get_withdrawn_ids():
    ids, reasons, years, offset = {}, {}, {}, 0
    while True:
        r = requests.get(f"{API}/drug_warning.json",
                         params={"warning_type": "Withdrawn", "limit": 1000, "offset": offset}, timeout=60)
        r.raise_for_status(); d = r.json()
        for w in d["drug_warnings"]:
            cid = w.get("parent_molecule_chembl_id") or w.get("molecule_chembl_id")
            if not cid: continue
            ids.setdefault(cid, True)
            if w.get("warning_description"): reasons.setdefault(cid, w["warning_description"][:120])
            if w.get("warning_year"): years.setdefault(cid, w["warning_year"])
        offset += 1000
        if offset >= d["page_meta"]["total_count"]: break
    return ids, reasons, years

def resolve_molecules(chembl_ids):
    rows, ids = [], list(chembl_ids)
    for i in range(0, len(ids), 50):
        batch = ids[i:i + 50]
        r = requests.get(f"{API}/molecule.json",
                         params={"molecule_chembl_id__in": ",".join(batch), "limit": 50}, timeout=60)
        r.raise_for_status()
        for m in r.json()["molecules"]:
            struct = m.get("molecule_structures") or {}
            ik = struct.get("standard_inchi_key")
            rows.append({"chembl_id": m["molecule_chembl_id"], "name": m.get("pref_name"),
                         "inchikey": ik, "connectivity": ik[:14] if ik else None})
        time.sleep(0.1)
    return pd.DataFrame(rows)

def main():
    os.makedirs(os.path.join(ROOT, "data", "dili"), exist_ok=True)
    print("fetching withdrawn drug_warnings from ChEMBL ...")
    ids, reasons, years = get_withdrawn_ids()
    print(f"unique withdrawn molecules: {len(ids)}")
    df = resolve_molecules(ids)
    df["withdrawn_year"] = df.chembl_id.map(years); df["reason"] = df.chembl_id.map(reasons)
    df = df.dropna(subset=["connectivity"]).drop_duplicates("connectivity")
    out = os.path.join(ROOT, "data", "dili", "withdrawn_chembl.csv")
    df.to_csv(out, index=False)
    print(f"wrote {out}: {len(df)} withdrawn drugs with structures")

    coh = pd.read_csv(os.path.join(ROOT, "master_cohort.csv"), dtype=str)
    comb = set(pd.read_parquet(os.path.join(ROOT, "data", "signatures", "combined_logfc.parquet")).index)
    wset = set(df.connectivity)
    ov = set(coh.connectivity) & wset
    print(f"\noverlap with rat cohort (672): {len(ov)} withdrawn compounds")
    print(f"  ...of which have expression (256): {len(ov & comb)}")

if __name__ == "__main__":
    main()
