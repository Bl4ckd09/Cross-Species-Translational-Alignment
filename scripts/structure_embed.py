#!/usr/bin/env python3
"""
structure_embed.py — SWAPPABLE structure encoder (SMILES -> vector).

Default: ECFP4 (Morgan r=2, 2048-bit) via RDKit — the standard, strong Tox21 structure
baseline, no torch needed. ChemBERT is a drop-in alternative behind the same featurize()
contract (blocked here only by this Intel-Mac's torch 2.2 ceiling; enable when available).

    feats = featurize(connectivities)   # -> DataFrame [connectivity x 2048], cached
"""
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "features"); os.makedirs(CACHE, exist_ok=True)

def _smiles_map():
    dmk = pd.read_csv(os.path.join(ROOT, "data", "dm_keyed.csv"), dtype=str).dropna(subset=["smiles"])
    return dmk.drop_duplicates("connectivity").set_index("connectivity")["smiles"].to_dict()

def featurize(connectivities, kind="ecfp4", n_bits=2048, radius=2):
    cache = os.path.join(CACHE, f"structure_{kind}_{n_bits}.parquet")
    if os.path.exists(cache):
        df = pd.read_parquet(cache)
        if set(connectivities).issubset(df.index):
            return df.loc[list(connectivities)]
    c2smi = _smiles_map()
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    rows, idx = [], []
    for c in connectivities:
        smi = c2smi.get(c)
        m = Chem.MolFromSmiles(smi) if smi else None
        if m is None:
            rows.append(np.zeros(n_bits, dtype="float32"))   # unresolved -> zero vector
        else:
            rows.append(gen.GetFingerprintAsNumPy(m).astype("float32"))
        idx.append(c)
    df = pd.DataFrame(np.vstack(rows), index=idx,
                      columns=[f"ecfp_{i}" for i in range(n_bits)])
    df.index.name = "connectivity"
    df.to_parquet(cache)
    return df

if __name__ == "__main__":
    sig = pd.read_parquet(os.path.join(ROOT, "data", "signatures", "drugmatrix_liver_logfc.parquet"))
    f = featurize(list(sig.index))
    print("structure features:", f.shape, "| nonzero bits/mol median:",
          int(np.median((f.values > 0).sum(axis=1))))
