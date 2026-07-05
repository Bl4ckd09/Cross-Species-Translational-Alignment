#!/usr/bin/env python3
"""
structure_embed.py — SWAPPABLE structure encoder (SMILES -> vector).

Two interchangeable backends behind ONE contract:  featurize(connectivities, kind=...)

  kind="ecfp4"   (default) : ECFP4 Morgan r=2, 2048-bit, RDKit. No torch. Standard, strong
                             Tox21 structure baseline. Used on the Intel-Mac dev box.
  kind="chembert"          : ChemBERTa embeddings (768-d, mean-pooled), frozen. This is the
                             Baseline-1 encoder the plan specifies. Needs torch+transformers.

--------------------------------------------------------------------------------------------
TEAMMATES — to use ChemBERT instead of ECFP:
  1) install (on a torch-capable machine — Apple-Silicon/Linux/cloud):
        pip install "torch>=2.1" "transformers>=4.30" "rdkit"
     (Intel-mac note: torch caps at 2.2.2, so there pin  "transformers==4.44.2"  to match.)
  2) flip ONE line in scripts/run_experiment.py:   CFG["struct_kind"] = "chembert"
     (or run:  python -c "from structure_embed import featurize; featurize(conns, kind='chembert')")
  Both backends cache to data/features/, so you featurize once and reuse.
--------------------------------------------------------------------------------------------
"""
import os
import numpy as np
import pandas as pd

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "features"); os.makedirs(CACHE, exist_ok=True)

# ChemBERTa checkpoint. Alternatives: "DeepChem/ChemBERTa-77M-MLM", "DeepChem/ChemBERTa-10M-MTR".
CHEMBERT_MODEL = "seyonec/ChemBERTa-zinc-base-v1"


def _smiles_map():
    """connectivity -> SMILES, merged across DrugMatrix + TG-GATEs keyed tables."""
    m = {}
    for fn in ("dm_keyed.csv", "tggates_keyed.csv"):
        p = os.path.join(ROOT, "data", fn)
        if os.path.exists(p):
            d = pd.read_csv(p, dtype=str).dropna(subset=["smiles"])
            for _, r in d.drop_duplicates("connectivity").iterrows():
                m.setdefault(r["connectivity"], r["smiles"])
    return m


# ------------------------------- ECFP4 backend --------------------------------
def _ecfp4(connectivities, c2smi, n_bits=2048, radius=2):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    rows = []
    for c in connectivities:
        smi = c2smi.get(c)
        mol = Chem.MolFromSmiles(smi) if smi else None
        rows.append(np.zeros(n_bits, "float32") if mol is None
                    else gen.GetFingerprintAsNumPy(mol).astype("float32"))
    return np.vstack(rows), [f"ecfp_{i}" for i in range(n_bits)]


# --------------------------- ECFP-counts + physchem ---------------------------
# The stronger classical Tox21 structure baseline: substructure *counts* (not just
# presence/absence) concatenated with interpretable physicochemical descriptors.
# Pairs naturally with an L1 head that selects the handful of bits/descriptors that
# actually carry each assay. Standardisation is left to the downstream (in-fold) scaler.
_PHYSCHEM = [
    ("MolWt",            "MolWt"),
    ("MolLogP",          "MolLogP"),
    ("TPSA",             "TPSA"),
    ("NumHDonors",       "NumHDonors"),
    ("NumHAcceptors",    "NumHAcceptors"),
    ("NumRotatableBonds","NumRotatableBonds"),
    ("NumAromaticRings", "NumAromaticRings"),
    ("FractionCSP3",     "FractionCSP3"),
    ("NumHeteroatoms",   "NumHeteroatoms"),
    ("RingCount",        "RingCount"),
    ("NumSaturatedRings","NumSaturatedRings"),
    ("NumAliphaticRings","NumAliphaticRings"),
    ("HeavyAtomCount",   "HeavyAtomCount"),
    ("NHOHCount",        "NHOHCount"),
    ("NOCount",          "NOCount"),
    ("LabuteASA",        "LabuteASA"),
    ("BertzCT",          "BertzCT"),
    ("qed",              "qed"),
]

def _ecfp_counts(connectivities, c2smi, n_bits=2048, radius=2):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    descs = {name: getattr(Descriptors, fn) for name, fn in _PHYSCHEM}
    fp_rows, ph_rows = [], []
    for c in connectivities:
        smi = c2smi.get(c)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            fp_rows.append(np.zeros(n_bits, "float32"))
            ph_rows.append(np.zeros(len(_PHYSCHEM), "float32"))
            continue
        fp_rows.append(gen.GetCountFingerprintAsNumPy(mol).astype("float32"))
        vals = []
        for name in descs:
            try:
                v = float(descs[name](mol))
            except Exception:
                v = 0.0
            vals.append(0.0 if (v != v or v in (float("inf"), float("-inf"))) else v)
        ph_rows.append(np.array(vals, "float32"))
    mat = np.hstack([np.vstack(fp_rows), np.vstack(ph_rows)])
    cols = [f"cnt_{i}" for i in range(n_bits)] + [f"phys_{name}" for name, _ in _PHYSCHEM]
    return mat, cols


# ------------------------------ ChemBERT backend ------------------------------
def _chembert(connectivities, c2smi, model_name=CHEMBERT_MODEL, batch_size=64, max_len=256):
    """Frozen ChemBERTa mean-pooled embeddings. Requires torch + transformers."""
    import torch
    from transformers import AutoTokenizer, AutoModel
    from rdkit import Chem

    device = ("cuda" if torch.cuda.is_available()
              else "mps" if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available()
              else "cpu")
    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(device).eval()

    # canonicalise SMILES; missing -> None (filled with zeros afterwards)
    smis, keep = [], []
    for c in connectivities:
        s = c2smi.get(c)
        m = Chem.MolFromSmiles(s) if s else None
        smis.append(Chem.MolToSmiles(m) if m is not None else None)
        keep.append(m is not None)

    dim = model.config.hidden_size
    out = np.zeros((len(connectivities), dim), "float32")
    idx = [i for i, k in enumerate(keep) if k]
    with torch.no_grad():
        for b in range(0, len(idx), batch_size):
            sel = idx[b:b + batch_size]
            enc = tok([smis[i] for i in sel], padding=True, truncation=True,
                      max_length=max_len, return_tensors="pt").to(device)
            h = model(**enc).last_hidden_state                      # (B, T, H)
            mask = enc["attention_mask"].unsqueeze(-1).float()      # (B, T, 1)
            pooled = (h * mask).sum(1) / mask.sum(1).clamp(min=1)   # mean-pool over real tokens
            out[sel] = pooled.cpu().numpy().astype("float32")
    n_missing = len(connectivities) - len(idx)
    if n_missing:
        print(f"[chembert] {n_missing} compounds had no SMILES -> zero vector")
    return out, [f"cb_{i}" for i in range(dim)]


# --------------------------------- dispatch -----------------------------------
def featurize(connectivities, kind="ecfp4", **kw):
    connectivities = list(connectivities)
    tag = {"ecfp4": "ecfp4_2048", "ecfp_counts": "ecfp_counts_2048_physchem",
           "chembert": CHEMBERT_MODEL.split("/")[-1]}[kind]
    cache = os.path.join(CACHE, f"structure_{tag}.parquet")
    if os.path.exists(cache):
        df = pd.read_parquet(cache)
        if set(connectivities).issubset(df.index):
            return df.loc[connectivities]

    c2smi = _smiles_map()
    if kind == "ecfp4":
        mat, cols = _ecfp4(connectivities, c2smi, **kw)
    elif kind == "ecfp_counts":
        mat, cols = _ecfp_counts(connectivities, c2smi, **kw)
    elif kind == "chembert":
        mat, cols = _chembert(connectivities, c2smi, **kw)
    else:
        raise ValueError(f"unknown structure encoder: {kind}")

    df = pd.DataFrame(mat, index=connectivities, columns=cols)
    df.index.name = "connectivity"
    df.to_parquet(cache)
    return df


if __name__ == "__main__":
    import sys
    kind = sys.argv[1] if len(sys.argv) > 1 else "ecfp4"
    sig = pd.read_parquet(os.path.join(ROOT, "data", "signatures", "drugmatrix_liver_logfc.parquet"))
    f = featurize(list(sig.index), kind=kind)
    print(f"[{kind}] structure features:", f.shape)
