#!/usr/bin/env python3
"""
run_structure_rich.py — does a STRONGER structure baseline move the needle?

Two worries about the headline result (fusion 0.766 vs ECFP4-binary structure 0.757):
  (a) was the structure arm just under-tuned?  A richer representation (substructure
      *counts* + physicochemical descriptors) with an L1 head that selects the few
      bits/descriptors that matter is the classical *strong* Tox21 structure baseline.
  (b) if structure is stronger, does gene expression still add — and still on SR?

So we run, in the identical leakage-safe repeated-stratified-CV harness, five heads:

  struct_pca_L2   L2 on richer-structure -> PCA-128        (drop-in for the 0.757 baseline;
                                                            same dim as GE-100 -> fair fusion)
  struct_raw_L1   L1 on standardised RAW richer structure  (best-shot ceiling probe; L1 does
                                                            the feature selection PCA can't)
  expr_L2         L2 on GE -> PCA-100
  fusion_pca_L2   L2 on [richer-PCA-128 , GE-100]          (GE-adds question, comparable)
  fusion_raw_L1   L1 on [raw richer , GE-100]              (does GE add on the STRONGEST struct?)

Everything (both PCAs, all scalers) is fit on the TRAIN FOLD ONLY; heads train on measured
labels only; pooled out-of-fold predictions; fixed seeds. Same folds/repeats/seed as the
main experiment so numbers are directly comparable.

Writes: data/results/structure_rich.csv (+ .txt summary with SR/NR split and deltas).
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score

import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure_embed import featurize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIG  = os.path.join(ROOT, "data", "signatures")
RES  = os.path.join(ROOT, "data", "results"); os.makedirs(RES, exist_ok=True)

ASSAYS = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma",
          "SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]
SR = ["SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]
NR = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma"]

CFG = dict(n_pca_expr=100, n_pca_struct=128, folds=5, repeats=10, C=1.0, seed0=1000,
           signatures="drugmatrix_liver_logfc.parquet", labels="labels.csv")
APPROACHES = ["struct_pca_L2", "struct_raw_L1", "expr_L2", "fusion_pca_L2", "fusion_raw_L1"]


def load(cfg=CFG):
    sig = pd.read_parquet(os.path.join(SIG, cfg["signatures"]))
    lab = pd.read_csv(os.path.join(SIG, cfg["labels"])).set_index("connectivity").loc[sig.index]
    st  = featurize(list(sig.index), kind="ecfp_counts").loc[sig.index]
    Y = lab[ASSAYS].apply(pd.to_numeric, errors="coerce").values.astype(float)
    return sig.index.to_numpy(), sig.values.astype("float32"), st.values.astype("float32"), Y


def head(penalty, C):
    # sklearn>=1.8: penalty is set via l1_ratio (1.0=L1, 0.0=L2). True L1 needs the saga solver.
    if penalty == "l1":
        return LogisticRegression(max_iter=8000, class_weight="balanced", C=C,
                                  l1_ratio=1.0, solver="saga", tol=1e-3)
    return LogisticRegression(max_iter=3000, class_weight="balanced", C=C)  # L2 (lbfgs default)


def fit_masked(mk_clf, Xtr, ytr):
    m = ~np.isnan(ytr)
    if ytr[m].sum() < 2 or (1 - ytr[m]).sum() < 2:
        return None
    clf = mk_clf()
    clf.fit(Xtr[m], ytr[m].astype(int))
    return clf


def evaluate(cfg, data):
    conns, GE, ST, Y = data
    N = len(conns)
    strat = pd.qcut(np.nansum(Y == 1, axis=1), 4, labels=False, duplicates="drop")
    rec = {ap: {a: {"auc": [], "ap": []} for a in ASSAYS} for ap in APPROACHES}
    for r in range(cfg["repeats"]):
        skf = StratifiedKFold(cfg["folds"], shuffle=True, random_state=cfg["seed0"] + r)
        oof = {ap: np.full((N, len(ASSAYS)), np.nan) for ap in APPROACHES}
        for tr, te in skf.split(np.zeros(N), strat):
            # --- GE block: standardise+PCA(100) in-fold, then z-score PCs ---
            gsc = StandardScaler().fit(GE[tr])
            gp  = PCA(n_components=min(cfg["n_pca_expr"], len(tr) - 1), random_state=0)
            Zge_tr = gp.fit_transform(gsc.transform(GE[tr])); Zge_te = gp.transform(gsc.transform(GE[te]))
            gz = StandardScaler().fit(Zge_tr)
            Zge_tr, Zge_te = gz.transform(Zge_tr), gz.transform(Zge_te)

            # --- structure RAW: standardise in-fold (no PCA) -> L1 works on real bits ---
            ssc = StandardScaler().fit(ST[tr])
            Sraw_tr, Sraw_te = ssc.transform(ST[tr]), ssc.transform(ST[te])

            # --- structure PCA-128: standardise+PCA in-fold, then z-score PCs ---
            sp = PCA(n_components=min(cfg["n_pca_struct"], len(tr) - 1), random_state=0)
            Zst_tr = sp.fit_transform(Sraw_tr); Zst_te = sp.transform(Sraw_te)
            sz = StandardScaler().fit(Zst_tr)
            Zst_tr, Zst_te = sz.transform(Zst_tr), sz.transform(Zst_te)

            feats = {
                "struct_pca_L2": ("l2", Zst_tr, Zst_te),
                "struct_raw_L1": ("l1", Sraw_tr, Sraw_te),
                "expr_L2":       ("l2", Zge_tr, Zge_te),
                "fusion_pca_L2": ("l2", np.hstack([Zst_tr, Zge_tr]), np.hstack([Zst_te, Zge_te])),
                "fusion_raw_L1": ("l1", np.hstack([Sraw_tr, Zge_tr]), np.hstack([Sraw_te, Zge_te])),
            }
            for ap, (pen, Ftr, Fte) in feats.items():
                for j, a in enumerate(ASSAYS):
                    clf = fit_masked(lambda pen=pen: head(pen, cfg["C"]), Ftr, Y[tr, j])
                    if clf is None: continue
                    mte = ~np.isnan(Y[te, j])
                    if mte.any():
                        oof[ap][te[mte], j] = clf.predict_proba(Fte[mte])[:, 1]
        for ap in APPROACHES:
            for j, a in enumerate(ASSAYS):
                yj, pj = Y[:, j], oof[ap][:, j]
                m = ~np.isnan(yj) & ~np.isnan(pj)
                if yj[m].sum() >= 2 and (1 - yj[m]).sum() >= 2:
                    rec[ap][a]["auc"].append(roc_auc_score(yj[m], pj[m]))
                    rec[ap][a]["ap"].append(average_precision_score(yj[m], pj[m]))
    return rec


def macro(rec, ap, assays=ASSAYS):
    return float(np.mean([np.mean(rec[ap][a]["auc"]) for a in assays]))


def run():
    data = load(CFG)
    print("loaded:", len(data[0]), "compounds |", data[2].shape[1], "structure features "
          "(2048 counts + physchem)\n")
    rec = evaluate(CFG, data)

    def ci(x):
        x = np.array(x); return (x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else (np.nan, np.nan)

    rows = []
    for j, a in enumerate(ASSAYS):
        row = {"assay": a, "panel": "SR" if a in SR else "NR",
               "n_active": int(np.nansum(data[3][:, j] == 1))}
        for ap in APPROACHES:
            m, e = ci(rec[ap][a]["auc"]); row[f"{ap}_AUC"] = round(m, 3)
        rows.append(row)
    tab = pd.DataFrame(rows)

    # macro + panel rows
    def add_summary(name, assays):
        r = {"assay": name, "panel": "", "n_active": ""}
        for ap in APPROACHES:
            r[f"{ap}_AUC"] = round(macro(rec, ap, assays), 3)
        return r
    tab = pd.concat([tab,
                     pd.DataFrame([add_summary("MACRO", ASSAYS),
                                   add_summary("SR_mean", SR),
                                   add_summary("NR_mean", NR)])], ignore_index=True)
    tab.to_csv(os.path.join(RES, "structure_rich.csv"), index=False)

    # deltas (per repeat, paired)
    def dstats(ap_a, ap_b, assays):
        ds = []
        for a in assays:
            x = np.array(rec[ap_a][a]["auc"]); y = np.array(rec[ap_b][a]["auc"])
            n = min(len(x), len(y)); ds.append((x[:n] - y[:n]).mean())
        return np.mean(ds)

    lines = [
        "STRUCTURE-RICH probe  (N=177 DrugMatrix liver, %dx%d CV, C=%.1f)" % (CFG["repeats"], CFG["folds"], CFG["C"]),
        "features: 2048 substructure COUNTS + 18 physchem descriptors",
        "",
        "MACRO ROC-AUC:",
        *[f"  {ap:16s} {macro(rec, ap):.3f}" for ap in APPROACHES],
        "",
        "reference (binary ECFP4 -> PCA-128, L2):  structure 0.757  fusion 0.766",
        "",
        "== Q(a) did a richer structure baseline raise the ceiling? ==",
        f"  struct_raw_L1 - struct_pca_L2  (macro) = {macro(rec,'struct_raw_L1')-macro(rec,'struct_pca_L2'):+.3f}",
        f"  best structure macro AUC              = {max(macro(rec,'struct_pca_L2'), macro(rec,'struct_raw_L1')):.3f}",
        "",
        "== Q(b) does GE still add — and still on SR? ==",
        f"  fusion_pca_L2 - struct_pca_L2   macro={dstats('fusion_pca_L2','struct_pca_L2',ASSAYS):+.3f}  "
        f"SR={dstats('fusion_pca_L2','struct_pca_L2',SR):+.3f}  NR={dstats('fusion_pca_L2','struct_pca_L2',NR):+.3f}",
        f"  fusion_raw_L1 - struct_raw_L1   macro={dstats('fusion_raw_L1','struct_raw_L1',ASSAYS):+.3f}  "
        f"SR={dstats('fusion_raw_L1','struct_raw_L1',SR):+.3f}  NR={dstats('fusion_raw_L1','struct_raw_L1',NR):+.3f}",
    ]
    txt = "\n".join(lines)
    with open(os.path.join(RES, "structure_rich.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(RES, "structure_rich_config.json"), "w") as f:
        json.dump({"config": CFG, "n_compounds": int(len(data[0])),
                   "n_struct_features": int(data[2].shape[1]), "assays": ASSAYS}, f, indent=2)

    pd.set_option("display.width", 220, "display.max_columns", 30)
    print(tab.to_string(index=False))
    print("\n" + txt)
    print(f"\nwrote {RES}/structure_rich.csv  and  structure_rich.txt")


if __name__ == "__main__":
    run()
