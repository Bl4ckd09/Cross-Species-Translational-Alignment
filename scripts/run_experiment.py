#!/usr/bin/env python3
"""
run_experiment.py — the controlled comparison: does gene expression add to structure?

Approaches (identical everything else, isolating the effect of the GE block):
  structure_only  (Baseline 1) : head on structure features (ECFP4 -> PCA-128)
  expr_only       (ablation)   : head on GE features        (logFC  -> PCA-100)
  fusion          (Approach A) : head on [structure-128 , GE-100]

Leakage protocol (followed exactly):
  - repeated STRATIFIED K-fold (stratified on per-compound active-count); each compound is a
    test compound exactly once per repeat -> pooled out-of-fold (OOF) predictions.
  - every data-dependent transform (both PCAs + their scalers) is fit on the TRAIN FOLD ONLY.
  - per-assay heads train on MEASURED labels only (masking); predict only measured test labels.
  - fixed seeds; splits persisted. Metrics: per-assay ROC-AUC + AUPRC, macro, with CIs across
    repeats, plus paired delta (fusion - structure).

Note: single-source (DrugMatrix liver, N=177) first pass -> no ComBat needed yet.
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score

import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from expr_embed import make_embedder
from structure_embed import featurize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIG  = os.path.join(ROOT, "data", "signatures")
RES  = os.path.join(ROOT, "data", "results"); os.makedirs(RES, exist_ok=True)
ASSAYS = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma",
          "SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]

CFG = dict(n_pca_expr=100, n_pca_struct=128, folds=5, repeats=10, C=1.0, seed0=1000)

def load():
    sig = pd.read_parquet(os.path.join(SIG, "drugmatrix_liver_logfc.parquet"))
    lab = pd.read_csv(os.path.join(SIG, "labels.csv")).set_index("connectivity").loc[sig.index]
    st  = featurize(list(sig.index)).loc[sig.index]
    Y = lab[ASSAYS].apply(pd.to_numeric, errors="coerce").values.astype(float)  # 0/1/nan
    return sig.index.to_numpy(), sig.values.astype("float32"), st.values.astype("float32"), Y

def fit_head(Xtr, ytr, C):
    m = ~np.isnan(ytr)
    if ytr[m].sum() < 2 or (1 - ytr[m]).sum() < 2:
        return None
    clf = LogisticRegression(max_iter=3000, class_weight="balanced", C=C)
    clf.fit(Xtr[m], ytr[m].astype(int))
    return clf

def run():
    conns, GE, ST, Y = load()
    N = len(conns)
    active_count = np.nansum(Y == 1, axis=1)
    strat = pd.qcut(active_count, 4, labels=False, duplicates="drop")
    approaches = ["structure_only", "expr_only", "fusion"]

    # per-repeat, per-approach, per-assay AUC/AUPRC
    rec = {ap: {a: {"auc": [], "ap": []} for a in ASSAYS} for ap in approaches}
    splits_log = []

    for r in range(CFG["repeats"]):
        seed = CFG["seed0"] + r
        skf = StratifiedKFold(CFG["folds"], shuffle=True, random_state=seed)
        oof = {ap: np.full((N, len(ASSAYS)), np.nan) for ap in approaches}
        for tr, te in skf.split(np.zeros(N), strat):
            splits_log.append({"repeat": r, "seed": seed, "n_train": len(tr), "n_test": len(te)})
            ge = make_embedder("pca", n_components=CFG["n_pca_expr"])
            Zge_tr, Zge_te = ge.fit_transform(GE[tr]), ge.transform(GE[te])
            se = make_embedder("pca", n_components=CFG["n_pca_struct"])
            Zst_tr, Zst_te = se.fit_transform(ST[tr]), se.transform(ST[te])
            feats = {
                "structure_only": (Zst_tr, Zst_te),
                "expr_only":      (Zge_tr, Zge_te),
                "fusion":         (np.hstack([Zst_tr, Zge_tr]), np.hstack([Zst_te, Zge_te])),
            }
            for ap, (Ftr, Fte) in feats.items():
                # z-score all PCs (fit on train) so structure & GE blocks are balanced
                sc = StandardScaler().fit(Ftr)
                Ftr, Fte = sc.transform(Ftr), sc.transform(Fte)
                for j, a in enumerate(ASSAYS):
                    clf = fit_head(Ftr, Y[tr, j], CFG["C"])
                    if clf is None: continue
                    mte = ~np.isnan(Y[te, j])
                    if mte.any():
                        oof[ap][te[mte], j] = clf.predict_proba(Fte[mte])[:, 1]
        # per-repeat metrics from pooled OOF
        for ap in approaches:
            for j, a in enumerate(ASSAYS):
                yj, pj = Y[:, j], oof[ap][:, j]
                m = ~np.isnan(yj) & ~np.isnan(pj)
                if yj[m].sum() >= 2 and (1 - yj[m]).sum() >= 2:
                    rec[ap][a]["auc"].append(roc_auc_score(yj[m], pj[m]))
                    rec[ap][a]["ap"].append(average_precision_score(yj[m], pj[m]))

    # ---- assemble results table ----
    def ci(x):
        x = np.array(x); return (x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else (np.nan, np.nan)
    prevalence = {a: np.nanmean(Y[:, j] == 1) for j, a in enumerate(ASSAYS)}
    rows = []
    for j, a in enumerate(ASSAYS):
        row = {"assay": a, "n_active": int(np.nansum(Y[:, j] == 1)), "prevalence": round(prevalence[a], 3)}
        for ap in approaches:
            am, ae = ci(rec[ap][a]["auc"]); pm, pe = ci(rec[ap][a]["ap"])
            row[f"{ap}_AUC"] = round(am, 3); row[f"{ap}_AUC_ci"] = round(ae, 3)
            row[f"{ap}_AUPRC"] = round(pm, 3)
        # paired delta fusion - structure (per repeat)
        d = np.array(rec["fusion"][a]["auc"]) - np.array(rec["structure_only"][a]["auc"])
        row["dAUC_fusion_minus_struct"] = round(d.mean(), 3)
        row["dAUC_ci"] = round(1.96 * d.std(ddof=1) / np.sqrt(len(d)), 3) if len(d) > 1 else np.nan
        row["fusion_wins_frac"] = round((d > 0).mean(), 2)
        rows.append(row)
    tab = pd.DataFrame(rows)
    # macro row
    macro = {"assay": "MACRO", "n_active": tab.n_active.sum()}
    for ap in approaches:
        macro[f"{ap}_AUC"] = round(tab[f"{ap}_AUC"].mean(), 3)
        macro[f"{ap}_AUPRC"] = round(tab[f"{ap}_AUPRC"].mean(), 3)
    macro["dAUC_fusion_minus_struct"] = round(tab.dAUC_fusion_minus_struct.mean(), 3)
    tab = pd.concat([tab, pd.DataFrame([macro])], ignore_index=True)

    tab.to_csv(os.path.join(RES, "results_table.csv"), index=False)
    with open(os.path.join(RES, "run_config.json"), "w") as f:
        json.dump({"config": CFG, "n_compounds": int(N), "assays": ASSAYS,
                   "n_splits_total": len(splits_log)}, f, indent=2)

    pd.set_option("display.width", 200, "display.max_columns", 30)
    show = ["assay", "n_active", "prevalence", "structure_only_AUC", "expr_only_AUC",
            "fusion_AUC", "dAUC_fusion_minus_struct", "dAUC_ci", "fusion_wins_frac"]
    print("\n=== RESULTS: per-assay ROC-AUC (OOF, mean over %d repeats x %d folds) ===" % (CFG["repeats"], CFG["folds"]))
    print(tab[show].to_string(index=False))
    print(f"\nMACRO AUC  structure={macro['structure_only_AUC']}  expr_only={macro['expr_only_AUC']}  "
          f"fusion={macro['fusion_AUC']}  |  mean dAUC(fusion-struct)={macro['dAUC_fusion_minus_struct']}")
    print(f"wrote {RES}/results_table.csv")

if __name__ == "__main__":
    run()
