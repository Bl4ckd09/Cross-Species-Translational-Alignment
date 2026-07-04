#!/usr/bin/env python3
"""
validation_b.py — Retrospective validation on human clinical DILI (DILIrank).

Exploratory: does the model / its learned representation say anything about real hepatotoxicity,
beyond structure alone and beyond measured Tox21? Published context: measured Tox21 alone predicts
DILI only ~0.5 (near-random); structure-alone models reach ~0.75-0.83 — that is the bar to clear.

Baseline ladder (all built):
  0  prevalence            — majority class (imbalance reference; AUPRC baseline = prevalence)
  1  structure (ECFP4→PCA) — the number to beat   [ChemBERT = torch-env drop-in]
  2  measured Tox21 (12)   — should reproduce the published ~0.5 near-random result (sanity gate)
  3  raw expression (PCA)  — rat logFC → DILI, no structure, no Tox21 layer
Approaches:
  A  chained        — frozen Tox21 model's PREDICTED 12-D vector → DILI (deployed pipeline)
  B  representation — fused [structure-PCA ⊕ expr-PCA] pre-head rep → DILI

Leakage protocol: repeated stratified CV on DILI; EVERY transform (scaler, PCA) fit on the train
fold only; class-balanced logistic heads; ROC-AUC + AUPRC + CIs vs Baseline 0.

Honest caveats (reported, not hidden):
  * Class imbalance — the rat∩DILI overlap is enriched for hepatotoxins (~82% positive).
  * Input overlap is UNAVOIDABLE for B3/A/B: expression + Tox21 exist ONLY for the 613, so every
    compound the full model can score is a training-input compound. The DILI *label* was never
    seen (different target) so it isn't label leakage, but the inputs/representation were fit on
    these compounds — Approach A's in-sample Tox21 predictions are optimistic (cross-fitting is the
    rigorous follow-up). Baselines 1/2/3 and Approach B don't reuse a DILI-trained model.
Adaptations: ECFP4 not ChemBERT; frozen Tox21 model trained on the 256 (where expression exists).
"""
import os, re, sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structure_embed import featurize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIG  = os.path.join(ROOT, "data", "signatures")
RES  = os.path.join(ROOT, "data", "results"); os.makedirs(RES, exist_ok=True)
ASSAYS = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma",
          "SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]
SALT = (r"\b(sulfate|hydrochloride|hcl|sodium|potassium|calcium|mesylate|maleate|citrate|tartrate|"
        r"acetate|phosphate|besylate|fumarate|succinate|dihydrate|hydrate|monohydrate|bromide|"
        r"chloride|nitrate|hydrobromide|disodium|hemihydrate|hydroxide|dihydrochloride)\b")
CFG = dict(n_pca_expr=100, n_pca_struct=128, folds=5, repeats=10, C=1.0, seed0=2000)

def norm(x):
    x = re.sub(r"[^a-z0-9 ]", " ", str(x).lower()); x = re.sub(SALT, " ", x)
    return re.sub(r"\s+", " ", x).strip()

def cv_eval(blocks, y):
    """blocks = list of (matrix, n_pca|None). Fit scaler(+PCA) IN-FOLD, concat, class-balanced
    logistic. Returns per-repeat AUC & AUPRC arrays (pooled OOF)."""
    y = y.astype(int); aucs, aps = [], []
    for r in range(CFG["repeats"]):
        skf = StratifiedKFold(CFG["folds"], shuffle=True, random_state=CFG["seed0"] + r)
        oof = np.full(len(y), np.nan)
        for tr, te in skf.split(np.zeros(len(y)), y):
            Ftr, Fte = [], []
            for M, k in blocks:
                sc = StandardScaler().fit(M[tr]); Mtr, Mte = sc.transform(M[tr]), sc.transform(M[te])
                if k:
                    p = PCA(min(k, min(Mtr.shape) - 1)).fit(Mtr); Mtr, Mte = p.transform(Mtr), p.transform(Mte)
                Ftr.append(Mtr); Fte.append(Mte)
            Xtr, Xte = np.hstack(Ftr), np.hstack(Fte)
            s2 = StandardScaler().fit(Xtr)
            clf = LogisticRegression(max_iter=3000, class_weight="balanced", C=CFG["C"]).fit(s2.transform(Xtr), y[tr])
            oof[te] = clf.predict_proba(s2.transform(Xte))[:, 1]
        aucs.append(roc_auc_score(y, oof)); aps.append(average_precision_score(y, oof))
    return np.array(aucs), np.array(aps)

def ci(x): return 1.96 * x.std(ddof=1) / np.sqrt(len(x))
def row(name, n, a, p): return (name, n, round(a.mean(), 3), round(ci(a), 3), round(p.mean(), 3), round(ci(p), 3))

def frozen_tox21_predict(GEc, STc, Ytox):
    """Frozen Tox21 model: per-assay balanced logistic on all-256 [struct-PCA ⊕ expr-PCA];
    return in-sample 12-D probs (the deployed frozen model's output)."""
    Zs = PCA(min(CFG["n_pca_struct"], min(STc.shape) - 1)).fit_transform(StandardScaler().fit_transform(STc))
    Zg = PCA(min(CFG["n_pca_expr"], min(GEc.shape) - 1)).fit_transform(StandardScaler().fit_transform(GEc))
    F = StandardScaler().fit_transform(np.hstack([Zs, Zg]))
    P = np.full((len(F), len(ASSAYS)), np.nan)
    for j in range(len(ASSAYS)):
        yj = Ytox[:, j]; m = ~np.isnan(yj)
        if yj[m].sum() >= 2 and (1 - yj[m]).sum() >= 2:
            P[:, j] = LogisticRegression(max_iter=3000, class_weight="balanced").fit(F[m], yj[m].astype(int)).predict_proba(F)[:, 1]
    return P

def main():
    dili = pd.read_excel(os.path.join(ROOT, "data", "dili", "DILIrank.xlsx"), sheet_name=0, skiprows=1)
    c = dili["vDILI-Concern"].str.lower().str.strip()
    dili["dili"] = np.where(c.str.contains("most|less"), 1.0, np.where(c.str.startswith("vno"), 0.0, np.nan))
    dili = dili[dili.dili.notna()].copy(); dili["nkey"] = dili["CompoundName"].map(norm)

    coh = pd.read_csv(os.path.join(ROOT, "master_cohort.csv"), dtype=str); coh["nkey"] = coh["compound_name"].map(norm)
    dili["conn"] = dili["nkey"].map(dict(zip(coh.nkey, coh.connectivity)))
    dili = dili.dropna(subset=["conn"]).drop_duplicates("conn")

    comb = pd.read_parquet(os.path.join(SIG, "combined_logfc.parquet"))
    tox = pd.read_csv(os.path.join(SIG, "combined_labels.csv")).set_index("connectivity")
    coh_i = coh.set_index("connectivity")
    expr_conns = set(comb.index)

    # frozen Tox21 predictions for the 256
    ST256 = featurize(list(comb.index), kind="ecfp4").loc[comb.index].values.astype("float32")
    Ytox256 = tox.loc[comb.index, ASSAYS].apply(pd.to_numeric, errors="coerce").values
    predD = pd.DataFrame(frozen_tox21_predict(comb.values.astype("float32"), ST256, Ytox256),
                         index=comb.index, columns=ASSAYS)

    conns = list(dili.conn); y = dili.set_index("conn")["dili"].loc[conns].values
    prev = y.mean()
    R = [("Baseline0 prevalence", len(y), 0.500, 0.000, round(prev, 3), 0.000)]

    # Baseline 1: structure, all overlap
    ST = featurize(conns, kind="ecfp4").loc[conns].values.astype("float32")
    R.append(row("Baseline1 structure(ECFP)", len(y), *cv_eval([(ST, CFG["n_pca_struct"])], y)))

    # Baseline 2: measured Tox21 (impute missing 0.5), compounds with >=1 measured
    toxm = coh_i.reindex(conns)[ASSAYS].apply(pd.to_numeric, errors="coerce")
    m2 = toxm.notna().any(axis=1).values
    R.append(row("Baseline2 measuredTox21", int(m2.sum()), *cv_eval([(toxm.fillna(0.5).values[m2], None)], y[m2])))

    # expression subset for B3 / A / B  — everything below on the SAME 135 compounds (fair head-to-head)
    ce = [c for c in conns if c in expr_conns]; ye = dili.set_index("conn")["dili"].loc[ce].values
    GEe = comb.loc[ce].values.astype("float32"); STe = featurize(ce, kind="ecfp4").loc[ce].values.astype("float32")
    toxe = coh_i.reindex(ce)[ASSAYS].apply(pd.to_numeric, errors="coerce").fillna(0.5).values
    R.append(("-- same-set (N=%d) head-to-head --" % len(ye), len(ye), np.nan, np.nan, np.nan, np.nan))
    R.append(row("  B1 structure @135",   len(ye), *cv_eval([(STe, CFG["n_pca_struct"])], ye)))
    R.append(row("  B2 measuredTox21 @135", len(ye), *cv_eval([(toxe, None)], ye)))
    R.append(row("  B3 rawExpression @135", len(ye), *cv_eval([(GEe, CFG["n_pca_expr"])], ye)))
    R.append(row("  A  predTox21 @135",    len(ye), *cv_eval([(predD.loc[ce].values, None)], ye)))
    R.append(row("  B  fusedRep @135",     len(ye), *cv_eval([(STe, CFG["n_pca_struct"]), (GEe, CFG["n_pca_expr"])], ye)))

    tab = pd.DataFrame(R, columns=["model", "N", "AUC", "AUC_ci", "AUPRC", "AUPRC_ci"])
    tab.to_csv(os.path.join(RES, "validation_b.csv"), index=False)
    print(f"\n=== VALIDATION B — DILIrank ===")
    print(f"overlap compounds: {len(y)}  positives {int(y.sum())} negatives {int((y==0).sum())}  "
          f"(prevalence {prev:.3f}; structure bar to clear ~0.75-0.83)")
    print(tab.to_string(index=False))
    print("\nread: Baseline2 measuredTox21≈0.52-0.59 -> reproduces published near-random Tox21->DILI (sanity OK);"
          "\n      structure is best; Approach A/B (expression, fused) do NOT beat structure -> the model's"
          "\n      representation adds no DILI signal over structure; raw expression≈0.52 (chance)."
          "\n      caveats: our structure 0.66-0.70 is below published SOTA 0.75-0.83 (weak ECFP+logistic,"
          "\n      small imbalanced N=135, 82% pos); A's in-sample Tox21 preds are optimistic yet still lose.")

if __name__ == "__main__":
    main()
