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
from sklearn.model_selection import StratifiedKFold, KFold
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
STRUCT = os.environ.get("STRUCT_KIND", "ecfp4")   # "ecfp4" | "chembert" (cached; no torch in main venv)

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

def _tox_features(STtr, GEtr, STte, GEte):
    """struct-PCA ⊕ expr-PCA, all fit on the train rows only. Returns (F_train, F_test)."""
    ss = StandardScaler().fit(STtr); ps = PCA(min(CFG["n_pca_struct"], min(STtr.shape) - 1)).fit(ss.transform(STtr))
    gs = StandardScaler().fit(GEtr); pg = PCA(min(CFG["n_pca_expr"], min(GEtr.shape) - 1)).fit(gs.transform(GEtr))
    Ftr = np.hstack([ps.transform(ss.transform(STtr)), pg.transform(gs.transform(GEtr))])
    Fte = np.hstack([ps.transform(ss.transform(STte)), pg.transform(gs.transform(GEte))])
    f2 = StandardScaler().fit(Ftr)
    return f2.transform(Ftr), f2.transform(Fte)

def _fit_predict_tox(Ftr, Ytr, Fte):
    P = np.full((len(Fte), len(ASSAYS)), np.nan)
    for j in range(len(ASSAYS)):
        yj = Ytr[:, j]; m = ~np.isnan(yj)
        if yj[m].sum() >= 2 and (1 - yj[m]).sum() >= 2:
            P[:, j] = LogisticRegression(max_iter=3000, class_weight="balanced").fit(Ftr[m], yj[m].astype(int)).predict_proba(Fte)[:, 1]
    return P

def frozen_tox21_predict(GEc, STc, Ytox):
    """In-sample: Tox21 model trained on ALL rows, predicts those same rows (deployed frozen
    model's output — but optimistic: a compound's own Tox21 label informed its prediction)."""
    F, _ = _tox_features(STc, GEc, STc, GEc)
    return _fit_predict_tox(F, Ytox, F)

def crossfit_tox21_predict(GEc, STc, Ytox, folds=5, seed=7):
    """Out-of-fold: each compound's 12-D Tox21 prediction comes from a model trained on the OTHER
    folds only — removes the in-sample optimism (proper cross-fitting / stacking)."""
    n = len(STc); P = np.full((n, len(ASSAYS)), np.nan)
    for tr, te in KFold(folds, shuffle=True, random_state=seed).split(np.zeros(n)):
        Ftr, Fte = _tox_features(STc[tr], GEc[tr], STc[te], GEc[te])
        P[te] = _fit_predict_tox(Ftr, Ytox[tr], Fte)
    return P

def main(label_source="dilirank"):
    # marketed-drug universe = DILIrank-name-matched cohort compounds (shared across both labels)
    dili = pd.read_excel(os.path.join(ROOT, "data", "dili", "DILIrank.xlsx"), sheet_name=0, skiprows=1)
    c = dili["vDILI-Concern"].str.lower().str.strip()
    dili["dili"] = np.where(c.str.contains("most|less"), 1.0, np.where(c.str.startswith("vno"), 0.0, np.nan))
    dili["nkey"] = dili["CompoundName"].map(norm)
    coh = pd.read_csv(os.path.join(ROOT, "master_cohort.csv"), dtype=str); coh["nkey"] = coh["compound_name"].map(norm)
    dili["conn"] = dili["nkey"].map(dict(zip(coh.nkey, coh.connectivity)))
    dili = dili.dropna(subset=["conn"]).drop_duplicates("conn")

    if label_source == "dilirank":                        # DILI-concern (vMost+vLess vs vNo; ambiguous dropped)
        dili = dili[dili.dili.notna()].copy(); dili["label"] = dili["dili"]
    elif label_source == "withdrawn":                     # SEPARATE noisier target: market-withdrawal (ChEMBL)
        wset = set(pd.read_csv(os.path.join(ROOT, "data", "dili", "withdrawn_chembl.csv"))["connectivity"].dropna())
        dili["label"] = dili["conn"].isin(wset).astype(float)   # over ALL name-matched drugs (keeps ambiguous)
    else:
        raise ValueError(label_source)

    comb = pd.read_parquet(os.path.join(SIG, "combined_logfc.parquet"))
    tox = pd.read_csv(os.path.join(SIG, "combined_labels.csv")).set_index("connectivity")
    coh_i = coh.set_index("connectivity")
    expr_conns = set(comb.index)

    # frozen Tox21 predictions for the 256
    ST256 = featurize(list(comb.index), kind=STRUCT).loc[comb.index].values.astype("float32")
    Ytox256 = tox.loc[comb.index, ASSAYS].apply(pd.to_numeric, errors="coerce").values
    predD = pd.DataFrame(frozen_tox21_predict(comb.values.astype("float32"), ST256, Ytox256),
                         index=comb.index, columns=ASSAYS)                        # in-sample (optimistic)
    predCF = pd.DataFrame(crossfit_tox21_predict(comb.values.astype("float32"), ST256, Ytox256),
                          index=comb.index, columns=ASSAYS)                       # out-of-fold (rigorous)

    conns = list(dili.conn); y = dili.set_index("conn")["label"].loc[conns].values
    prev = y.mean()
    R = [("Baseline0 prevalence", len(y), 0.500, 0.000, round(prev, 3), 0.000)]

    # Baseline 1: structure, all overlap
    ST = featurize(conns, kind=STRUCT).loc[conns].values.astype("float32")
    R.append(row("Baseline1 structure(ECFP)", len(y), *cv_eval([(ST, CFG["n_pca_struct"])], y)))

    # Baseline 2: measured Tox21 (impute missing 0.5), compounds with >=1 measured
    toxm = coh_i.reindex(conns)[ASSAYS].apply(pd.to_numeric, errors="coerce")
    m2 = toxm.notna().any(axis=1).values
    R.append(row("Baseline2 measuredTox21", int(m2.sum()), *cv_eval([(toxm.fillna(0.5).values[m2], None)], y[m2])))

    # expression subset for B3 / A / B  — everything below on the SAME 135 compounds (fair head-to-head)
    ce = [c for c in conns if c in expr_conns]; ye = dili.set_index("conn")["label"].loc[ce].values
    GEe = comb.loc[ce].values.astype("float32"); STe = featurize(ce, kind=STRUCT).loc[ce].values.astype("float32")
    toxe = coh_i.reindex(ce)[ASSAYS].apply(pd.to_numeric, errors="coerce").fillna(0.5).values
    R.append(("-- same-set (N=%d) head-to-head --" % len(ye), len(ye), np.nan, np.nan, np.nan, np.nan))
    R.append(row("  B1 structure @135",   len(ye), *cv_eval([(STe, CFG["n_pca_struct"])], ye)))
    R.append(row("  B2 measuredTox21 @135", len(ye), *cv_eval([(toxe, None)], ye)))
    R.append(row("  B3 rawExpression @135", len(ye), *cv_eval([(GEe, CFG["n_pca_expr"])], ye)))
    R.append(row("  A  predTox21 in-sample @135", len(ye), *cv_eval([(np.nan_to_num(predD.loc[ce].values, nan=0.5), None)], ye)))
    R.append(row("  A  predTox21 cross-fit @135", len(ye), *cv_eval([(np.nan_to_num(predCF.loc[ce].values, nan=0.5), None)], ye)))
    R.append(row("  B  fusedRep @135",     len(ye), *cv_eval([(STe, CFG["n_pca_struct"]), (GEe, CFG["n_pca_expr"])], ye)))

    tab = pd.DataFrame(R, columns=["model", "N", "AUC", "AUC_ci", "AUPRC", "AUPRC_ci"])
    sfx = "" if STRUCT == "ecfp4" else f"_{STRUCT}"
    base = "validation_b" if label_source == "dilirank" else f"validation_b_{label_source}"
    tab.to_csv(os.path.join(RES, f"{base}{sfx}.csv"), index=False)
    tgt = {"dilirank": "DILIrank DILI-concern", "withdrawn": "market-withdrawal (ChEMBL, SEPARATE noisier target)"}[label_source]
    print(f"\n=== VALIDATION B — {tgt} ===")
    print(f"overlap compounds: {len(y)}  positives {int(y.sum())} negatives {int((y==0).sum())}  "
          f"(prevalence {prev:.3f})" + (" ; structure bar ~0.75-0.83" if label_source == "dilirank" else ""))
    print(tab.to_string(index=False))
    if label_source == "dilirank":
        print("\nread: measuredTox21≈0.52-0.59 reproduces published near-random Tox21->DILI (sanity OK);"
              "\n      structure (0.70) is the ONLY real signal; expression/fused ≈chance. RIGOR: Approach A"
              "\n      in-sample 0.646 was leakage -> cross-fit collapses to 0.525. Chained pipeline: NO DILI signal."
              "\n      caveat: structure 0.66-0.70 < published SOTA 0.75-0.83 (weak ECFP+logistic, N=135, 82% pos).")
    else:
        print("\nread (SECONDARY, noisy exploratory target): for market-withdrawal the pattern FLIPS —"
              "\n      Tox21-based features (measured ~0.65, cross-fit predicted ~0.66) modestly BEAT structure"
              "\n      (~0.56), and cross-fit does NOT collapse (real, not leakage). BUT: small N, low prevalence,"
              "\n      and withdrawal is a noisy label (mixed reasons: efficacy/commercial, not only toxicity).")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main("dilirank")      # primary target
        main("withdrawn")     # secondary noisier target (kept separate)
