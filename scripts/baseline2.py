#!/usr/bin/env python3
"""
baseline2.py — Baseline 2 (strong structure reference) + robustness check.

The plan's Baseline 2 = "structure-alone at its best" (a strong published SMILES->Tox21 model,
e.g. chemprop, retrained on our split). chemprop is a D-MPNN needing torch>=2.1 — deferred to
the teammates' torch env. Here we use a strong NON-torch stand-in: gradient-boosted trees
(HistGradientBoosting) on the ECFP structure features, same leakage-safe CV.

Two questions answered:
  (a) Is our logistic structure baseline weak?  -> compare structure_only under logistic vs GBM.
  (b) Does GE still add on top of a STRONGER nonlinear head? -> fusion vs structure under GBM.
Output: data/results/baseline2.csv
"""
import os
import numpy as np
import pandas as pd
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_experiment import load, evaluate, ASSAYS, CFG, RES

SR = ["SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53"]

def macro(rec, ap):        return np.mean([np.mean(rec[ap][a]["auc"]) for a in ASSAYS])
def panel(rec, ap, keys):  return np.mean([np.mean(rec[ap][a]["auc"]) for a in keys])
def sr_dauc(rec):          return np.mean([np.mean(rec["fusion"][a]["auc"]) -
                                           np.mean(rec["structure_only"][a]["auc"]) for a in SR])

def main():
    data = load()
    rows = []
    for model in ["logistic", "gbm"]:
        cfg = {**CFG, "model": model, "repeats": 6}
        rec, _ = evaluate(cfg, data=data)
        rows.append({
            "head": "logistic (Baseline1)" if model == "logistic" else "GBM (Baseline2)",
            "structure_AUC": round(macro(rec, "structure_only"), 4),
            "fusion_AUC":    round(macro(rec, "fusion"), 4),
            "macro_dAUC":    round(macro(rec, "fusion") - macro(rec, "structure_only"), 4),
            "SR_structure":  round(panel(rec, "structure_only", SR), 4),
            "SR_fusion":     round(panel(rec, "fusion", SR), 4),
            "SR_dAUC":       round(sr_dauc(rec), 4),
        })
        print(f"[{model:8s}] structure={rows[-1]['structure_AUC']:.4f}  fusion={rows[-1]['fusion_AUC']:.4f}  "
              f"macroΔ={rows[-1]['macro_dAUC']:+.4f}  SRΔ={rows[-1]['SR_dAUC']:+.4f}")
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(RES, "baseline2.csv"), index=False)
    print("\n" + tab.to_string(index=False))
    # interpretation
    log_s = tab.loc[0, "structure_AUC"]; gbm_s = tab.loc[1, "structure_AUC"]
    print(f"\n(a) structure-at-its-best: logistic {log_s:.3f} vs GBM {gbm_s:.3f} "
          f"-> baseline is {'near-ceiling (GE gain is not just a weak-baseline artifact)' if abs(gbm_s-log_s)<0.02 else 'model-sensitive'}")
    print(f"(b) GE still adds under GBM on SR panel: ΔAUC={tab.loc[1,'SR_dAUC']:+.4f} "
          f"({'YES — robust to head choice' if tab.loc[1,'SR_dAUC']>0 else 'no'})")
    print(f"\nwrote {RES}/baseline2.csv")

if __name__ == "__main__":
    main()
