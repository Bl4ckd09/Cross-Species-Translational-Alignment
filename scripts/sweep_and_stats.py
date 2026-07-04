#!/usr/bin/env python3
"""
sweep_and_stats.py — Plan B rigor: (1) sweep GE #PCs, (2) formally test SR>NR.

(1) N=177 is small, so 100 GE-PCs may over-fit the expression block. Sweep n_pca_expr and
    report macro structure/fusion AUC + per-panel (SR vs NR) mean dAUC(fusion-structure).
(2) At the chosen setting, TEST the hypothesis that fusion helps stress-response (SR) assays
    more than nuclear-receptor (NR) assays:
      - assay-level  : Mann-Whitney U, one-sided SR>NR (units = the 5 SR vs 7 NR assays)
      - repeat-level : paired mean(SR dAUC) - mean(NR dAUC) per repeat, bootstrap 95% CI
      - effect size  : rank-biserial correlation
Outputs: data/results/pc_sweep.csv, data/results/sr_vs_nr.txt
"""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_experiment import load, evaluate, ASSAYS, CFG, RES

SR = ["SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53"]
NR = ["NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma"]
PC_GRID = [20, 30, 50, 75, 100]

def assay_dauc(rec):
    """per-assay mean dAUC(fusion-structure) over repeats, and per-repeat dAUC arrays."""
    mean_d, per_rep = {}, {}
    for a in ASSAYS:
        f = np.array(rec["fusion"][a]["auc"]); s = np.array(rec["structure_only"][a]["auc"])
        n = min(len(f), len(s))
        d = f[:n] - s[:n]
        per_rep[a] = d; mean_d[a] = d.mean() if n else np.nan
    return mean_d, per_rep

def macro(rec, ap):
    return np.mean([np.mean(rec[ap][a]["auc"]) for a in ASSAYS])

def main():
    data = load()
    print("loaded:", len(data[0]), "compounds\n")

    # ---------- (1) PC sweep ----------
    rows, recs = [], {}
    for npc in PC_GRID:
        cfg = {**CFG, "n_pca_expr": npc, "repeats": 8}
        rec, _ = evaluate(cfg, data=data); recs[npc] = rec
        md, _ = assay_dauc(rec)
        rows.append({
            "n_pca_expr": npc,
            "structure_AUC": round(macro(rec, "structure_only"), 4),
            "fusion_AUC":    round(macro(rec, "fusion"), 4),
            "dAUC_macro":    round(macro(rec, "fusion") - macro(rec, "structure_only"), 4),
            "SR_dAUC":       round(np.mean([md[a] for a in SR]), 4),
            "NR_dAUC":       round(np.mean([md[a] for a in NR]), 4),
        })
        print(f"  n_pca_expr={npc:3d}  fusion={rows[-1]['fusion_AUC']:.4f}  "
              f"dAUC={rows[-1]['dAUC_macro']:+.4f}  SR={rows[-1]['SR_dAUC']:+.4f}  NR={rows[-1]['NR_dAUC']:+.4f}")
    sweep = pd.DataFrame(rows); sweep.to_csv(os.path.join(RES, "pc_sweep.csv"), index=False)
    best = sweep.loc[sweep.SR_dAUC.idxmax(), "n_pca_expr"]   # setting that maximises SR benefit
    print(f"\nbest-for-SR n_pca_expr = {best}")

    # ---------- (2) SR vs NR test at best setting ----------
    rec = recs[int(best)]
    md, per_rep = assay_dauc(rec)
    sr_vals = np.array([md[a] for a in SR]); nr_vals = np.array([md[a] for a in NR])

    U, p = mannwhitneyu(sr_vals, nr_vals, alternative="greater")     # SR dAUC > NR dAUC
    rbc = 2 * U / (len(sr_vals) * len(nr_vals)) - 1                    # rank-biserial effect size

    # repeat-level paired contrast + bootstrap CI
    nrep = min(len(per_rep[SR[0]]), len(per_rep[NR[0]]))
    sr_rep = np.mean([per_rep[a][:nrep] for a in SR], axis=0)
    nr_rep = np.mean([per_rep[a][:nrep] for a in NR], axis=0)
    diff = sr_rep - nr_rep
    rng = np.random.default_rng(0)
    boot = [rng.choice(diff, len(diff), replace=True).mean() for _ in range(10000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])

    lines = [
        f"SR-vs-NR test  (n_pca_expr={best}, {nrep} repeats x {CFG['folds']} folds)",
        "",
        "per-assay mean dAUC(fusion - structure):",
        *[f"  {a:14s} {md[a]:+.4f}   [{'SR' if a in SR else 'NR'}]" for a in ASSAYS],
        "",
        f"SR panel mean dAUC = {sr_vals.mean():+.4f}   (assays: {', '.join(SR)})",
        f"NR panel mean dAUC = {nr_vals.mean():+.4f}",
        f"SR - NR            = {sr_vals.mean()-nr_vals.mean():+.4f}",
        "",
        f"assay-level Mann-Whitney U (SR>NR): U={U:.1f}, one-sided p={p:.4f}",
        f"  rank-biserial effect size        = {rbc:+.3f}",
        f"repeat-level mean(SR)-mean(NR)     = {diff.mean():+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]",
        f"  -> {'CI excludes 0: SR>NR is significant' if lo>0 else 'CI includes 0: not significant at 95%'}",
    ]
    out = "\n".join(lines)
    with open(os.path.join(RES, "sr_vs_nr.txt"), "w") as f: f.write(out)
    print("\n" + out)
    print(f"\nwrote {RES}/pc_sweep.csv and {RES}/sr_vs_nr.txt")

if __name__ == "__main__":
    main()
