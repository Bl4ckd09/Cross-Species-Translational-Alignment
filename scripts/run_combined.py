#!/usr/bin/env python3
"""
run_combined.py — the Plan A payoff: re-run the controlled comparison on the ComBat-merged
N=256 set (DrugMatrix 177 + TG-GATEs-new 79) and test whether the SR benefit GREW vs N=177.

Prereq: scripts/combat_merge.py has written combined_logfc.parquet + combined_labels.csv.
Prints the per-assay table + SR-vs-NR test at N=256, and the N=177 -> N=256 deltas.
Writes data/results/results_combined.csv.
"""
import os
import numpy as np
import pandas as pd
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_experiment import evaluate, load, ASSAYS, CFG, RES

SR = ["SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53"]
NR = ["NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma"]

def summ(rec):
    macro = lambda ap: np.mean([np.mean(rec[ap][a]["auc"]) for a in ASSAYS])
    d = {a: np.mean(rec["fusion"][a]["auc"]) - np.mean(rec["structure_only"][a]["auc"]) for a in ASSAYS}
    return dict(struct=macro("structure_only"), expr=macro("expr_only"), fusion=macro("fusion"),
                SR=np.mean([d[a] for a in SR]), NR=np.mean([d[a] for a in NR]), per_assay=d)

def main():
    from scipy.stats import mannwhitneyu
    cfg256 = {**CFG, "signatures": "combined_logfc.parquet", "labels": "combined_labels.csv", "repeats": 10}
    print("=== N=256 (combined, post-ComBat) ===")
    rec256, (conns, Y) = evaluate(cfg256, data=load(cfg256))
    s256 = summ(rec256)

    print("--- N=177 (DrugMatrix only, reference) ---")
    rec177, _ = evaluate({**CFG, "repeats": 10}, data=load(CFG))
    s177 = summ(rec177)

    rows = []
    for a in ASSAYS:
        rows.append({"assay": a, "panel": "SR" if a in SR else "NR",
                     "dAUC_177": round(s177["per_assay"][a], 4),
                     "dAUC_256": round(s256["per_assay"][a], 4),
                     "grew": s256["per_assay"][a] > s177["per_assay"][a]})
    tab = pd.DataFrame(rows); tab.to_csv(os.path.join(RES, "results_combined.csv"), index=False)

    sr_v = np.array([s256["per_assay"][a] for a in SR]); nr_v = np.array([s256["per_assay"][a] for a in NR])
    U, p = mannwhitneyu(sr_v, nr_v, alternative="greater")

    print("\n=== per-assay GE benefit: N=177 -> N=256 ===")
    print(tab.to_string(index=False))
    print(f"\nMACRO fusion AUC : 177={s177['fusion']:.4f}  ->  256={s256['fusion']:.4f}")
    print(f"SR-panel dAUC    : 177={s177['SR']:+.4f}  ->  256={s256['SR']:+.4f}  "
          f"({'GREW' if s256['SR']>s177['SR'] else 'shrank'})")
    print(f"NR-panel dAUC    : 177={s177['NR']:+.4f}  ->  256={s256['NR']:+.4f}")
    print(f"SR-vs-NR @256    : Mann-Whitney U={U:.1f}, one-sided p={p:.4f}")
    print(f"\nwrote {RES}/results_combined.csv")

if __name__ == "__main__":
    main()
