#!/usr/bin/env python3
"""
build_tggates_signatures.py — collapse RMA'd TG-GATEs liver CELs into one logFC vector/compound.

Input : data/expression/tggates_liver_rma.tsv   (from scripts/rma_tggates.R, probes x samples)
        data/_raw/tggates_targets.csv            (sample -> compound/conn/dose/time/role)
Output: data/signatures/tggates_liver_logfc.parquet  (connectivity x 31099 probes, float32)

Collapse (mirrors DrugMatrix): logFC = treated - TIMEPOINT-matched control mean; dose = high
(already selected in the target list); time = aggregate across available timepoints (mean);
replicates = mean.
Note: TG-GATEs controls are matched on timepoint (vehicle-specific matching not in the SDRF) —
a documented approximation; ComBat + downstream PCA absorb residual baseline.
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPR = os.path.join(ROOT, "data", "expression")
SIG  = os.path.join(ROOT, "data", "signatures"); os.makedirs(SIG, exist_ok=True)

def build():
    expr = pd.read_csv(os.path.join(EXPR, "tggates_liver_rma.tsv"), sep="\t", index_col=0)
    expr.columns = [c.strip() for c in expr.columns]
    man = pd.read_csv(os.path.join(ROOT, "data", "_raw", "tggates_targets.csv"))
    man["sample"] = man["cel"].str.replace(r"\.CEL$", "", regex=True, case=False).str.strip()
    man = man[man["sample"].isin(expr.columns)].copy()
    print(f"RMA matrix: {expr.shape[0]} probes x {expr.shape[1]} samples | manifest matched: {len(man)}")

    # timepoint-matched control baseline
    ctl = man[man.role == "control"]
    base_by_t = {t: expr[g["sample"].tolist()].mean(axis=1) for t, g in ctl.groupby("time_n")}
    glob = expr[ctl["sample"].tolist()].mean(axis=1)
    def baseline(t): return base_by_t.get(t, glob)

    # per-sample logFC (treated), then collapse per compound
    trt = man[(man.role == "treated") & man.conn.astype(str).str.len().gt(0)].copy()
    sigs, conns = [], []
    for conn, g in trt.groupby("conn"):
        if not conn or conn == "nan": continue
        # per-timepoint mean logFC, then mean across timepoints
        per_t = []
        for t, gt in g.groupby("time_n"):
            lfc = expr[gt["sample"].tolist()].sub(baseline(t).values, axis=0).mean(axis=1)
            per_t.append(lfc.values)
        sigs.append(np.mean(per_t, axis=0)); conns.append(conn)

    sig = pd.DataFrame(np.vstack(sigs), index=conns, columns=expr.index).astype("float32")
    sig.index.name = "connectivity"
    sig.to_parquet(os.path.join(SIG, "tggates_liver_logfc.parquet"))
    print(f"wrote tggates_liver_logfc.parquet: {sig.shape[0]} compounds x {sig.shape[1]} probes")
    print(f"logFC range {np.nanmin(sig.values):.2f}..{np.nanmax(sig.values):.2f} "
          f"median|logFC| {np.nanmedian(np.abs(sig.values)):.3f}")

if __name__ == "__main__":
    build()
