#!/usr/bin/env python3
"""
build_signatures.py — collapse per-sample expression into ONE logFC vector per compound.

Collapse decision (all config-driven; see CONFIG). Defaults + rationale:
  signature = logFC (treated - matched vehicle), matched on VEHICLE + TIMEPOINT within
              the same platform  -> isolates the drug effect, cancels baseline (free 1st-pass
              batch correction).
  tissue    = liver              -> already all-liver in GSE57815.
  dose      = high (max mg/kg per compound) -> strongest, most reliable signal.
  time      = aggregate across timepoints (mean | maxabs) -> the two DBs use different time
              grids, so aggregating keeps signatures comparable across sources.
  replicate = mean               -> standard noise reduction.
  level     = probe              -> DrugMatrix-only first pass; switch to 'gene' (GPL1355
              collapse) when TG-GATEs joins and cross-source gene alignment is needed.

Also runs a SENSITIVITY check: default (aggregated) vs (high-dose / single latest timepoint),
reports per-compound correlation so the collapse choice is documented, not silent.

Output (data/signatures/):
  <src>_logfc.parquet     compounds (connectivity) x probes, float32   [the signature matrix]
  labels.csv              connectivity x 12 Tox21 assays (0/1/blank), aligned to covered cpds
  signature_report.txt    coverage + sensitivity summary
"""
import os, re
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPR = os.path.join(ROOT, "data", "expression")
SIG  = os.path.join(ROOT, "data", "signatures"); os.makedirs(SIG, exist_ok=True)

CONFIG = dict(
    source     = "drugmatrix_liver",
    dose       = "high",     # 'high' (max mg/kg) | 'all'
    time_agg   = "mean",     # 'mean' | 'maxabs'
    level      = "probe",    # 'probe' | 'gene' (gene needs GPL1355 map; added at TG-GATEs merge)
    control    = "vehicle_time",  # match control on vehicle+time, graded fallback
)
ASSAYS = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma",
          "SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]

def _days(s):
    m = re.match(r"([\d.]+)\s*d", str(s).strip().replace(" ", ""))
    return float(m.group(1)) if m else np.nan

def _mgkg(s):
    m = re.match(r"([\d.]+)", str(s).strip())
    return float(m.group(1)) if m else np.nan

def build(cfg=CONFIG):
    expr = pd.read_parquet(os.path.join(EXPR, f"{cfg['source']}_expr.parquet"))  # probes x samples
    man  = pd.read_csv(os.path.join(EXPR, f"{cfg['source']}_manifest.csv"))
    man["day"]  = man["time"].map(_days)
    man["mgkg"] = man["dose"].map(_mgkg)

    # ---- control baselines: mean log-intensity per (vehicle, day), with fallbacks ----
    ctl = man[man.is_control]
    def ctl_mean(sub): return expr[sub["sample"].tolist()].mean(axis=1) if len(sub) else None
    by_vt = {k: ctl_mean(g) for k, g in ctl.groupby(["vehicle", "day"])}
    by_v  = {k: ctl_mean(g) for k, g in ctl.groupby("vehicle")}
    glob  = ctl_mean(ctl)
    def baseline(veh, day):
        return by_vt.get((veh, day), by_v.get(veh, glob))

    # ---- per-sample logFC for treated samples ----
    trt = man[~man.is_control & man.in_labelled].copy()
    logfc = {}   # sample -> Series(probe)
    for _, r in trt.iterrows():
        logfc[r["sample"]] = expr[r["sample"]].values - baseline(r["vehicle"], r["day"]).values
    L = pd.DataFrame(logfc, index=expr.index)   # probes x treated-samples

    # ---- collapse to one vector per compound ----
    def collapse_compound(conn, variant):
        s = trt[trt.connectivity == conn]
        if cfg["dose"] == "high" and s["mgkg"].notna().any():
            s = s[s["mgkg"] == s["mgkg"].max()]
        if variant == "single_latest":                       # sensitivity variant
            s = s[s["day"] == s["day"].max()]
        sub = L[s["sample"].tolist()]
        # replicates+dose -> mean within day, then aggregate across days
        per_day = sub.T.groupby(s.set_index("sample")["day"]).mean().T   # probes x days
        if cfg["time_agg"] == "maxabs" and variant == "default":
            idx = per_day.abs().values.argmax(axis=1)
            return per_day.values[np.arange(len(per_day)), idx]
        return per_day.mean(axis=1).values

    conns = sorted(trt.connectivity.unique())
    default = np.vstack([collapse_compound(c, "default")       for c in conns])
    single  = np.vstack([collapse_compound(c, "single_latest") for c in conns])
    sig = pd.DataFrame(default, index=conns, columns=expr.index).astype("float32")
    sig.index.name = "connectivity"

    # ---- align labels ----
    coh = pd.read_csv(os.path.join(ROOT, "master_cohort.csv"), dtype=str).set_index("connectivity")
    labels = coh.loc[conns, ASSAYS]
    labels.to_csv(os.path.join(SIG, "labels.csv"))
    sig.to_parquet(os.path.join(SIG, f"{cfg['source']}_logfc.parquet"))

    # ---- sensitivity: default vs single-latest, per-compound correlation ----
    cors = [np.corrcoef(default[i], single[i])[0, 1] for i in range(len(conns))]
    rep = (f"signature source: {cfg['source']}  |  config: {cfg}\n"
           f"compounds: {len(conns)}  |  probes: {sig.shape[1]}\n"
           f"logFC range: {np.nanmin(default):.2f} .. {np.nanmax(default):.2f} "
           f"(median|logFC| {np.nanmedian(np.abs(default)):.3f})\n"
           f"sensitivity (aggregated-time vs single-latest-timepoint):\n"
           f"  per-compound corr  median={np.nanmedian(cors):.3f}  "
           f"min={np.nanmin(cors):.3f}  frac>0.8={np.mean(np.array(cors)>0.8):.2f}\n"
           f"label coverage: {labels.notna().sum().sum()} cells, "
           f"{(labels=='1').sum().sum()} actives across 12 assays\n")
    with open(os.path.join(SIG, "signature_report.txt"), "w") as f: f.write(rep)
    print(rep)
    print(f"wrote {SIG}/{cfg['source']}_logfc.parquet  ({sig.shape[0]} x {sig.shape[1]})")

if __name__ == "__main__":
    build()
