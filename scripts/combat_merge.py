#!/usr/bin/env python3
"""
combat_merge.py — align TG-GATEs onto DrugMatrix (ComBat) and build the combined N=256 set.

- Both sources are Rat230-2 (31099 probesets) -> shared probe space, no ortholog step.
- ComBat with batch = data source, REF batch = DrugMatrix (so the N=177 DM result stays fixed
  and TG-GATEs is calibrated onto it).
- VALIDATE on the ~42 compounds present in both: correlation of DM vs TG signature should rise
  after correction. Report before/after.
- Final signature per compound: DrugMatrix for the 177 it has; ComBat-corrected TG-GATEs for the
  new compounds. -> data/signatures/combined_logfc.parquet + combined_labels.csv
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIG  = os.path.join(ROOT, "data", "signatures")
ASSAYS = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD","NR-PPAR-gamma",
          "SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]

def mean_shared_corr(dm, tg, shared):
    return float(np.mean([np.corrcoef(dm.loc[c], tg.loc[c])[0, 1] for c in shared]))

def main():
    dm = pd.read_parquet(os.path.join(SIG, "drugmatrix_liver_logfc.parquet"))
    tg = pd.read_parquet(os.path.join(SIG, "tggates_liver_logfc.parquet"))
    probes = dm.columns.intersection(tg.columns)
    dm, tg = dm[probes], tg[probes]
    print(f"DrugMatrix {dm.shape} | TG-GATEs {tg.shape} | shared probes {len(probes)}")

    shared = sorted(set(dm.index) & set(tg.index))
    new_tg = [c for c in tg.index if c not in set(dm.index)]
    print(f"overlap compounds (ComBat validation): {len(shared)} | new from TG-GATEs: {len(new_tg)}")

    # ---- ComBat (ref = DrugMatrix) ----
    from inmoose.pycombat import pycombat_norm
    X = pd.concat([dm, tg], axis=0)                       # (samples, probes), DM rows then TG rows
    batch = np.array([0] * len(dm) + [1] * len(tg))       # 0=DrugMatrix(ref), 1=TG-GATEs
    corr_before = mean_shared_corr(dm, tg, shared)
    Xc = pycombat_norm(X.T.values, batch, ref_batch=0).T  # correct in (probes,samples) space
    Xc = pd.DataFrame(Xc, index=X.index, columns=X.columns)
    tg_c = Xc.iloc[len(dm):]                              # corrected TG-GATEs rows
    corr_after = mean_shared_corr(dm, tg_c.loc[shared] if set(shared) <= set(tg_c.index) else tg_c, shared)
    print(f"\nComBat validation (mean DM-vs-TG signature corr on {len(shared)} shared compounds):")
    print(f"  before = {corr_before:+.4f}   after = {corr_after:+.4f}   "
          f"({'improved' if corr_after > corr_before else 'no improvement'})")

    # ---- combined signature: DM (unchanged) + corrected TG for NEW compounds ----
    combined = pd.concat([dm, tg_c.loc[new_tg]], axis=0).astype("float32")
    combined.index.name = "connectivity"
    combined.to_parquet(os.path.join(SIG, "combined_logfc.parquet"))

    coh = pd.read_csv(os.path.join(ROOT, "master_cohort.csv"), dtype=str).set_index("connectivity")
    labels = coh.loc[combined.index, ASSAYS]
    labels.to_csv(os.path.join(SIG, "combined_labels.csv"))
    src = pd.Series(["DrugMatrix"] * len(dm) + ["TG-GATEs"] * len(new_tg), index=combined.index, name="source")
    src.to_csv(os.path.join(SIG, "combined_source.csv"))
    print(f"\nwrote combined_logfc.parquet: {combined.shape[0]} compounds x {combined.shape[1]} probes")
    print(f"  DrugMatrix {len(dm)} + TG-GATEs-new {len(new_tg)} = {len(combined)}")
    print(f"  actives/assay:", {a: int((labels[a] == '1').sum()) for a in ASSAYS})

if __name__ == "__main__":
    main()
