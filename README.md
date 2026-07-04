# Toxicogenomics fusion cohort — TG-GATEs + DrugMatrix × Tox21

Goal: build a training substrate for detecting **subtle / pre-histopathological
toxicity signatures** in animal transcriptome data, with **mechanism-of-toxicity
labels** attached. This directory contains the compound-level linkage layer:
every compound that has rat in-vivo perturbation data cross-referenced to Tox21
mechanism assays via standardized chemical identifiers.

## The deliverable: `master_cohort.csv`

One row per unique compound (keyed on InChIKey connectivity) across the two rat
transcriptome resources, annotated with data-source flags and Tox21 mechanism
labels.

| column | meaning |
|---|---|
| `compound_name` | display name (TG-GATEs name preferred, else DrugMatrix) |
| `inchikey` | full standardized InChIKey |
| `connectivity` | first 14 chars — the salt/stereo-insensitive join key |
| `in_tggates` / `tggates_liver` / `tggates_kidney` | present in Open TG-GATEs (rat in vivo) |
| `in_drugmatrix` | present in DrugMatrix (rat, multi-tissue) |
| `cross_platform_replicate` | in **both** transcriptome sources — use to separate biological signal from batch effects |
| `in_tox21` | has Tox21 assay data |
| `tox21_assays_labelled` / `tox21_assays_active` | of 12 assays, how many are measured / positive |
| `NR-AR … SR-p53` | the 12 Tox21 mechanism assay calls (0/1/blank) |

### Headline numbers

| | compounds |
|---|---|
| Unique compounds (TG-GATEs ∪ DrugMatrix) | 672 |
| — in TG-GATEs | 160 |
| — in DrugMatrix | 624 |
| — cross-platform replicates (both) | 112 |
| **With Tox21 mechanism labels** | **613** |
| TG-GATEs **liver** + Tox21 | 141 |

Tox21 coverage on shared compounds is dense (~10 of 12 assays measured per
compound), so this is not a missing-data swamp.

## Tox21 mechanism assays (the 12 labels)

Nuclear-receptor panel: `NR-AR`, `NR-AR-LBD`, `NR-AhR`, `NR-Aromatase`, `NR-ER`,
`NR-ER-LBD`, `NR-PPAR-gamma`. Stress-response panel: `SR-ARE` (Nrf2 oxidative
stress), `SR-ATAD5` (genotoxicity), `SR-HSE` (heat-shock/proteotoxic),
`SR-MMP` (mitochondrial membrane potential), `SR-p53` (DNA damage).

**Caveat:** Tox21 assays are *human cell-based in vitro*; the transcriptomes are
*rat in vivo*. Treat Tox21 as an orthogonal mechanism **prior**, not as ground
truth about the rat.

## Results — does gene expression add to structure?

Controlled comparison on the **177** cohort compounds that have DrugMatrix rat-liver expression
(single-source v1), predicting the 12 Tox21 endpoints as a multi-task problem with masked labels.
Every arm is identical except the feature block, so the contrast isolates *"does gene expression
add?"*. **Leakage-safe** protocol: repeated stratified 5-fold × 10 (50 splits, pooled out-of-fold);
every transform (PCA, scalers) fit on the training fold only. Full table:
[`data/results/results_table.csv`](data/results/results_table.csv).

| Arm | Features | Macro ROC-AUC | Macro AUPRC |
|---|---|---|---|
| Structure only (baseline) | ECFP4-2048 → PCA-128 | 0.757 | 0.538 |
| Expression only (ablation) | logFC → PCA-100 | 0.679 | 0.374 |
| **Fusion (structure + expression)** | **[struct-128, GE-100]** | **0.766** | **0.548** |

**Headline: fusion beats structure by ΔAUC +0.009 macro — but the gain is not uniform.** It
concentrates in the **stress-response (SR)** assays (and PPAR-γ / ER), and is neutral-to-negative on
the receptor-**binding** endpoints where structure already saturates:

| Assay | Structure AUC | Fusion AUC | ΔAUC | Fusion wins |
|---|---|---|---|---|
| SR-p53 | 0.700 | 0.749 | **+0.048** | 100% |
| NR-PPAR-γ | 0.764 | 0.812 | **+0.048** | 90% |
| SR-MMP | 0.778 | 0.819 | **+0.042** | 80% |
| NR-ER | 0.641 | 0.667 | +0.027 | 90% |
| SR-ARE | 0.635 | 0.652 | +0.017 | 90% |
| SR-HSE | 0.538 | 0.551 | +0.013 | 60% |
| SR-ATAD5 | 0.829 | 0.834 | +0.005 | 40% |
| NR-AhR | 0.819 | 0.813 | −0.006 | 30% |
| NR-AR | 0.795 | 0.785 | −0.010 | 50% |
| NR-AR-LBD | 0.899 | 0.881 | −0.017 | 10% |
| NR-Aromatase | 0.760 | 0.733 | −0.027 | 20% |
| NR-ER-LBD | 0.921 | 0.894 | −0.028 | 10% |

**The SR>NR pattern.** Rat-liver transcriptional response adds most where the toxic mechanism is a
*cellular stress program* — oxidative stress (SR-ARE), DNA damage (SR-p53), mitochondrial/heat-shock
(SR-MMP/HSE) — signals a static structure fingerprint can't see. It adds little on *ligand-binding*
endpoints (ER/AR-LBD, aromatase) that structure already determines well. That mechanism-split is the
main scientific takeaway.

**Caveats — this is a first-pass signal, not the final number:**
- **Single source, 177 / 613 compounds** — DrugMatrix liver only; TG-GATEs not yet fused (would raise coverage and enable the cross-source ComBat validation on the 112 shared compounds).
- **Structure baseline = ECFP4 + logistic regression**, not the planned frozen-ChemBERT + MLP head. ECFP4 is a strong Tox21 baseline, but a stronger structure arm could shrink the apparent expression gain, so the magnitude is provisional.
- No batch correction is applied (unnecessary for a single platform).

## Files

Data:
- `master_cohort.csv` — the fusion table (above)
- `data/tox21_keyed.csv` — 7,586 unique Tox21 compounds, standardized + keyed, with 12 assays
- `data/tggates_keyed.csv` — 170 TG-GATEs compounds keyed (160 small molecules; rest are biologics/mixtures with no structure)
- `data/dm_keyed.csv` — 641 DrugMatrix chemicals keyed (635 resolved)
- `data/shared_tggates_tox21.csv`, `data/dm_shared_tox21.csv` — per-source intersections with assay calls
- `data/open_tggates_main.csv` — source compound list (LSDB Archive)
- `data/chemicals.Rds` — source DrugMatrix treatment table (combspk/Complete-DrugMatrix)

Scripts (run from `data/`):
- `scripts/key_tox21.py` — standardize + InChIKey the Tox21 SMILES  → `tox21_keyed.csv`
- `scripts/resolve_tggates.py` — TG-GATEs names → InChIKeys via PubChem → `tggates_keyed.csv`
- `scripts/resolve_dm.py` — DrugMatrix names → InChIKeys via PubChem → `dm_keyed.csv`
- `scripts/intersect.py` — print all overlap statistics
- `scripts/build_master.py` — assemble `master_cohort.csv`

## Reproducing the data (large files not in this repo)

The raw + processed expression matrices are intentionally **not committed** (size + third-party
redistribution terms). The repo ships the code, the linkage cohort, the small derived tables and
the results table; everything else regenerates from **one public download**. Full instructions and
the list of datasets to fetch are in **[DATA.md](DATA.md)**.

```bash
pip install -r requirements.txt
mkdir -p data/_raw
curl -L -o data/_raw/GSE57815_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57815/matrix/GSE57815_series_matrix.txt.gz
python scripts/retrieve_expression.py && python scripts/build_signatures.py && python scripts/run_experiment.py
```

The modeling scripts (`retrieve_expression`, `build_signatures`, `expr_embed`, `structure_embed`,
`run_experiment`) sit alongside the curation scripts listed above; `run_experiment.py` produces
`data/results/results_table.csv` (structure vs expression vs fusion, leakage-safe).

## Method notes

- **Standardization**: RDKit `Cleanup` → `LargestFragmentChooser` (desalt) →
  `Uncharger` (neutralize) → InChI → InChIKey. Applied identically to all three
  sources so keys are comparable.
- **Matching**: on both the full InChIKey and the 14-char connectivity block.
  Connectivity matching recovered ~7% more real compounds (salt/stereo form
  mismatches, e.g. ketoconazole, rifampicin, etoposide).
- **Name resolution**: PubChem PUG-REST `name → SMILES`. Note PubChem renamed
  its output field to `SMILES`/`ConnectivitySMILES` (old `CanonicalSMILES`/
  `IsomericSMILES` are gone). Calls go through `curl` because the sandbox TLS
  proxy uses a cert Python's `ssl` module does not trust.

## Provenance / sources

- Open TG-GATEs — LSDB Archive: https://dbarchive.biosciencedbc.jp/en/open-tggates/
- DrugMatrix — NIEHS/NTP; compound table via https://github.com/combspk/Complete-DrugMatrix
- Tox21 — MoleculeNet (`tox21.csv`), EPA/NCATS Tox21 program
- PubChem PUG-REST — https://pubchem.ncbi.nlm.nih.gov/

## Next steps

1. Pull expression + histopathology for the **112 cross-platform compounds** to
   test the core hypothesis (transcriptome moves before histopathology fires)
   with built-in replication.
2. Add clinical failure labels (DILIrank, ClinTox, WITHDRAWN) keyed the same way.
3. Extend beyond liver using DrugMatrix's other tissues.
