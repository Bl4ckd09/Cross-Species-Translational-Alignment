# Cross-Species Translational Alignment — TG-GATEs + DrugMatrix × Tox21

Goal: build a training substrate for detecting **subtle / pre-histopathological
toxicity signatures** in animal transcriptome data, with **mechanism-of-toxicity
labels** attached. This directory contains the compound-level linkage layer:
every compound that has rat in-vivo perturbation data cross-referenced to Tox21
mechanism assays via standardized chemical identifiers.

## Background — the hackathon

Built at **Building an AI Scientist**, an AI-for-Drug-Discovery hackathon by **TernaryTx, future.bio,
Pluto House & Anthropic** (3–5 July 2026, 50Y Soho Square, London —
[agenda](https://drive.google.com/file/d/1ldM0mJ_SW8jphKHAufZkjvD9TsGCCQG9/view)). The brief: use
agentic AI to make real progress on a drug-discovery problem over one weekend, judged on innovation,
technical execution, scientific relevance, potential impact and presentation.

Our question: **does rat in-vivo gene expression add anything to chemical structure when predicting
Tox21 mechanism-of-toxicity outcomes?** We built the full pipeline, got an encouraging first signal,
then stress-tested it — and report where it held and where it didn't (see [Results](#results--does-gene-expression-add-to-structure)).

- Code: <https://github.com/Bl4ckd09/Cross-Species-Translational-Alignment>
- Data + results: <https://huggingface.co/datasets/Marcolini/cross-species-translational-alignment>

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

We ran the controlled comparison (identical everything except the feature block: structure only /
expression only / fusion), **leakage-safe** — repeated stratified CV, every data-dependent transform
(PCA, scalers, ComBat) fit on the training fold only — then stress-tested the headline three ways.
**The finding did not replicate.**

| Stage | Setup | Fusion vs structure | SR ΔAUC | SR-vs-NR |
|---|---|---|---|---|
| First run | N=177, DrugMatrix liver, **logistic** head | +0.009 macro (0.766 vs 0.757) | **+0.025** | p=0.074 |
| Baseline 2 | N=177, **GBM** head (chemprop stand-in) | **−0.010** macro | −0.016 | — |
| Plan A | **N=256** (+TG-GATEs, ComBat-merged) | 0.752 macro (**−0.014 vs N=177**) | **+0.001** | p=0.38 |

At N=177 with a linear head, fusion adds a small benefit concentrated in the **stress-response (SR)**
assays (SR-p53, SR-MMP, PPAR-γ) and neutral-to-negative on receptor-**binding** endpoints — an
appealing "expression sees stress programs that structure can't" story. But that signal **reverses
under a stronger nonlinear head** (Baseline 2, gradient-boosted trees) and **washes out when the
second data source is added** (Plan A: SR ΔAUC +0.025 → +0.001; SR-vs-NR p=0.074 → 0.38; fusion macro
0.766 → 0.752). The ComBat alignment itself worked — 79 TG-GATEs compounds merged cleanly onto the
DrugMatrix reference — the *effect* just isn't robust.

**Honest conclusion:** the apparent SR-specific gain was largely a **small-sample / single-source /
linear-head artifact**. More data *and* a stronger classifier each remove it. That only becomes
visible if you actually run the robustness checks — which is the point.

Per-stage numbers:
[`results_table.csv`](data/results/results_table.csv) (N=177 three-arm),
[`baseline2.csv`](data/results/baseline2.csv) (logistic vs GBM head),
[`sr_vs_nr.txt`](data/results/sr_vs_nr.txt) (SR-vs-NR test),
[`pc_sweep.csv`](data/results/pc_sweep.csv) (PCA-component sweep),
[`results_combined.csv`](data/results/results_combined.csv) (per-assay N=177 → N=256).

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

Modeling & experiment scripts:
- `scripts/retrieve_expression.py` — parse GSE57815 → DrugMatrix liver expression matrix + manifest
- `scripts/build_signatures.py` — per-compound logFC signatures (treated − vehicle), configurable collapse
- `scripts/structure_embed.py` — SMILES → ECFP4 fingerprints (swappable `featurize()`; ChemBERT drop-in)
- `scripts/expr_embed.py` — expression → PCA latent (swappable `embed()`; future rat→human drop-in)
- `scripts/run_experiment.py` — the controlled structure / expr / fusion comparison → `results_table.csv`
- `scripts/baseline2.py` — logistic vs GBM head robustness check → `baseline2.csv`
- `scripts/sweep_and_stats.py` — PCA-component sweep + formal SR-vs-NR test → `pc_sweep.csv`, `sr_vs_nr.txt`
- `scripts/rma_tggates.R` — RMA-normalise the TG-GATEs liver CELs (R + `affy`)
- `scripts/build_tggates_signatures.py` — TG-GATEs logFC signatures
- `scripts/combat_merge.py` — ComBat-align TG-GATEs onto DrugMatrix → `combined_logfc.parquet`
- `scripts/run_combined.py` — N=256 combined comparison, N=177 vs N=256 → `results_combined.csv`

## Reproducing

Large matrices (raw GEO, RMA output, feature parquets) are **not committed** — size + third-party
redistribution terms. Everything regenerates from the scripts plus public downloads; the full
dataset list is in **[DATA.md](DATA.md)**.

**N=177 (DrugMatrix only) — the first run:**
```bash
pip install -r requirements.txt
mkdir -p data/_raw
curl -L -o data/_raw/GSE57815_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57815/matrix/GSE57815_series_matrix.txt.gz
python scripts/retrieve_expression.py     # -> data/expression/drugmatrix_liver_expr.parquet
python scripts/build_signatures.py         # -> data/signatures/drugmatrix_liver_logfc.parquet + labels.csv
python scripts/run_experiment.py           # -> data/results/results_table.csv
python scripts/baseline2.py                # -> data/results/baseline2.csv       (logistic vs GBM head)
python scripts/sweep_and_stats.py          # -> data/results/pc_sweep.csv + sr_vs_nr.txt
```

**N=256 (add TG-GATEs) — the robustness expansion.** Same code, different signature file. Needs the
E-MTAB-799 liver CELs + R with `affy` (see [DATA.md](DATA.md)):
```bash
Rscript scripts/rma_tggates.R              # RMA -> data/expression/tggates_liver_rma.tsv
python scripts/build_tggates_signatures.py # -> data/signatures/tggates_liver_logfc.parquet
python scripts/combat_merge.py             # ComBat -> data/signatures/combined_logfc.parquet + combined_labels.csv
python scripts/run_combined.py             # -> data/results/results_combined.csv (N=177 vs N=256)
```

**Same pipeline for both N.** `run_experiment.py` is parameterised by which signature file its config
points at — `drugmatrix_liver_logfc.parquet` + `labels.csv` for N=177, `combined_logfc.parquet` +
`combined_labels.csv` for N=256. The two file sets are **separate**: running the N=256 pipeline does
**not** overwrite the N=177 inputs or `results_table.csv`, so you can re-run N=177 any time with just
`python scripts/run_experiment.py`.

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

The single-source signal didn't survive, so the next moves sharpen the test rather than declare
victory:
1. **Stronger structure arm** — swap ECFP4 → ChemBERT (a torch-env drop-in behind the same
   `featurize()` contract) to confirm the null holds against a stronger structure encoder.
2. **More data** — repeat-dose TG-GATEs and DrugMatrix's other tissues, beyond the liver single-dose set.
3. **Harder target** — move from Tox21 mechanism priors toward clinical failure labels (DILIrank,
   ClinTox, WITHDRAWN), where in-vivo transcriptomics may carry signal structure genuinely lacks.
