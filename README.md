---
pretty_name: Cross-Species Translational Alignment (TG-GATEs + DrugMatrix × Tox21)
license: cc-by-4.0
tags:
  - toxicology
  - toxicogenomics
  - drug-discovery
  - cheminformatics
  - tox21
  - dili
  - transcriptomics
  - cross-species
task_categories:
  - tabular-classification
size_categories:
  - n<1K
configs:
  - config_name: master_cohort
    data_files: master_cohort.csv
  - config_name: results_tox21
    data_files: data/results/results_table.csv
  - config_name: results_dili
    data_files: data/results/validation_b.csv
---

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
- Visual overview: [`overview.html`](overview.html) — a self-contained field report (hypothesis → methods → every result → conclusions); open locally to view rendered.

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
(PCA, scalers, ComBat) fit on the training fold only — then stress-tested the headline four ways.
**The headline held only under a linear head at N=177; pooling a second dataset revealed the real
bottleneck is cross-dataset comparability, not sample size.**

| Stage | Setup | Fusion macro | SR ΔAUC | ComBat r | SR-vs-NR |
|---|---|---|---|---|---|
| First run | N=177, DrugMatrix liver, **logistic** head | 0.766 (vs 0.757) | **+0.025** | — | p=0.074 |
| Baseline 2 | N=177, **GBM** head (chemprop stand-in) | 0.723 (**−0.010**) | −0.016 | — | — |
| Plan A · single-dose | N=256, +TG-GATEs **hours** (mismatched), ComBat | 0.752 | +0.001 | 0.38 | p=0.38 |
| Plan A · repeat-dose | N=256, +TG-GATEs **days** (time-matched), ComBat | 0.756 | **+0.013** | **0.44** | p=0.27 |

At N=177 with a linear head, fusion adds a small benefit concentrated in the **stress-response (SR)**
assays (SR-p53, SR-MMP, PPAR-γ) and neutral-to-negative on receptor-**binding** endpoints — "expression
sees stress programs that structure can't." That signal **reverses under a stronger nonlinear head**
(Baseline 2, gradient-boosted trees, which overfits at N=177) and **washes out when a *mismatched*
second source is pooled** (single-dose, exposure = hours: SR +0.025 → +0.001, cross-dataset agreement
r=0.38). But the **fair, time-matched test** — TG-GATEs *repeat-dose* (exposure = days, bracketing
DrugMatrix's ≤7 d) — **recovers it partway** (SR → +0.013, agreement r=0.44). The recovered signal
*tracks* the agreement: 0.38 → 0.44 ⇒ +0.001 → +0.013.

**Honest conclusion:** the bottleneck is **cross-dataset comparability, not sample size.** More data
didn't help; better-*matched* data helped partially, in proportion to how comparable it was. Even
carefully harmonised, same-platform rat liver signatures for identical molecules top out at **r ≈ 0.44**
— a measured, in-rat preview of the animal→human translational gap. The clean **N=177** result stays
the primary evidence; the three-point pooling experiment (177 / hours-256 / days-256) is the supporting
story. **Full per-result detail and reasoning: [RESULTS.md](RESULTS.md).**

Per-stage numbers:
[`results_table.csv`](data/results/results_table.csv) (N=177 three-arm),
[`baseline2.csv`](data/results/baseline2.csv) (logistic vs GBM head),
[`sr_vs_nr.txt`](data/results/sr_vs_nr.txt) (SR-vs-NR test),
[`pc_sweep.csv`](data/results/pc_sweep.csv) (PCA-component sweep),
[`results_combined_singledose.csv`](data/results/results_combined_singledose.csv) +
[`results_combined_repeat.csv`](data/results/results_combined_repeat.csv) (per-assay N=177 → N=256,
hours vs days), pipeline logs in [`data/results/logs/`](data/results/logs/). Full narrative:
**[RESULTS.md](RESULTS.md)**.

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

Comparability — not data volume — is the bottleneck, so the next moves sharpen the test:
1. **Stronger structure arm** — *done*: swapped ECFP4 → ChemBERT (`STRUCT_KIND=chembert`). Frozen
   ChemBERT is actually *weaker* than ECFP for these tox tasks (Tox21 0.68 vs 0.76; DILI 0.59 vs
   0.70), so ECFP was the stronger baseline and the conclusions hold under both (see RESULTS.md §5b).
   A *fine-tuned* transformer remains the open follow-up.
2. **Better harmonisation** — go beyond ComBat (the r≈0.44 ceiling), e.g. tissue/time as covariates
   or a learned rat→rat alignment, before pooling. *(Done: repeat-dose day-matching, which lifted
   r 0.38 → 0.44.)*
3. **Harder target** — move from Tox21 mechanism priors toward clinical failure labels (DILIrank,
   DILIst, withdrawal lists), where in-vivo transcriptomics may carry signal structure genuinely
   lacks. *(This is the retrospective-validation plan now underway.)*
