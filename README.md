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
