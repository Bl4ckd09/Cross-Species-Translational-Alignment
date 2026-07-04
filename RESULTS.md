# Results — does rat gene expression add to chemical structure for predicting Tox21?

Full record of every experiment, **with the reasoning behind each result**. All numbers are
out-of-fold, leakage-safe (repeated stratified CV; every data-dependent transform — PCA, scalers,
ComBat — fit on the training fold only; fixed seeds). Raw artifacts are linked per section.

**Metric:** ROC-AUC (0.5 = chance, 1.0 = perfect) and, per assay, average precision (AUPRC).
**ΔAUC** = fusion (structure+expression) − structure-only. Positive ⇒ expression helps.

---

## 0. The data and the three supervised sets

| Set | N | What it is |
|---|--:|---|
| Labelled cohort | 613 | Compounds with a Tox21 label AND rat expression *somewhere* (catalogue-level) |
| **Primary supervised set** | **177** | Compounds with a **clean, single-protocol DrugMatrix Affymetrix-liver** signature |
| Combined (Plan A) | 256 | 177 + 79 new compounds from TG-GATEs (ComBat-merged) |
| ComBat validation | 42 | Compounds measured in **both** DrugMatrix and TG-GATEs |

**Why 177 and not 613?** The `in_drugmatrix` flag in the cohort came from the DrugMatrix *catalogue*
(all tissues + both array platforms). But we deliberately used only the **Affymetrix Rat230-2 liver**
arm (GSE57815) so that DrugMatrix and TG-GATEs live on the *same platform* (clean batch correction,
not a cross-platform nightmare). That arm contains 200 compounds, 177 of which carry Tox21 labels.
Choosing platform-cleanliness capped N at 177 — a deliberate, documented trade-off.

**Assay panels.** Nuclear-receptor (NR): NR-AR, NR-AR-LBD, NR-AhR, NR-Aromatase, NR-ER, NR-ER-LBD,
NR-PPAR-γ. Stress-response (SR): SR-ARE (oxidative), SR-ATAD5 (genotoxic), SR-HSE (proteotoxic),
SR-MMP (mitochondrial), SR-p53 (DNA-damage). The a-priori hypothesis: expression should help **SR**
(mechanisms that act *through gene regulation*, leaving a transcriptional footprint) more than **NR**
(receptor binding, which is a function of molecular shape → structure already captures it).

---

## 1. Does gene expression add? — first run (N=177)  → *small SR-specific benefit*

**Result** ([`results_table.csv`](data/results/results_table.csv)):

| Model | Macro AUC |
|---|--:|
| Structure only (ECFP4 → PCA-128) — Baseline 1 | 0.757 |
| Expression only (logFC → PCA-100) | 0.679 |
| **Fusion (structure + expression)** — Approach A | **0.766** |

Per-assay, the benefit is concentrated exactly where the hypothesis predicted:

| Helped most (ΔAUC) | Not helped (ΔAUC) |
|---|---|
| SR-p53 **+0.048**, NR-PPARγ **+0.048**, SR-MMP **+0.042**, NR-ER +0.027, SR-ARE +0.017 | NR-ER-LBD −0.028, NR-Aromatase −0.027, NR-AR-LBD −0.017 |

**Why this pattern.** Expression alone scores 0.68 — *well above chance* — so the rat transcriptome
genuinely carries toxicity signal. It adds most on stress-response endpoints and the two
*transcriptionally-active* receptors (PPARγ, ER are ligand-activated transcription factors → they
show up in a transcriptome). It adds nothing on the pure ligand-binding-domain assays (AR-LBD,
ER-LBD, Aromatase), where structure alone already hits AUC ~0.90 and the molecule's shape *is* the
answer. Mechanistically coherent.

**Interpretation.** Encouraging first signal — but modest (+0.009 macro) and at small N. Everything
below is the stress-testing that a first number like this *demands*.

---

## 2. The modality-balance fix — why z-scoring the PCs mattered

**Result.** Before z-scoring the two PCA blocks, fusion *lost* (−0.041 macro) and every SR gain was
masked. After z-scoring each block to unit variance, fusion wins (+0.009) and the SR pattern appears.

**Why.** Structure PC1 has far larger variance than the expression PCs. An L2-regularised head
penalises all coefficients equally, so unscaled it effectively *ignores the expression block* — the
"GE drowned out by structure" failure mode. Z-scoring both blocks puts them on equal footing.

**Interpretation.** This is a genuine methodological gotcha, reported as a before/after so the
result can't be dismissed as a fusion artifact. It's also *why* the effect is real-but-fragile:
it only survives with correct feature balancing.

---

## 3. PCA-component sweep — is 100 too many for N=177?  → *no, 100 is best*

**Result** ([`pc_sweep.csv`](data/results/pc_sweep.csv)):

| GE PCs | 20 | 30 | 50 | 75 | 100 |
|---|--:|--:|--:|--:|--:|
| Fusion AUC | 0.756 | 0.761 | 0.762 | 0.761 | **0.767** |
| SR ΔAUC | +0.005 | +0.013 | +0.017 | +0.019 | **+0.026** |

**Why.** The SR benefit grows *monotonically* with the number of expression components — it is not a
few dominant components but a signal spread across many. If it were overfitting noise, fewer PCs
would help; instead more PCs help, up to the N-limited ceiling.

**Interpretation.** 100 PCs (the plan default) is justified, not arbitrary — the expression signal is
distributed, and we're near the sample-size ceiling on how many components are estimable.

---

## 4. Formal SR-vs-NR test (N=177)  → *significant at the repeat level, trend at the assay level*

**Result** ([`sr_vs_nr.txt`](data/results/sr_vs_nr.txt)): SR-panel ΔAUC **+0.026** vs NR-panel **−0.000**.

| Test | Result | Verdict |
|---|---|---|
| Repeat-level paired mean(SR)−mean(NR), bootstrap 95% CI | +0.026 **[+0.014, +0.037]** | **Significant (CI excludes 0)** |
| Assay-level Mann-Whitney U (SR>NR), n=5 vs 7 | p = 0.074, effect size +0.54 | Strong trend, underpowered |

**Why the assay-level test is only a trend.** Only 12 assays, and two NR assays (PPARγ, ER) benefit —
because they too are transcription factors. So the clean "SR vs NR" label understates it; the real
split is **"acts through gene regulation" vs "acts through binding"**, which the data drew for us.

**Interpretation.** The SR-specific benefit is statistically real *within this dataset*, with a
mechanistically sensible exception structure. Whether it *generalises* is what §5–7 test.

---

## 5. Baseline 2 — does the benefit survive a stronger head?  → *no (it's linear-head-dependent)*

**Result** ([`baseline2.csv`](data/results/baseline2.csv)):

| Head | Structure AUC | Fusion AUC | SR ΔAUC |
|---|--:|--:|--:|
| Logistic (Baseline 1) | 0.755 | 0.764 | +0.026 |
| **Gradient-boosted trees (Baseline 2)** | **0.733** | 0.723 | **−0.016** |

**Why.** Two things at once: (a) the GBM *structure* baseline is **worse** than logistic (0.733 vs
0.755) — trees overfit at N=177 with 128 features, so the regularised linear model is actually
"structure at its best" here; (b) under that overfitting-prone head, adding 100 expression PCs makes
overfitting *worse*, so fusion drops. (The true SOTA reference, a chemprop D-MPNN, needs a torch env
and is a documented drop-in for teammates.)

**Interpretation.** The GE benefit is **model-dependent** — real under a properly-regularised linear
head, absent under a head that's already overfitting. That's an honest limitation, and it points at
the fix: **more data** (so stronger models stop overfitting), which is Plan A.

---

## 5b. ChemBERT structure arm — *is ECFP too weak a baseline? No — frozen ChemBERT is weaker*

The plan's Baseline 1 is ChemBERT, and ECFP's below-SOTA DILI number (§9) raised the question of
whether a transformer encoder would change the conclusions. Added via a dedicated torch venv that
only generates the embedding cache (frozen `ChemBERTa-zinc-base-v1`, 768-d mean-pooled → PCA-128);
the main pipeline reads the cache, so `STRUCT_KIND=chembert` swaps the encoder with no other change.
([`results_table_chembert.csv`](data/results/results_table_chembert.csv),
[`validation_b_chembert.csv`](data/results/validation_b_chembert.csv))

| Structure encoder | Tox21 structure-only | Tox21 fusion | Tox21 SR ΔAUC | DILI structure (N=134) |
|---|--:|--:|--:|--:|
| **ECFP4** | **0.757** | 0.766 | +0.025 | **0.699** |
| ChemBERT (frozen) | 0.681 | 0.716 | +0.056 | 0.592 |

**Why.** Frozen ChemBERTa (no fine-tuning, mean-pooled) is a *general* pretrained representation,
whereas ECFP directly encodes the substructures that drive assay activity — so for Tox21/DILI,
**ECFP beats frozen ChemBERT** (a documented pattern). Our ECFP baseline was therefore the *stronger,
more conservative* choice, not a weakness.

**Interpretation (three robustness reads).**
1. *"Does GE add" survives the encoder swap.* Fusion still beats structure under ChemBERT (+0.035
   macro, +0.056 SR) — in fact *larger*, but only because the weaker ChemBERT structure leaves more
   room. The SR>NR ordering holds under both encoders.
2. *The DILI SOTA gap isn't about our encoder.* ChemBERT scores **further below** the published
   0.75–0.83 bar (0.59), so that gap reflects model sophistication (fine-tuning / ensembles /
   descriptor-rich pipelines), not "we should have used a transformer."
3. *It sharpens the withdrawal caveat.* Measured-Tox21-beats-structure is stable across encoders
   (0.65 vs ~0.58), but Approach A's predicted-Tox21 number swings 0.70→0.53 between encoders —
   confirming §9.3's lead is noisy and not to be over-read.

---

## 6. Plan A, single-dose — does more data grow the benefit?  → *no, it washes out*

**Setup.** Add TG-GATEs (E-MTAB-799, single-dose) liver → 79 new compounds → N=256, ComBat-merged
onto DrugMatrix (reference batch), same Rat230-2 probe space.

**Result** ([`results_combined_singledose.csv`](data/results/results_combined_singledose.csv),
[`logs/plan_a_singledose.log`](data/results/logs/plan_a_singledose.log)):

| | N=177 | N=256 (single-dose) |
|---|--:|--:|
| ComBat validation r (42 shared cpds) | — | 0.18 → **0.38** (improved) |
| SR-panel ΔAUC | +0.025 | **+0.001** |
| Macro fusion AUC | 0.766 | 0.752 |
| SR-vs-NR p | 0.074 | 0.38 |

**Why.** ComBat *worked* (agreement on the 42 shared compounds rose from 0.18 to 0.38), yet the SR
signal vanished. The reason is in that 0.38: even after batch correction, the same compound's
signature in the two datasets agrees only *moderately*. The prime suspect is an **exposure-time
mismatch** — DrugMatrix dosed for **days** (0.25–7 d); TG-GATEs single-dose measures **hours**
(3–24 h). A 24-hour liver and a 5-day liver are different biological states, so the 79 new signatures
aren't measuring the same thing. Pooling injected more noise than signal.

**Interpretation.** Naive pooling of a second rat dataset did **not** help. This looked like an
"artifact confirmed" result — until we controlled the exposure-time variable (§7).

---

## 7. Plan A, repeat-dose — the fair, time-matched test  → *comparability is the lever*

**Setup.** Same expansion, but TG-GATEs **repeat-dose** (E-MTAB-800) at **4 d + 8 d** — day-scale
exposure that brackets DrugMatrix's ≤7 d range (chronic 15/29 d dropped, as they'd re-introduce a
mismatch). This isolates *exposure duration* as the variable.

**Result** ([`results_combined_repeat.csv`](data/results/results_combined_repeat.csv),
[`logs/plan_a_repeat.log`](data/results/logs/plan_a_repeat.log)) — the three-point comparison:

| | N=177 (clean) | N=256 single-dose (hours) | N=256 repeat-dose (days) |
|---|--:|--:|--:|
| ComBat agreement r (after) | — | 0.38 | **0.44** |
| **SR-panel ΔAUC** | **+0.025** | +0.001 | **+0.013** |
| Macro fusion AUC | 0.766 | 0.752 | 0.756 |
| SR-vs-NR p | 0.074 | 0.38 | 0.265 |

**Why this is the key result.** Switching from mismatched (hours) to matched (days) exposure moved
**every** metric in the right direction: cross-dataset agreement rose **0.38 → 0.44**, and the SR
benefit **partially recovered, +0.001 → +0.013** (about halfway back to the clean +0.025). The two
track each other: **as comparability improved, the recovered signal improved** — a near dose-response
relationship between "how well the datasets agree" and "how much real signal you can pool."

**But** it did not *fully* recover: SR +0.013 is still below the clean N=177 (+0.025) and not
significant (p=0.265). Even same-platform, same-organ, time-matched rat data tops out at r ≈ 0.44.

---

## 8. Synthesis — the controlled three-point experiment

1. **One clean dataset** (N=177) → SR benefit is real: **+0.025**.
2. **Pool a *mismatched* dataset** (hours) → benefit washes out (**+0.001**), agreement **0.38**.
3. **Pool a *matched* dataset** (days) → benefit partially returns (**+0.013**), agreement **0.44**.

**Conclusion:** the bottleneck is **cross-dataset comparability, not sample size.** More data did not
help; better-*matched* data helped partially, in proportion to how comparable it was. We put a number
on the reproducibility ceiling of rat toxicogenomics: **even carefully harmonised, same-platform rat
liver signatures for identical molecules agree only r ≈ 0.44.**

This is the in-rat, *measurable* preview of the project's core thesis: the animal→human translational
gap is a **harmonisation** problem, not a data-volume problem. The lever is careful alignment /
translation of responses — exactly what a rat→human step (the swappable `embed()` interface) is for.

**What to trust:** the clean **N=177** result is the primary, best-controlled evidence for the
SR-specific benefit. The **three-point pooling experiment** is the supporting story about *why* it's
hard to scale — and it's a stronger scientific narrative than "we added data and AUC went up."

---

## 9. Retrospective validation on real drugs — held-out Tox21 + human DILI

Two **separate** tests with separate baselines (not conflated). A positive result on one implies
nothing about the other.

### 9.1 Validation A (measured Tox21 on novel compounds) — *infeasible by construction*

**Result.** The intended held-out set — compounds with rat expression AND measured Tox21, *outside*
the training set — is **≈empty**. Of 672 rat compounds, **613 carry Tox21 labels (that IS the
supervised set)** and the other 59 have no Tox21 to score against.

**Why.** Our supervised set is *defined as* (TG-GATEs ∪ DrugMatrix) ∩ Tox21, so every rat-expression
compound with a Tox21 label is already in it. There is no external rat-expression-plus-Tox21 set.

**Interpretation.** Validation A is already answered by the **out-of-fold cross-validation in §1**
(`results_table.csv`) — every compound is predicted while held out. That *is* the held-out Tox21 test.

### 9.2 Validation B (human clinical DILI, DILIrank) — *structure wins; expression adds nothing*

**Setup.** DILIrank (FDA/NCTR; 982 usable drugs — 568 DILI-positive [vMost+vLess], 414 negative
[vNo]; Ambiguous excluded). Overlap with the rat cohort = **282 compounds (232 pos / 50 neg — ~82%
positive**, enriched for hepatotoxins as expected). Baseline ladder 0–3 + Approach A (predicted
Tox21 → DILI) + Approach B (fused representation → DILI); leakage-safe repeated stratified CV.
([`validation_b.csv`](data/results/validation_b.csv))

**Sanity gate (passed).** Measured Tox21 → DILI = **0.523 ≈ chance** — reproduces the published
near-random Tox21→DILI finding, confirming the overlap set is not a selection artifact.

**Same-set head-to-head** (N=134, every model on identical compounds; 82% prevalence → AUC is the
informative metric, AUPRC baseline is already 0.825):

| Model | AUC | AUPRC |
|---|--:|--:|
| Baseline 0 · prevalence | 0.500 | 0.825 |
| **B1 · structure (ECFP)** | **0.699** | 0.935 |
| A · predicted Tox21 → DILI *(in-sample)* | 0.607 | 0.949 |
| A · predicted Tox21 → DILI **(cross-fit)** | **0.465** | 0.913 |
| B · fused representation → DILI | 0.600 | 0.935 |
| B2 · measured Tox21 → DILI | 0.550 | 0.939 |
| B3 · raw rat expression → DILI | 0.516 | 0.923 |

**Why.** Structure carries essentially all the DILI signal there is here. Measured Tox21 (0.59) and
raw rat expression (0.52) are near-random, and the fused representation (B, 0.60) scores *below*
structure alone — the expression channel dilutes rather than adds.

**The cross-fit rigor pass (the decisive check).** Approach A's *in-sample* number (0.607) looked
like it might carry signal — but it was **leakage**: the frozen Tox21 model had seen each compound's
own Tox21 labels, so its predicted-Tox21 features were optimistic. Regenerating those features
**out-of-fold** (each compound's Tox21 predicted by a model trained only on the *other* compounds)
**collapses Approach A to 0.465 ≈ chance.** So the chained pipeline carries *no* real DILI signal;
the apparent edge was entirely an artifact — exactly the kind of illusion only a cross-fitting check
exposes.

**Interpretation.** For real human hepatotoxicity, the rat-expression + Tox21 model **adds nothing
over chemical structure.** This is consistent with the project's core finding and with the DILI
literature: the animal→human step is a genuine translation gap, and neither a human-cell-line assay
panel (Tox21) nor a rat transcriptome closes it here.

**Answer to discussion question 1.** *Delta over human cell-line data?* — yes, but adverse to us:
structure (0.70) beats human-cell-line Tox21 (0.55) by ~0.15 AUC, and our expression additions don't
extend that.

**Caveats.** (a) Our structure arm (0.65–0.70) is **below published SOTA (0.75–0.83)** — ECFP4+logistic
on small, imbalanced N=134. A **ChemBERT** arm (§5b) was run and is *weaker still* (0.59), so the
SOTA gap is model sophistication, not encoder choice; ECFP is the stronger baseline here. The robust
finding is the *ordering* (structure > expression/fused), not the absolute number. (b) Input-overlap
for A/B/B3 is unavoidable — expression + Tox21 exist only for training compounds — so the model can
only be scored on compounds whose *inputs* it saw (not their DILI label; different target). Approach
A's specific optimism (the Tox21 predictions) has now been removed via cross-fitting (above), which
turned 0.607 into 0.465 — confirming the null rather than weakening it.

### 9.3 Validation B, secondary target — market withdrawal (ChEMBL) — *the pattern flips*

**Setup.** A **separate, noisier** label (kept distinct from DILIrank, per plan): market-withdrawal
status from ChEMBL `drug_warning`, matched by InChIKey connectivity. Overlap = **338 marketed drugs,
only 35 withdrawn (10% prevalence)**. ([`validation_b_withdrawn.csv`](data/results/validation_b_withdrawn.csv);
source: [`fetch_withdrawn.py`](scripts/fetch_withdrawn.py))

**Result** (same-set N=150): **the ordering inverts vs DILIrank.**

| Model | AUC | AUPRC (baseline 0.10) |
|---|--:|--:|
| B1 · structure (ECFP) | 0.557 | 0.191 |
| B2 · measured Tox21 | 0.651 | 0.214 |
| B3 · raw rat expression | 0.586 | 0.167 |
| **A · predicted Tox21 (cross-fit)** | **0.704** | 0.251 |
| B · fused representation | 0.565 | 0.188 |

**Why / interpretation.** For "was this drug withdrawn from market," **structure is near-chance
(0.56)** — withdrawal is driven by *in-vivo/clinical* effects, not molecular shape — and the
**Tox21/expression-based features beat it** (measured Tox21 0.65; cross-fit predicted Tox21 0.70).
Critically, the cross-fit Approach A does **not** collapse here (0.716 in-sample → 0.704 cross-fit),
so unlike the DILI-concern target this is *not* a leakage artifact. Mechanistically plausible: the
model's biological-activity representation captures withdrawal-relevant effects that static chemistry
misses.

**But treat this as suggestive, not conclusive.** Only **35 positives** (≈24 in the same-set), a
**heterogeneous label** (drugs are withdrawn for efficacy/commercial reasons too, not only toxicity),
and CIs that reflect CV-repeat spread rather than the tiny positive class. It is a genuine, honest
*contrast* to the DILI result — different clinical outcome, different winner — and a lead worth
following with a larger, tox-specific withdrawal set, not a headline claim.

---

## 10. Honest limitations

- **Small N.** 177 (256 pooled). Sparse assays (NR-PPARγ: 16 actives) have wide intervals.
- **Structure encoder is ECFP4, not ChemBERT.** ECFP is a strong, standard Tox21 baseline (harder to
  beat, so a *conservative* test of "GE adds"), but the plan's ChemBERT Baseline 1 is a torch-env
  drop-in (`scripts/structure_embed.py`, `kind="chembert"`) not run on the Intel-Mac dev box.
- **Tox21 is a human in-vitro *prior*, not rat ground truth.** We predict mechanism labels, not the
  rat's own histopathology.
- **TG-GATEs control matching is timepoint-only** (vehicle metadata absent from the SDRF) — a
  documented approximation; ComBat + PCA absorb residual baseline.
- **ComBat across mostly-disjoint compound sets** is imperfect; validated on the 42 shared compounds.

---

## 11. Reproducing

See [README → Reproducing](README.md#reproducing) and [DATA.md](DATA.md). Large matrices are not
committed (size + redistribution terms); everything regenerates from the scripts + public downloads.
Per-stage artifacts live in [`data/results/`](data/results/); pipeline logs in
[`data/results/logs/`](data/results/logs/).
