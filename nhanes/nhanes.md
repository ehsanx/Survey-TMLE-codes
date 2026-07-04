# NHANES Application — Re-Analysis Plan (v1, 2026-06-07)

Plan for the real-data illustration of the cross-fitted survey-weighted TMLE paper (Biometric
Methodology). Supersedes the narrow "Step 4" note in `plan-phase3.md`. Goal: a **fully reproducible**,
**non-toy**, **full-estimator-suite** NHANES analysis, with everything (including raw data) under
`nhanes/`. Grounded in a 4-agent research pass (repo audit + recent-literature scan + `nhanesA` pipeline
review + covariate/DAG design); key external sources are listed in §9.

---

## 0. What changes vs. the current analysis

The current pipeline (`codes/nhanes1_save_analytic_data.R` → `nhanes2_run_tmle.R` → `nhanes3_run_glm.R`
→ `nhanes4_descriptives.R`) works but has five problems for a flagship methods paper:

1. **Only 3 single-fit arms.** `nhanes2` calls the `tmle` package directly and builds Non/Partially/
   Fully-Aware by hand; it does **not** call `R/estimators.R::run_estimators()`, so there is **no
   Fully-Aware-CF (cross-fitted) and no Fully-Aware-CV** arm — the two that the theory and simulation are
   about.
2. **Not offline-reproducible.** `nhanes1` re-downloads from CDC on every run and saves only the imputed
   analytic frame; raw NHANES is never cached. `nhanes/` is empty.
3. **A real data bug.** `nhanes1` requests `BPX_J` for 2017–2018, but CDC **renamed the BP table to
   `BPXO_J`** (oscillometric, variables `BPXOSY1/BPXODI1`) that cycle — so **blood pressure is silently
   all-`NA` for the entire 2017–2018 wave** and those people drop out. Must be handled per-cycle.
4. **Toy covariate set** (6 variables: age, race, education, PIR, BMI, macrosomia).
5. **Brittle/incorrect plumbing.** `nhanes3`/`nhanes4` read a **hard-coded path in a different repo**
   (`<repo-root>/...`); the subpopulation is created by **deleting rows before `svydesign`**
   (which biases SEs); and decoding relies on cross-cycle label-string matching.

This plan fixes all five.

---

## 1. Report the FULL estimator suite — not only Fully-Aware-CF

**Decision (answers the author's question):** the application reports **all five arms side by side**, with
**Fully-Aware-CF as the primary** estimate. Reporting only CF would discard the contrast that
makes the illustration persuasive.

| Arm | Weights | Variance | Nuisances | Role in the application |
|---|---|---|---|---|
| **Non-Aware** | none | i.i.d. | single-fit SL | shows weights matter (point estimate moves) |
| **Partially-Aware** | yes | i.i.d. | single-fit SL | shows design variance matters (CI too narrow) |
| **Fully-Aware** (single-fit) | yes | design (Eq 8) | single-fit SL | the Donsker-boundary arm (CI may be falsely narrow with rich SL) |
| **Fully-Aware-CV** | yes | design (Eq 8) | single-fit SL + cluster-aware CV folds | shows internal CV ≠ cross-fitting |
| **Fully-Aware-CF** ✅ | yes | design (Eq 8) | **PSU-level cross-fitted** SL | **primary** |
| *GLM benchmark* | 3 levels | design / i.i.d. | survey GLM (RD) | sanity check (`nhanes3` style) |

The applied story is the **echo of Figure 1**: the point estimate is stable while the **CI widens** as the
design is progressively respected (Non → Partially → Fully → CF), and single-fit SL / internal-CV can give
a misleadingly narrow CI that CF corrects. We also report **overlap/positivity diagnostics** (propensity
histograms by arm, realized min/max `g`, fraction near the truncation bound, weight-trimming sensitivity).

**Implementation note.** `run_estimators()` in `R/estimators.R` currently hard-codes the covariate
matrix as `W = c("L1","L2","L3","L4")` (from the sim). It must be **generalized to accept an arbitrary
covariate set** (a `W` formula/character vector argument) before it can run on the NHANES covariates.
Everything else (the `tmle()` targeting, the reconstructed-EIF + `svymean` design variance, `make_cf_folds`)
is reusable as-is.

**SL library for the application.** Report the five arms at one rich SL library (≈ ladder rung L3:
`glm + gam + earth + ranger`). Optionally add a **Web-Appendix-F "applied echo of Figure 1"**: single-fit
Fully-Aware vs Fully-Aware-CF across the same L1→L4 ladder on the NHANES data, to show the single-fit CI
narrowing as the library crosses the Donsker boundary.

---

## 2. Choice of example — is GDM → hypertension the right illustration?

**Assessment of the current example (GDM history → later-life hypertension, ever-pregnant women).**
Defensible but reads as a **toy / under-sells the method** for these reasons: (a) **rarity/positivity** —
self-reported GDM is ~6–9% of ever-pregnant women, so the treated cell is small and `g(W)` is pushed toward
0, undermining the very overlap-stability story the paper wants to showcase; (b) **temporality** — NHANES is
cross-sectional and the "confounders" (current BMI, PIR) are **contemporaneous/post-exposure**, so the
adjustment set is not temporally coherent; (c) **exposure misclassification** — decades-old self-recall;
(d) **noisy outcome** — a single `BPXSY1/BPXDI1` reading at ≥140/90 misclassifies treated-but-controlled
hypertensives as "No"; (e) **thin subpopulation** strains the clustered (Eq 8) variance; (f) **low
novelty** — a well-trodden association.

**Recommendation: switch the lead example to a high-overlap, contemporary question, and (optionally) keep a
*hardened* GDM → HTN as a deliberate "limited-overlap" secondary.** A large, well-overlapped example lets
all five arms be compared cleanly and gives a real design effect to display — exactly what the application
needs. Ranked candidates from the recent (2022–2025) NHANES literature (all binary→binary, ATE on the
risk-difference scale, NHANES 2007–2018; sources in §9):

1. **★ Short sleep (<7 h vs 7–9 h) → obesity (measured BMI ≥30).** *Top methods pick.* Both exposure
   (~30–40% short sleepers) and outcome (~37% obese) common → textbook overlap; large sample exercises the
   Super Learner; objective MEC outcome; rich confounders (diet energy, activity, SES, comorbidity); prior
   NHANES sleep–obesity papers **explicitly lack survey weighting** (the gap we fill). Caveat: cross-sectional
   temporality / reverse causation (state plainly).
2. **E-cigarette use (ever-vaped vs never, ~19%) → hypertension (≥130/80 OR self-report Dx OR
   antihypertensive use).** *Most topical (2025 policy interest).* Use ever-vs-never for overlap (current-vs-
   never ~5% strains positivity like GDM, keep as sensitivity). Wrinkle: confounding by combustible smoking.
3. **Food insecurity (USDA HFSSM ≥3 affirmatives) → depression (PHQ-9 ≥10).** *Best "social-determinants /
   absolute-scale" framing*; validated instruments both sides, ~14–19% exposed, large; the risk-difference is
   the policy-relevant quantity for a disparities exposure.
4. **Physical activity (≥150 min/wk MVPA) → depression (PHQ-9 ≥10).** Optional **paired** example to show a
   **protective / negative** RD (nice contrast). Excellent overlap.
5. *Sensitivity-only:* ultra-processed food → obesity (topical, but continuous exposure needs a justified
   cut and the **`WTDRD1`** dietary-day weight — a good teaching point) and PFAS → hypertension/gout
   (maximal topicality, forces the **`WTSA2YR` subsample weight**, but dichotomized biomarker + smaller
   subsample + reverse-causation → best as a vignette).

**Decision (author, 2026-06-07): analyze ALL FOUR examples, then choose the lead.** Rather than commit up
front, the pipeline runs the **full 5-arm suite + overlap diagnostics on all four** examples below, and the
author selects the lead (and an optional secondary) after comparing results. The rankings above remain
decision-support: short sleep→obesity and food insecurity→depression should give the cleanest overlap;
e-cig→HTN is the most topical; hardened GDM→HTN is the deliberate limited-overlap case. The
cross-example comparison is itself informative (how overlap, the Non→CF CI-widening, and single-fit-vs-CF
divergence vary across applications).

### 2.1 The four example configurations

All four are NHANES 2007–2018, **MEC-based so all use `WTMEC2YR/6`**, binary exposure → binary outcome,
estimand = population ATE (risk difference). The pipeline is parameterized by an **example registry**; each
entry supplies the extra components, exposure/outcome derivations, subpopulation, and confounder set.

| # | Exposure | Outcome | Subpopulation | Extra components | Overlap (≈exposed) |
|---|---|---|---|---|---|
| **E1** | Short sleep <7h vs 7–9h (`SLD010H`/`SLD012`) | Obesity, `BMXBMI`≥30 | adults ≥20 | SLQ, BMX | excellent (~35%) |
| **E2** | Food insecurity, HFSSM ≥3 (`FSDAD`/`FSD032*`) | Depression, PHQ-9 ≥10 (`DPQ010–090`) | adults ≥20 | FSQ, DPQ | good (~15%) |
| **E3** | Ever-vaped vs never (`SMQ900`/ECQ; current as sens.) | HTN: SBP≥130 or DBP≥80 OR `BPQ020` OR `BPQ040A` | adults ≥20 | SMQ/ECQ, BPX/BPXO, BPQ | good ever (~19%); strained current (~5%) |
| **E4** | GDM history `RHQ162` (hardened) | HTN: ≥130/80 OR Dx OR meds | ever-pregnant (`RHQ131`) non-diabetic women ≥20 | RHQ, DIQ, BPX/BPXO, BPQ | strained (~6–9%) |

*Exact variable names/availability per cycle must be confirmed during `00_download` via
`nhanesSearchVarName()` (the sleep-hours item wording changes at 2015–16; e-cig items are most complete in
2015–18; `BPX→BPXO_J` and `BPQ` per §3.2). Every example carries the **shared core confounders** (age, sex
where applicable, race/ethnicity, education, PIR, survey cycle) plus the **example-specific** confounders in
§4. **The lead/secondary choice is deferred to after the all-four run (see §8).***

---

## 3. Reproducible `nhanesA` pipeline

### 3.1 Folder layout (everything under `nhanes/`)

```
nhanes/
  nhanes.md                  ← this plan
  raw/                       ← immutable raw NHANES snapshot, one .rds per table (e.g. DEMO_E.rds)
  raw/codebooks/             ← cached nhanesTranslate() code→label maps per table
  raw/manifest.csv           ← table, nrow, ncol, download timestamp, nhanesA version, digest hash
  analytic/<example>_analytic.rds     ← per-example merged+harmonized+derived frame (E1..E4)
  analytic/<example>_imputed.rds      ← per-example single-imputation completed frame (+ design columns)
  results/<example>/estimators_all_arms.rds   ← run_estimators() (5 arms) + GLM benchmark, per example
  results/<example>/positivity_diagnostics.rds ← g(W) by arm, min/max g, weight-trim sensitivity
  results/<example>/table1_descriptives.rds    ← survey-weighted Table 1 (by exposure)
  results/cross_example_summary.csv   ← all four examples side by side (estimate/SE/CI/df by arm) for the choice
  figures/<example>/                  ← per-example forest plot of arms + propensity overlap plot
  R/                                  ← the rewritten, example-parameterized scripts (00_download … 06_compare)
  R/examples.R                        ← the EXAMPLE REGISTRY (E1..E4 configs: components, derivations, covariates)
```

Outputs that the manuscript consumes (Figure 2, Table 2) are copied to `outputs/figures/` and
`results/` as needed. **Commit `nhanes/raw/`** (or attach as a release artifact) so collaborators reproduce
without re-hitting CDC.

### 3.2 Download + cache (the fix for offline reproducibility)

- Use **`nhanesA::nhanes(tbl, translated = FALSE, cleanse_numeric = TRUE)`** — raw integer codes are
  **stable across cycles** (decode later), and `cleanse_numeric` pushes 7/77/9/99 "Refused/Don't-know"
  sentinels to `NA` in continuous covariates.
- **Download-if-missing cache** to `nhanes/raw/`:
  ```r
  raw_dir <- here::here("nhanes","raw"); dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  get_table <- function(tbl) {
    f <- file.path(raw_dir, paste0(tbl, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    x <- nhanesA::nhanes(tbl, translated = FALSE, cleanse_numeric = TRUE); saveRDS(x, f); x
  }
  ```
- **Per-component table names** (do NOT apply one suffix list to all components):
  ```r
  tbl_name <- function(comp, suf) if (comp == "BPX" && suf == "J") "BPXO_J" else paste0(comp, "_", suf)
  ```
  Suffixes `_E _F _G _H _I _J` = 2007-08 … 2017-18. **`BPX → BPXO_J`** in 2017-18 (variables
  `BPXOSY1/BPXODI1`) — harmonize to common SBP/DBP after download. **`RIDRETH3`** (adds Non-Hispanic Asian)
  exists only from `_G` (2011-12) on; `_E/_F` have `RIDRETH1` only → `coalesce(RIDRETH3, RIDRETH1)` or use
  `RIDRETH1`'s 4-level scheme throughout for cross-cycle comparability.
- **Provenance manifest** (`raw/manifest.csv`): table, nrow, ncol, timestamp, `packageVersion("nhanesA")`,
  `digest::digest()` hash. Pin `nhanesA` with **`renv`** (current CRAN 1.4.1, 2025-10-14). For maximum
  determinism, optionally point `NHANES_TABLE_BASE` at the CCB-HMS NHANES Docker/SQL container.

### 3.3 Merge, decode, derive

- Merge components **within a cycle** by `SEQN` (`full_join`), add a `CYCLE`/`SDDSRVYR` column, then
  `bind_rows` across cycles.
- **Decode** from cached codebooks (`nhanesTranslate(tbl, cols, data=)`) or explicit integer recodes (e.g.
  `DMDEDUC2`: 1=<9th … 5=college grad+; 7/9→NA). This removes the brittle string-`case_when` matching.
- Derive exposure/outcome/covariates per the chosen example (§2/§4). Set explicit factor reference levels
  (`No < Yes`) before imputation/estimation.

### 3.4 Survey weights & design (the correctness core)

- **Weight:** analysis uses MEC variables (BP, BMI) → **`WTMEC2YR`** (not the interview weight). For the
  six pooled cycles, **`WTMEC_POOLED = WTMEC2YR / 6`** (CDC rule: divide the 2-yr weight by the number of
  pooled 2-yr cycles). The denominator must equal the number of cycles actually present; if a covariate
  forces dropping a cycle, choose a cross-cycle-available variable instead, or adjust the denominator.
  *(UPF→obesity would use `WTDRD1`; PFAS would use `WTSA2YR` — never mix weight types in one design.)*
- **Design:**
  ```r
  options(survey.lonely.psu = "adjust")
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_POOLED,
                   nest = TRUE, data = dat_all_MEC)   # nest=TRUE is mandatory (SDMVPSU repeats across strata)
  ```
- **Subpopulation = subset the DESIGN, not the data** (domain estimation). Build `des` on **all
  MEC-examined records** (`WTMEC2YR > 0`), add a logical `inpop` indicator (e.g. female & age≥20 &
  ever-pregnant & not prevalent-diabetic & non-missing exposure/outcome), then `des_sub <- subset(des,
  inpop)`. Deleting rows first removes strata/PSUs the variance estimator needs and biases SEs downward.
  For the cross-fitted arm, fit nuisances on the domain but compute the Eq-8 variance over the design with
  the domain indicator. Verify design df and min PSUs/stratum remain adequate after restriction.

---

## 4. Covariates — a defensible, non-toy set

**Principle:** condition only on **pre-exposure** common causes; prefer measured **proxies** of strong
unmeasured confounders (with sensitivity analyses) over omission; **never adjust descendants of the exposure
or outcome**; collapse high-cardinality covariates to protect positivity under a rare exposure; keep design
variables (strata/PSU/weight/cycle) distinct from substantive confounders.

**Worked case — GDM → hypertension** (the richer set generalizes to any example):

*Include (core minimal sufficient set):* age `RIDAGEYR`, race/ethnicity `RIDRETH1`/`RIDRETH3`, education
`DMDEDUC2`, poverty-income ratio `INDFMPIR`, survey cycle `SDDSRVYR` (secular trend), BMI `BMXBMI`
(*proxy for unmeasured pre-pregnancy adiposity — partial-mediator caveat, sensitivity analysis*), number of
pregnancies `RHQ160`, parity `RHQ171`, age at first birth `RHD180`, family history of diabetes `MCQ300c`
(best available genetic proxy — NHANES has **no** family-HTN variable), smoking `SMQ020`+`SMQ040`,
menopausal status `RHQ031`.

*Richer set adds:* alcohol `ALQ*`, physical activity `PAQ*`, diet quality HEI-2015 (`DR1TOT`+`FPED`, **uses
`WTDRD1` → secondary/sensitivity model to protect n/positivity**), age at last birth `RHD190`, family
history of early CVD `MCQ300a`. → ~14–16 covariates across demographic / SES / anthropometric / reproductive
/ familial / behavioral domains — decisively past the 5-variable toy.

*EXCLUDE (mediators / colliders / post-treatment — over-adjustment bias):* subsequent/current type-2
diabetes `DIQ010` (mediator; use **only as a cohort exclusion**, never a covariate), post-pregnancy weight
gain/Δ-BMI (mediator), antihypertensive meds/treated BP `BPQ040A/BPQ050A` (**belongs in the outcome
definition**), prevalent CVD `MCQ160*` (collider), fasting glucose/HbA1c/insulin `LBXGLU/LBXGH/LBXIN`
(metabolic mediators), same-pregnancy preeclampsia/gestational HTN (collider), post-index OC/HT use.

*Effect modifiers (report stratum-specific RDs):* race/ethnicity, BMI category, family history of diabetes,
age / time-since-pregnancy, menopausal status, survey period.

*Positivity:* plot `g(W)` by arm; report realized min/max `g` and fraction near the bound (`run_estimators`
already records `g_*_min/max`); collapse parity/number-of-pregnancies to ordered categories and bin age;
rely on the cross-fitted propensity (less overfit in sparse cells); report a common-support-trimmed
sensitivity analysis alongside the full-sample estimate.

**Outcome hardening (all examples with a BP outcome):** define hypertension as **averaged SBP ≥130 or DBP
≥80 OR self-reported diagnosis (`BPQ020`) OR antihypertensive use (`BPQ040A`)** — not a single ≥140/90
reading — so treated-controlled hypertensives are not misclassified.

### 4.1 Example-specific confounders (beyond the shared core: age, sex where applicable, race/ethnicity, education, PIR, cycle)

- **E1 sleep→obesity:** physical activity (`PAQ*`), total energy intake (`DR1TKCAL`), alcohol (`ALQ*`),
  smoking (`SMQ020/040`), depression (`DPQ` — a *confounder* of sleep→obesity, not a mediator here),
  marital status, chronic conditions affecting sleep. Stated caveat: reverse causation (obesity→short sleep).
- **E2 food insecurity→depression:** marital status, employment, household size, health insurance
  (`HIQ011`), chronic-disease burden (`MCQ`), BMI (SES core already covered by PIR/education).
- **E3 e-cig→hypertension:** **combustible smoking** (`SMQ020/040`, pack-years `SMD650`) is the key
  confounder/wrinkle; plus BMI, alcohol, physical activity, diabetes (`DIQ010`, a *pre-exposure* confounder
  of HTN here), SES. Sensitivity: restrict to never-combustible-smokers.
- **E4 GDM→hypertension:** the detailed pre-exposure reproductive/metabolic set above (parity, number of
  pregnancies, age at first birth, family history of diabetes, menopausal status, BMI proxy), with the
  strict mediator/collider EXCLUDE list (subsequent diabetes, weight gain, lipids/glucose, preeclampsia).

Each example's include/exclude lists live in the registry; **verify no post-exposure variable enters**
(never condition on descendants of the exposure or outcome). Note `DIQ010` is a legitimate pre-exposure
confounder in E3 but a **mediator (exclude / use only as a cohort exclusion) in E4** — the same variable
plays different causal roles across examples, which the registry must encode per example.

---

## 5. Analysis outputs

Produced **for each of E1–E4**, plus a **cross-example comparison** (`cross_example_summary.csv` + figure)
to support the lead/secondary choice. Once the lead is chosen, its forest plot/table become the manuscript's
Figure 2 / Table 2; a secondary example (if kept) goes to the Web Appendix.

- **Forest plot (→ manuscript Figure 2 for the lead):** ATE (risk difference) with design-aware 95% CIs
  across the five arms (+ GLM benchmark), showing the CI widening Non → Partially → Fully → CF.
  Black-on-white, alt text.
- **Table (→ manuscript Table 2 for the lead):** estimate, SE, 95% CI, and design df by arm; primary CF
  headline.
- **Cross-example comparison (decision aid):** the four examples' arm estimates/CIs side by side, with
  overlap summaries (min/max `g`, exposed n, design df) — the basis for §8 decision 1.
- **Overlap diagnostic figure (Web):** propensity-score histograms by exposure arm, single-fit vs CF.
- **Web "applied echo of Figure 1" (Web F):** single-fit Fully-Aware vs Fully-Aware-CF across the L1→L4
  ladder on the NHANES data.
- **Web Appendix E:** survey-weighted Table 1 (descriptives by exposure), variable coding, imputation note.
- **Caveats to state:** cross-sectional temporality / possible reverse causation; single imputation
  under-propagates variance (note multiple imputation with Rubin's rules as the principled extension);
  positivity/overlap; calibrated `WTMEC2YR` treated as `1/π` makes Eq 8 conservative.

---

## 6. Scripts to write (replacing `codes/nhanes{1,2,3,4}*.R`)

All under `nhanes/`, sourcing `R/config.R` for paths and `R/estimators.R` for the engine. The
scripts are **example-parameterized**: each loops over (or is invoked per) the four registry entries E1–E4.

- `examples.R` — the **example registry**: a list of four configs, each with `{id, components, subpop_fn,
  exposure_fn(df), outcome_fn(df), covariates_include, covariates_exclude_note, weight = "WTMEC2YR"}`.
  Single source of truth for what differs across E1–E4.
- `00_download_cache.R` — union of all components needed across E1–E4 (DEMO, BMX, BPX/BPXO, BPQ, RHQ, DIQ,
  SLQ, FSQ, DPQ, SMQ/ECQ, PAQ, ALQ, MCQ, HIQ, DR1TOT); per-component/-cycle `get_table()` + manifest;
  handles `BPXO_J` and `RIDRETH3` availability; caches codebooks. **Idempotent / offline after first run.**
- `01_build_analytic.R` — for each example: merge, decode (integer recodes), apply `exposure_fn`/`outcome_fn`
  /covariates, build the eligibility indicator (**no row deletion**); save `analytic/<example>_analytic.rds`.
- `02_impute.R` — per example: single imputation (`mice`, documented seed/methods);
  save `analytic/<example>_imputed.rds`. *(Stretch: `m>1` MI + Rubin's rules on Eq 8.)*
- `03_run_estimators.R` — per example: build `svydesign` (full MEC, `nest=TRUE`), `subset()` to the domain,
  call the **generalized** `run_estimators()` for all five arms + the GLM benchmark; save `results/<example>/`.
- `04_diagnostics.R` — per example: descriptives (Table 1), positivity/overlap (`g(W)` histograms, min/max,
  trim sensitivity), design-df / min-PSU checks; save `results/<example>/`.
- `05_figures.R` — per example: forest plot of arms + overlap plot → `nhanes/figures/<example>/`.
- `06_compare.R` — assemble `results/cross_example_summary.csv` (all four × all arms) + a cross-example
  comparison figure, to support the lead/secondary choice (§8).
- `run_all.R` — sources `examples.R` then 00→06 over E1–E4.

**First refactor (blocking):** generalize `R/estimators.R::run_estimators()` to accept an arbitrary
covariate set `W` (currently hard-coded `L1..L4`); add a thin NHANES adapter mapping each example's analytic
frame to the `obs` contract (`Y, A, W…, strata=SDMVSTRA, cluster=SDMVPSU, weight=WTMEC_POOLED`).

---

## 7. Compute & cost

Five arms × one rich SL library on n≈10k, with PSU-level cross-fitting, is modest per example (minutes to
low-tens-of-minutes locally; the SL+CF arm dominates). **Four examples ≈ 4× that** — still comfortably
local; parallelize across examples if desired. The optional ladder echo (single vs CF across L1–L4) is the
expensive add-on; run it once, for whichever example becomes the lead, and cache. No cluster needed.

---

## 8. Open decisions for the author

1. **Example — decided: analyze ALL FOUR (E1–E4), then choose** the lead (and optional secondary) **after
   comparing** `results/cross_example_summary.csv` + the per-example overlap diagnostics. The decision
   criteria: overlap quality, how cleanly the Non→CF CI-widening shows, whether single-fit/CV diverges from
   CF, sample size after restriction, and subject-matter novelty.
2. **SL library for the application** (single rich library for the main forest plots vs. adding the full
   L1→L4 ladder echo in Web F for the chosen lead).
3. **Imputation:** keep single imputation (with caveat) or upgrade to multiple imputation + Rubin's rules.
4. **Commit `nhanes/raw/` to git** (≈tens of MB) vs. attach as a release artifact / regenerate on demand.

---

## Implementation progress (2026-06-07)

- ✅ **Refactor (blocking item) done & validated.** `R/estimators.R::run_estimators()` now takes
  `W_cols` (arbitrary covariates), `nest`, and `inpop` (sub-population). Sub-population variance is computed
  by **building the design on the full sample and `subset()`-ing to the domain** (per EpiMethods
  surveydata8) — smoke-tested: in domain mode the clustered arms' df = (#PSU−#strata) of the **full**
  design, not the reduced sub-sample. Simulation path (default args) is structurally unchanged.
  `make_cf_folds` fixed to key on the (stratum, PSU) pair so it is correct for NHANES `SDMVPSU` (repeats
  across strata); identical to before when ids are globally unique.
- ✅ **Registry + downloader written.** `nhanes/examples.R` (E1–E4 configs) and
  `nhanes/00_download_cache.R` (download-if-missing → `nhanes/raw/<TABLE>.rds` + `manifest.csv`,
  `translated=FALSE`, BPX→`BPXO_J` handled).
- ✅ **Data layer de-risked.** Validation download of 6 tables succeeded against live CDC (nhanesA 1.4);
  **`BPXO_J` confirmed to carry `BPXOSY1/BPXODI1`** (the legacy `BPX_J` request lost 2017–18 BP). Full
  union download (16 components × 6 cycles) run to cache `nhanes/raw/`.
- ✅ **All 96 tables cached** (16 components × 6 cycles, 0 missing) + columns verified against real data:
  **E3 e-cig items (`SMQ900/905`) exist only 2015–18 → E3 uses 2 cycles, pooled weight `WTMEC2YR/2`**
  (E1/E2/E4 = all six, `/6`); sleep var `SLD010H` (2007–14) → `SLD012` (2015–18), harmonize + clean 77/99;
  race from **`RIDRETH1`** (`RIDRETH3` absent 2007–10); BP `BPXSY*`→`BPXOSY*` at 2017–18; E2/E4 vars all
  present in 6 cycles. Per-example `cycles` field now in `examples.R`.
- ⏭ **Next:** `01_build_analytic.R` (per-example exposure/outcome/covariate derivations from cached raw;
  harmonize BP & sleep; pooled weight per `length(cycles)`; build the `inpop` indicator WITHOUT deleting
  rows) → `02_impute` → `03_run_estimators` (full design + `subset` to domain, 5 arms) → `04_diagnostics`
  → `05_figures` → `06_compare`.

## 9. Key sources

- nhanesA: `nhanes()`/`nhanesTranslate()`/`nhanesSearchVarName()` refs (rdrr.io/cran/nhanesA), the
  *Introducing nhanesA* vignette, CRAN news (v1.0 Docker revision), CRAN 1.4.1 landing page, and the
  *Database* 2024 paper (PMC11020206).
- CDC NHANES **Weighting tutorial** and **Analytic Guidelines** (use `WTMEC2YR`; pooled weight = `/`number
  of 2-yr cycles; combine `SDMVSTRA`/`SDMVPSU`; subpopulation must not delete records):
  wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx ; .../analyticguidelines.aspx.
- Lumley, survey "Estimates in subpopulations" (domain) vignette — subset the design, not the data.
- Example literature (2022–2025): sleep–obesity (BMC Public Health 2025, PMC12010617); e-cig–hypertension
  (PMC11574958, 2024); food insecurity–depression (Nutrients 2022, PMC9370686); physical activity–
  depression (PMC11224452, 2024); UPF (J Nutrition 2022, NOVA-NHANES); PFAS–gout (Front Public Health 2025,
  PMC11847820).
