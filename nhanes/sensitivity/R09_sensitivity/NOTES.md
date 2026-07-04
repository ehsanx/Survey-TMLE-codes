# R09_sensitivity — Causal-validity sensitivity analyses (LOCAL run)

## Purpose

Pre-empt an **identification-grounds major revision** by providing reviewer-ready
re-runs of the **primary** NHANES estimator on **modified adjustment
sets**, plus an empirical **design-validity check**. Each sensitivity = the same
primary pipeline run on a different covariate set, reported as the
**ATE / CI delta vs the main (full-adjustment) analysis**.

Four deliverables:

| id | what | covariate change | why |
|----|------|------------------|-----|
| (a) `E4_noBMI` | E4 GDM → hypertension, **BMI omitted** | drop `bmi` | BMI is a *pre-pregnancy adiposity PROXY*; if it sits on the GDM→HTN path it is a **mediator** (drop), not a confounder. Pre-specified mediator-vs-confounder sensitivity. |
| (b) `E2_noBMI` | E2 food insecurity → depression, **BMI as mediator** | drop `bmi` | BMI may **mediate** food-insecurity → depression; show estimate with BMI dropped vs the confounder-treatment in the main analysis. |
| (c) `E1_noPHQ` | E1 short sleep → obesity, **PHQ-9 dropped** | drop `phq` | Alternate adjustment: depression (PHQ-9) is treated as a confounder in the main analysis; show robustness to dropping it. |
| (d) design check | **SDMVSTRA cross-cycle non-overlap** | none (tabulation) | Verify no `SDMVSTRA` value spans >1 NHANES cycle in the pooled design → the pooled stratification / `nest=TRUE` is valid. A tabulation, **not** an estimator run. |

**Certified primary** = `LIB = c("SL.glm","SL.earth","SL.glmnet")` (the
near-Donsker library used by `nhanes/03_run_estimators.R`). Headline arm =
**Fully-Aware-CF**; all five arms are kept. Seed `20260607` (== `03_run_estimators.R`).

This is a **LOCAL** run: a covariate-set swap loop + one tabulation. It runs on a
**laptop in a few minutes — NO SLURM needed**. An optional `submit.slurm` is
included only if you prefer to run it on ARC for uniformity.

## Non-destructive design

- SOURCES the canonical engine `R/estimators.R` **read-only** (exposes
  `run_estimators` + building blocks). Edits **nothing** in `codes/` or `nhanes/`.
- Reads the already-imputed frames `nhanes/analytic/<id>_imputed.rds`. The
  adjustment-set swap is done purely by overriding the `attr(obs,"covs")`
  attribute (the same hook `01_build_analytic.R` sets and `03` consumes) — the
  underlying data columns are left untouched. Helper: `set_covs()` in
  `sensitivity_helpers.R`.
- All new logic lives in **`sensitivity_helpers.R`** (this folder).
- Outputs go to this run's **own** folders and never clobber the locked headline
  results (`results/`, `nhanes/results/`).
- **The MAIN comparison is recomputed here** under the SAME library/seed (it is
  cached per example and reused across that example's cells), so the delta is
  apples-to-apples — it is NOT differenced against the locked Table-2 CSV (whose
  seed / library could differ).

## Files

- `run.R` — the driver (loops the 3 sensitivity cells, runs the design check, writes CSVs + manifest).
- `sensitivity_helpers.R` — `set_covs`, `encode_domain`, `run_primary`, `ci_delta`, `strata_cycle_check`.
- `submit.slurm` — OPTIONAL convenience submit (single task, no array; 32 cpus, 32G, 1h).
- `NOTES.md` — this file.

## How to smoke-test (do this first, locally)

```bash
# from the repo root <repo-root>
SMOKE=1 Rscript nhanes/arc_runs/R09_sensitivity/run.R
```

SMOKE mode runs only cell **(a) E4-no-BMI** plus the **(d) design check**, with a
**GLM-only** library (`SL.glm`). Wall time ~**3–6 min** on a laptop (it runs both
the full-covs MAIN and the no-BMI SENSITIVITY 5-arm `tmle` suites on the ~11k-row
E4 domain). It writes `*_SMOKE.csv` files and a `sens_E4_noBMI.rds`; delete those
after validating. (PowerShell: `$env:SMOKE='1'; Rscript ...`.)

Assumes R 4.4/4.5 with the same packages as the rest of the project
(SuperLearner, tmle, survey, surveyCV, earth, glmnet — all confirmed present).

## How to run (full, LOCAL — the intended mode)

```bash
# from the repo root
Rscript nhanes/arc_runs/R09_sensitivity/run.R
```

Runs all three sensitivity cells (E4-no-BMI, E2-no-BMI, E1-no-PHQ) on the full
3-learner library + the design check. Expect ~**15–30 min** on a laptop (six
5-arm `tmle` suites total: a MAIN + a SENSITIVITY fit per example).

## How to submit (OPTIONAL, ARC)

```bash
# from the SLURM submit dir / scratch (repo root)
sbatch nhanes/arc_runs/R09_sensitivity/submit.slurm
```

Single task (no array). Mirrors the proven ARC environment
(`YOUR_ALLOCATION`, gcc/9.4.0 + r/4.4.0, `R_LIBS`). Resources: 32 cpus, 32G, 1h
(generous — the job is light and largely serial).

## Expected output

Per-cell RDS → `nhanes/nhanes_output/arc_runs/R09/sens_<cell>.rds`
Manifest → `nhanes/nhanes_output/arc_runs/R09/manifest/manifest.rds`
Summary CSVs → `results/`:

- **`R09_sensitivity_delta.csv`** — *the headline deliverable*. One row per
  (cell × arm): `b_main`, `lcl_main`, `ucl_main`, `signif_main`, `b_sens`,
  `lcl_sens`, `ucl_sens`, `signif_sens`, `d_b` (= b_sens − b_main), `d_lcl`,
  `d_ucl`, `rel_d_b`, **`conclusion_flip`** (did the CI-excludes-0 conclusion
  change?). Inspect the **Fully-Aware-CF** rows.
- `R09_sensitivity_arms.csv` — all five arms, main vs sensitivity, every cell
  (`example,label,cell,adjust,method,b,se,df,lcl,ucl,signif`) — matches the
  locked Table-2 column conventions.
- `R09_design_check.csv` — per example: `n_strata`, `n_cycles`,
  `max_cycles_spanned`, `n_strata_multicycle`, **`overlap_clean`** (full design)
  and `dom_*` (domain-only).

**Key numbers to read:** in `R09_sensitivity_delta.csv`, the Fully-Aware-CF `d_b`
and `conclusion_flip`; in `R09_design_check.csv`, `overlap_clean` must be `TRUE`
for all four examples.

## Aggregation

None needed — `run.R` is a single local pass that writes the final CSVs directly.
(For consistency with other runs you could `rbind` the `sens_*.rds` `delta`
tables, but the CSVs already contain the full result.)

## Smoke-gate decision rule (STOP-and-report vs proceed)

After the smoke run, inspect `results/R09_design_check_SMOKE.csv` and the
E4-no-BMI delta:

1. **Design check** — every example must show `overlap_clean = TRUE`
   (`max_cycles_spanned = 1`). **If any `overlap_clean = FALSE`** → a `SDMVSTRA`
   stratum spans >1 cycle in the pooled design: **STOP and report** — the pooled
   stratification (and `nest=TRUE`) would be mis-specified; the analytic build
   needs a cycle-aware stratum re-keying before any estimate is trusted. This is
   a hard gate.
   *(Validated 2026-06-07: all four examples CLEAN, `max_cycles_spanned = 1`.)*
2. **Sanity of the swap** — the E4-no-BMI cell must produce finite, in-range
   ATEs for all arms and a Fully-Aware-CF SE within ~2× the main-fit SE. **If the
   no-BMI fit returns `NaN`/non-finite estimates, an `Inf` design df, or an SE
   that explodes (>5× the main SE)** → the encoder/positivity broke on the
   reduced covariate set: **STOP and report** rather than launching the full run.
3. Otherwise (design clean + finite, comparable estimates) → **proceed** to the
   full local run. A large `|d_b|` or `conclusion_flip = TRUE` is a *substantive
   finding to report*, not a reason to stop — that is exactly the
   mediator-vs-confounder signal the reviewers want quantified.

## Dependencies

- `R/config.R` (seeds; sourced transitively is not required — `run.R` sets
  its own seed), `R/estimators.R` (the engine), `R/learners.R` (only if
  you swap in ladder libraries; the default LIB needs no custom wrappers).
- `nhanes/analytic/{E1,E2,E4}_imputed.rds` for the sensitivity cells and
  `{E1,E2,E3,E4}_imputed.rds` for the design check (produced by
  `nhanes/01_build_analytic.R` → `02_impute.R`; already present).
- R packages: SuperLearner, tmle, survey, surveyCV, earth, glmnet (+ mice/ranger
  only listed in the manifest, not required to run).

## Caveats

- **Single imputation.** Uses the existing `*_imputed.rds` (mice m=1), so SEs are
  the single-imputation SEs (under-stated, same caveat as the headline Table 2).
  The MI extension lives in R06; this run is a covariate-set sensitivity, not an
  imputation sensitivity, so single-imputation is the right matched comparison.
- **Dropping a mediator does not, by itself, identify the total effect** — it
  changes the estimand's adjustment set. The deliverable is the *magnitude* of
  the resulting ATE/CI shift (does the conclusion survive treating BMI / PHQ-9
  as off the adjustment set?), which is what an identification reviewer asks for.
- **Positivity.** Removing a covariate generally *improves* overlap (smaller
  design matrix), but the run still reports the FA single-fit g-range
  (`g_fa_main`, `g_fa_sens` in the per-cell RDS) so any overlap change is visible.
- The design check is a **necessary** condition for the pooled stratification, not
  a full survey-design audit; it confirms strata are cycle-nested, which is the
  specific concern (`nest=TRUE` validity) raised for the pooled multi-cycle design.
