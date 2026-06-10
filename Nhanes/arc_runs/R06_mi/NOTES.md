# R06_mi — Multiple imputation (m=40) headline NHANES Table-2 intervals

## Purpose
The locked NHANES analysis (`Nhanes/R/02_impute.R`) imputes missing covariates
with **single imputation (mice m=1)**, which understates variance (a stated
caveat in the paper). This enhancement run replaces it with **multiple
imputation, m=40**, for the **certified primary estimator** and combines the
per-imputation estimates/SEs by **Rubin's rules** to produce MI-based Table-2
intervals on **all four examples**. The headline run uses **m=40**, following
White, Royston & Wood (2011): m >= the percentage of incomplete cases (E4 = 29%).

- **Certified primary** = the 3-learner near-Donsker library
  `LIB = c("SL.glm","SL.earth","SL.glmnet")` — the **same `LIB` as
  `Nhanes/R/03_run_estimators.R`** (the L2-equivalent rung). The reported
  headline arm is **Fully-Aware-CF**, but all five arms returned by
  `run_estimators` are pooled so the full Table-2 column is available with MI.
- We deliberately **do NOT MI the full L1–L4 ladder** — that would multiply the
  expensive deep rungs (L3/L4) by 10. Only the certified-primary `LIB` is run.

## What it does (per example = per array task)
1. **Impute m=40** completed datasets with mice, reusing `02_impute.R`'s exact
   setup (same SEED=20260607, maxit=5, A/Y as predictors-never-imputed,
   character→factor) but with `m=40` instead of `m=1`
   (`impute_m()` in `mi_helpers.R`).
2. **Run `run_estimators(... LIB, inpop=domain, nest=TRUE)`** on each completed
   dataset (one-hot encoding the imputed domain covariates exactly as
   `03_run_estimators.R` does, via `encode_domain()`). Each imputation gets its
   own deterministic CF/SL split seed.
3. **Pool per arm by Rubin's rules** (`rubin_pool()`): within variance
   `Ubar=mean(se^2)`, between variance `B=var(b)`, total `T=Ubar+(1+1/M)B`,
   pooled `se=sqrt(T)`, **t reference with Barnard–Rubin (1999) df**, plus FMI
   and RIV. (This reproduces what `mice::pool()` computes.)
4. CI = `b ± qt(0.975, df_BR)·se`; `signif = CI excludes 0`.

## NON-DESTRUCTIVE guarantees
- Only **sources** the canonical engine `codes/estimators.R` (read-only). All new
  logic (`impute_m`, `rubin_pool`, `encode_domain`) lives in this folder's
  `mi_helpers.R`, reusing the exposed `run_estimators` building block.
- Edits **nothing** in `codes/*.R` or `Nhanes/R/*.R`.
- Does **not** read or overwrite `Nhanes/analytic/*_imputed.rds` (the locked
  single-imputation frames). It reads only `*_analytic.rds` and imputes fresh
  in-memory.
- Outputs go to this run's own folders (see below); the locked
  `Nhanes/nhanes_output/results/*` and `Nhanes/results/*` are untouched.

## Files
- `run.R`        — the driver (array task = example; SMOKE mode built in).
- `mi_helpers.R` — `impute_m()`, `rubin_pool()`, `encode_domain()` (the new logic).
- `submit.slurm` — ARC array job (1–4), 32 cpus, 64G, 3 h, mirrors `nhanes.slurm`.
- `NOTES.md`     — this file.

## How to SMOKE-TEST (locally, ~2–6 min, before submitting)
Runs the two **borderline** examples (E2 significant, E3 null) with **M=2**:
```powershell
# Windows (R 4.5) — from anywhere
$env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/codes"
$env:SMOKE="1"; $env:MI_M="2"; $env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/Nhanes/arc_runs/R06_mi/run.R"
```
```bash
# Linux / ARC login node (R 4.4)
export REPO_ROOT=$PWD SIM_CODE=$PWD/codes SMOKE=1 MI_M=2 SLURM_CPUS_PER_TASK=2
Rscript Nhanes/arc_runs/R06_mi/run.R
```
SMOKE writes `results/arc/R06_mi_summary_SMOKE.csv` and
`results/arc/R06_mi_primary_SMOKE.csv` and per-example
`Nhanes/nhanes_output/arc_runs/R06/mi_{E2,E3}.rds`. The `_SMOKE` suffix means it
never collides with the full-run CSVs.

## How to SUBMIT (full run on ARC)
From the repo root (SLURM submit dir / scratch):
```bash
sbatch Nhanes/arc_runs/R06_mi/submit.slurm
```
4 array tasks (one per example), M=10. After it finishes, tighten walltime with
`seff <jobid>`.

## Expected output
**Per-example RDS** → `Nhanes/nhanes_output/arc_runs/R06/mi_E{1,2,3,4}.rds`
each containing: `pooled` (per-arm MI estimate/se/df/CI + `b_imp_sd`,
`within_se`, `between_sd`, `riv`, `fmi`, `signif`), `per_imp` (the M×arm raw
estimates), `imputed_covs`, `n_domain`, `A_prev`, `Y_prev`.

**Summary CSVs** → `results/arc/`:
- `R06_mi_summary.csv`  — all arms × all examples (MI-pooled).
- `R06_mi_primary.csv`  — the certified-primary (Fully-Aware-CF) MI Table-2 column.

**Manifest** → `Nhanes/nhanes_output/arc_runs/R06/manifest/manifest_*.rds`
(R version, package versions, git rev, seeds, timing).

Key numbers to inspect: for each example, the **Fully-Aware-CF** row's
`b / se / lcl / ucl / df / fmi` and `signif`. Compare `se` and CI width against
the locked single-imputation Table-2 (Fully-Aware-CF rows in
`Nhanes/nhanes_output/results/nhanes_results_summary.csv`, the L2_smooth rung,
which is the same `LIB`). **MI SE should be slightly LARGER** (between-imputation
variance added) and CI slightly WIDER.

## Aggregation
The driver writes the combined CSVs directly (no separate aggregate step needed
when run per-task: each task appends only its own example's rows to a
per-example RDS, and the combined CSV is rebuilt by re-running the driver with
all examples, OR assemble across tasks by reading the four
`mi_E*.rds` and `rbind`-ing their `$pooled`). The per-example RDS are the
canonical artifacts; the CSVs are the convenience tables matching the locked
Table-2 column names (`example,label,method,b,se,df,lcl,ucl` + MI extras).

## Decision rule (smoke-gate)
**Smoke E2 and E3 first** (the borderline examples — E2 is the significant one,
E3 the null one). Proceed to the full 4-example submit only if the smoke result
is qualitatively consistent with single imputation:

- **PROCEED** if, for the certified primary (Fully-Aware-CF):
  - **E2 MI CI still excludes 0** (`signif = TRUE`, effect stays significant), and
  - **E3 MI CI still includes 0** (`signif = FALSE`, effect stays null), and
  - MI SEs are **≥** the single-imputation SEs (within ~1.5×) — i.e. MI **widens**
    intervals modestly, as expected (between-imputation variance), and
  - `fmi` is sane (0 ≤ fmi < ~0.5) and the run completes without estimator errors.

- **STOP-and-report** (do NOT launch the full run; flag for review) if any of:
  - E2 flips to **non-significant** or E3 flips to **significant** under MI
    (would change a paper conclusion — investigate before scaling up);
  - MI SE comes out **smaller** than single-imputation SE (Rubin pooling should
    not shrink variance — points to a pooling bug or degenerate between-variance);
  - `fmi ≥ ~0.5` or `df_BR` collapses to a tiny value (imputation model unstable);
  - any imputation throws (mice non-convergence, `run_estimators` error, or
    one-hot column mismatch across imputations).

Reasoning: at M=2 the point estimates and significance pattern should already
match single imputation closely (only the variance inflates a little). A flip in
significance or a SE that shrinks is the signal that something is wrong with the
MI/pooling wiring, and scaling to M=10 × 4 examples would just burn ARC time.

## Dependencies
- Engine: `codes/config.R` (seeds; sourced indirectly is not required — the
  driver sets its own paths but reuses SEED=20260607), `codes/estimators.R`
  (`run_estimators` + building blocks).
- Data: `Nhanes/analytic/E{1,2,3,4}_analytic.rds` (must exist; produced by the
  upstream `Nhanes/R/01_*` build). NOT the `_imputed.rds` frames.
- R packages (same as nhanes.slurm): `SuperLearner, tmle, survey, surveyCV,
  earth, glmnet, gam, ranger, mice`. Verified locally: mice 3.18, earth 5.3.4,
  glmnet 4.1.10.

## Caveats / assumptions
- **Assumes the `_analytic.rds` frames carry `attr(obs,"covs")`, an `inpop`
  logical, and design columns `SDMVSTRA / SDMVPSU / WTMEC_POOLED / A / Y`** —
  verified for all four examples (E1–E4 all have missing covariates, so MI
  actually does work; if an example had no missingness `impute_m` returns M
  identical copies and `B=0`, so MI reduces to single imputation cleanly).
- **mice is run on the survey-design-IGNORED domain frame** (same as the locked
  `02_impute.R`): the imputation model does not itself incorporate the weights /
  PSU structure. This matches the paper's existing single-imputation procedure;
  MI here only adds Rubin's between-imputation variance on top of the same
  design-based within-imputation variance from `run_estimators`. (A
  survey-aware imputation model would be a further, separate extension.)
- **Barnard–Rubin df** uses the per-imputation design df reported by
  `run_estimators` (the Fully-Aware-CF design df ≈ #PSU − #strata, ~30–94) as
  the complete-data df, so the pooled df never exceeds the design df.
- Each imputation uses an independent CF/SL split seed, so the reported MI SE
  also absorbs split-induced variability across imputations (a mild extra
  conservatism, consistent with the run-design philosophy of `nhanes_arc.R`).
- The user must **verify the smoke CSV** (`R06_mi_primary_SMOKE.csv`) against the
  smoke-gate rule above before `sbatch`.
