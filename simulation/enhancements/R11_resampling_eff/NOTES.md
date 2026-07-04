# R11_resampling_eff — Resampling-SE check + IPW/svyglm efficiency benchmark

Web Appendix D. Self-contained, **non-destructive** ARC run package. Sources the
canonical engine (`R/config.R`, `dgp.R`, `estimators.R`, `diagnostics.R`,
`learners.R`) **read-only**; all new logic lives in `helpers_R11.R`. No canonical
file is edited; outputs are isolated and never clobber the locked results.

## Purpose
Two deliverables, both at rungs **L1_param** (run FIRST) and **L3_adaptive**,
both scenarios (`standard`, `R1`):

1. **(a) Resampling-SE check (defends Eq 8).** Corroborate the engine's
   design-**linearization (Taylor)** SE against a replication-based **survey
   jackknife (JKn, delete-one-PSU)** SE computed on the **SAME** influence
   function — for the single-fit **Fully-Aware** arm and the cross-fitted
   **Fully-Aware-CF** arm. Built with `survey::as.svrepdesign(type='JKn')` +
   `survey::withReplicates(rep, quote(sum(.weights*eif)/sum(.weights)))`. Because
   `.eif_from_tmle()` already centers the EIF (subtracts psi), the replicate SE of
   its weighted Hájek mean is the jackknife SE of psi-hat — the apples-to-apples
   partner of the Taylor SE. **Expect jk_lin_ratio ~ 1.**

2. **(b) Efficiency benchmark.** A survey-weighted **IPW / g-computation (AIPW)**
   baseline built only from `svyglm` (quasibinomial) nuisances, scored with the
   **same** design-EIF machinery (`.se_des`). Compare coverage and CI width to the
   cross-fitted TMLE. **Expect TMLE >= IPW efficiency (CI width ratio <= ~1) with
   valid coverage for both.**

500 reps/cell. L1 runs first; only trust/extend to L3 if L1 agrees.

## Files
- `run.R` — driver. Reuses `run_estimators()` for the TMLE arms, recomputes the
  CF-arm EIF from the exposed building blocks (`make_cf_folds`, `.sl`,
  `.eif_from_tmle`), jackknifes the FA and CF EIFs, and runs the IPW baseline.
- `helpers_R11.R` — NEW logic only: `jk_se_on_eif()` and `ipw_svyglm_ate()`.
- `aggregate.R` — combines per-task RDS into the summary CSV.
- `submit.slurm` — ARC array job (resources mirror `codes/sim.slurm`).

## Output column conventions
Per-rep rows match the engine: `method, b, se, df` (se = Taylor), plus `se_jk,
df_jk` (jackknife SE/df; populated only on Fully-Aware & Fully-Aware-CF), plus
`n, df_design, w_cv`. Methods: the five engine arms + `IPW-svyglm`.

`aggregate.R` writes a long CSV with a `kind` column:
- `kind="taylor"` — standard summary (`bias, emp_sd, mean_se, se_ratio, coverage,
  mcse_cov`) per method, plus `mean_ci_width` and `ci_width_ratio_vs_CF`
  (efficiency; <1 = narrower CI than Fully-Aware-CF).
- `kind="jack"` — Fully-Aware & Fully-Aware-CF only: `mean_se_taylor,
  mean_se_jack, jk_lin_ratio` (≈1 corroborates Eq 8), `se_ratio_taylor/jack`,
  `coverage_taylor/jack`.
Conventions follow `aggregate_sim.R`: `bias=mean(b)-Psi; emp_sd=sd(b);
se_ratio=mean(se)/emp_sd; coverage` via t df; `mcse_cov=sqrt(cov(1-cov)/nreps)`.

## How to SMOKE-test (local, ~1–3 min on a login node / laptop, R 4.4/4.5)
`SMOKE=1` forces a tiny run: ONE cell (standard / L1_param), **20 reps**, serial.
```
# PowerShell (Windows dev):
$env:SMOKE="1"; $env:SIM_CODE="<repo-root>/codes"; `
  $env:R11_OUT="<repo-root>/sim_output/arc_runs/R11_resampling_eff"; `
  Rscript <repo-root>/simulation/enhancements/R11_resampling_eff/run.R

# bash (ARC login node):
SMOKE=1 SIM_CODE="$PWD/codes" \
  R11_OUT="$PWD/sim_output/arc_runs/R11_resampling_eff" \
  Rscript simulation/enhancements/R11_resampling_eff/run.R

# then aggregate the smoke RDS into the CSV:
Rscript simulation/enhancements/R11_resampling_eff/aggregate.R
```
Inspect `results/R11_resampling_eff_summary.csv` (and the printed table).

### SMOKE-GATE decision rule (what means STOP-and-report vs proceed)
At **L1_param** (parametric nuisances — linearization and jackknife are
asymptotically identical), in the `kind="jack"` rows for **Fully-Aware**:
- **PROCEED** if `jk_lin_ratio` ∈ ~[0.9, 1.1] (Taylor ≈ jackknife) AND the
  jackknife SE is finite for ~all reps (no widespread `NA` in `se_jk`).
- **STOP-and-report** if `jk_lin_ratio` is systematically off (< 0.8 or > 1.25),
  or `se_jk` is `NA` for many reps, or the run errors. That would indicate the
  jackknife replicate design or the EIF-centering comparison is mis-specified —
  fix before spending the full-run budget (and definitely before L3, where any
  L1 discrepancy can only widen).
Note: 20 smoke reps give a noisy ratio; treat ~[0.9,1.1] as the smoke band and
expect tighter agreement at 500 reps. Also sanity-check that `IPW-svyglm`
produced a finite `b`/`se` and coverage in a plausible range.

## How to SUBMIT (full run on ARC)
From the repo root on ARC (= `$SLURM_SUBMIT_DIR` / scratch):
```
sbatch simulation/enhancements/R11_resampling_eff/submit.slurm
```
- Grid = rung × scenario × rep-chunk, `R11_N_REPS=500`, `R11_CHUNK=100` → 5
  chunks/cell → 2 rungs × 2 scenarios × 5 = **20 tasks** (`--array=1-20`).
  Driver orders rungs (L1_param, L3_adaptive): **tasks 1–10 = L1**, 11–20 = L3.
- **Cluster-side gate (optional):** submit `--array=1-10` (L1) first; once the L1
  `jk_lin_ratio` checks out, submit `--array=11-20` (L3).
- After all tasks finish, aggregate once:
  `Rscript simulation/enhancements/R11_resampling_eff/aggregate.R`

## Expected output
- Per-task RDS: `sim_output/arc_runs/R11_resampling_eff/R11_<scenario>_<rung>_chunk###.rds`
  and `manifest_<scenario>_<rung>_chunk###.rds` (seeds, sessionInfo, pkg versions,
  git rev — like `run_sim.R`).
- Summary CSV: `results/R11_resampling_eff_summary.csv` (+ `_combined.rds`).
- Key numbers to inspect:
  - `kind="jack"`: `jk_lin_ratio ≈ 1` for Fully-Aware AND Fully-Aware-CF at L1
    (and L3) — defends Eq 8; `coverage_jack ≈ coverage_taylor`.
  - `kind="taylor"`: `IPW-svyglm` vs `Fully-Aware-CF` — TMLE `coverage` near
    nominal and `ci_width_ratio_vs_CF >= 1` for IPW (TMLE no wider, i.e. ≥ as
    efficient), both with valid coverage.

## Aggregation
`aggregate.R` reads `R11_OUT`, builds the two-block long CSV described above, and
writes to `results/`. Feeds the Web Appendix D table comparing
linearization vs jackknife SE and TMLE vs IPW coverage/efficiency.

## Dependencies
R 4.4.0 (ARC) / 4.4–4.5 (dev). Packages already used by the engine: `SuperLearner,
tmle, survey (≥4.x for as.svrepdesign/withReplicates), surveyCV, earth, gam,
ranger, parallel`. No new packages. `R_LIBS` set in `submit.slurm` to the proven
ARC library.

## Caveats
- **Population is rebuilt per task** (deterministic from `POP_SEED=20260606`),
  identical to `run_sim.R`. `truth_M=2e6` MC truth adds ~1 min to each task's
  startup — unavoidable and harmless.
- **CF-arm EIF is recomputed** in `run.R` to jackknife it (the engine's
  `run_estimators` returns only the FA EIF). The recompute mirrors the engine's
  CF arm exactly (same `make_cf_folds(V=5)`, unweighted `.sl` per fold,
  `g_oof_bound=0.05`, weighted pooled `tmle` targeting), so the jackknifed EIF is
  the same influence function the engine scored. Per-fold RNG inside `.sl` runs
  under the same `clusterSetRNGStream`/`set.seed` as the engine call within the
  rep, so the reconstructed fit matches up to SuperLearner's internal CV draw;
  the jackknife/Taylor comparison is on that single reconstructed EIF, so this is
  internally consistent regardless.
- **IPW baseline is parametric and NOT cross-fitted** (svyglm), by design — it is
  the standard applied survey-IPW/AIPW competitor. It is doubly robust (AIPW
  form) but uses only main-effects `svyglm` nuisances; it is expected to be no
  more efficient than TMLE.
- **Jackknife cost:** `standard` has many PSUs (delete-one-PSU = many replicates),
  so the `withReplicates` pass is the per-rep bottleneck — hence `--time=03:00:00`
  (vs 1 h for the base sim). If a task times out, lower `R11_CHUNK` and raise the
  array count, or restrict the smoke/first submit to L1 only.
- `withReplicates` evaluates `quote(... eif ...)` in the design's data frame, so
  `eif` must be a column of the design data (it is). `jk_se_on_eif()` is wrapped
  in `tryCatch` → a failing rep yields `se_jk=NA` rather than aborting the chunk;
  widespread `NA` is a STOP signal (see gate).
- Coverage uses a **t** reference with the design df (`survey::degf`); the
  jackknife CI uses the jackknife df (`degf(rep_des)`), matching survey defaults.
