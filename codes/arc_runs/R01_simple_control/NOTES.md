# R01_simple_control — correctly-specified ('simple') L1 calibration control

## Purpose
A clean **pipeline-correctness control** that explains the +0.012..+0.021
cross-fit (and single-fit) bias seen at L1–L3 in the headline run
(`results/sim_full_summary.csv`).

The headline run builds the population with `model_type='complex'` — the
Kang–Schafer transforms in `codes/dgp.R::apply_L`. Under `complex`, the
learners see the transformed `L1..L4`, **not** the linear `C1..C4` that actually
drive `A` and `Y`, so even L1's GLM is misspecified and biased. This run re-runs
the **same engine and same five arms** on `model_type='simple'` (where `apply_L`
is the identity, so the parametric/smooth learners are *correctly specified* and
bias should be ~0), for **both scenarios** (`standard`, `R1`), **1000 reps**.
It ALSO re-runs L1 under `complex` so the matched complex/simple pair sits side
by side. Every output row carries a **`model_type` column** (the canonical
tables lack it).

Cells (`run.R` grid = 2 scenarios × 3 model_type/rung pairs = 6 cells):

| model_type | rung      | expectation                                   |
|------------|-----------|-----------------------------------------------|
| simple     | L1_param  | bias ~0, coverage ~0.94–0.95 (the control)    |
| simple     | L2_smooth | context (also correctly specified)            |
| complex    | L1_param  | reproduces the headline L1 bias (side by side)|

All five arms run per cell (Fully-Aware, Fully-Aware-CF, Partially-Aware,
Non-Aware; Fully-Aware-CV only emits for multi-learner libraries, i.e. L2).

## Non-destructive guarantees
- `run.R` **sources** `codes/{config,dgp,estimators,diagnostics,learners}.R`
  unchanged. No canonical file is edited.
- The only new logic is (a) looping `model_type` and building the population
  with the **matching** `model_type` (the headline driver hardcodes `'complex'`),
  and (b) drawing samples with that `model_type`. Both reuse `make_population` /
  `draw_sample` / `run_estimators` / `deff_clust` as-is.
- Outputs land in a **private** subtree and never clobber locked results:
  - per-task RDS + manifest → `sim_output/arc_runs/R01_simple_control/`
    (ARC) / `<DATA_ROOT>/arc_runs/R01_simple_control/` (local, via `ARC_OUT`).
  - summary CSV → `results/arc/R01_simple_control_summary.csv`.

## Files
- `run.R` — the driver (sources the engine; loops model_type+rung; tags output).
- `aggregate.R` — combines per-task RDS → `results/arc/R01_simple_control_summary.csv`
  (same coverage/bias/se_ratio/MCSE math as `aggregate_sim.R`, plus a
  `model_type` column and a printed simple/L1 PASS/FLAG verdict).
- `submit.slurm` — SLURM array job (mirrors `codes/sim.slurm` resources).

## Smoke-test (local, before submitting) — ~1–3 min
Runs ONE cell (standard, simple, L1), 20 reps, single chunk, ≤4 cores.

PowerShell:
```powershell
$env:SMOKE="1"; $env:SIM_CODE="<repo-root>/codes"
Rscript <repo-root>/codes/arc_runs/R01_simple_control/run.R
Remove-Item Env:SMOKE
```
bash / ARC login node:
```bash
cd /path/to/survey-tmle2
SMOKE=1 SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  Rscript codes/arc_runs/R01_simple_control/run.R
```
This writes one RDS to `<DATA_ROOT>/arc_runs/R01_simple_control/`
(`r01_standard_simple_L1_param_chunk001.rds`). To inspect the smoke result you
can aggregate the single file:
```bash
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  Rscript codes/arc_runs/R01_simple_control/aggregate.R
```

Requirements: R 4.4/4.5 with the same packages the engine uses
(`SuperLearner`, `tmle`, `survey`, `surveyCV`, `earth`, `gam`, `ranger`).

## Submit (ARC, full 1000-rep run)
From the repo root (= `SLURM_SUBMIT_DIR`):
```bash
sbatch codes/arc_runs/R01_simple_control/submit.slurm     # --array=1-60 (FULL)
```
For a quick on-cluster test first, edit the three TEST/FULL knobs in
`submit.slurm` (`--array=1-6`, `SIM_N_REPS=100`) and submit.

## Aggregate (after the array completes)
```bash
SIM_CODE="$SLURM_SUBMIT_DIR/codes" REPO_ROOT="$SLURM_SUBMIT_DIR" \
  DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output" \
  Rscript codes/arc_runs/R01_simple_control/aggregate.R
```
Pull `results/arc/R01_simple_control_summary.csv` back to the repo.

## Expected output
- Per-task: `sim_output/arc_runs/R01_simple_control/r01_<scenario>_<model_type>_<rung>_chunk###.rds`
  plus a matching `manifest_*.rds` (params, truth, pop ICC audit, seeds, package
  versions, git rev) — 60 of each for the full run.
- Summary: `results/arc/R01_simple_control_summary.csv` with columns
  **`scenario, model_type, rung, method, n_reps, Psi, bias, emp_sd, mean_se,
  se_ratio, coverage, mcse_cov, deff_clust, icc_eif`** — identical to
  `sim_full_summary.csv` plus the leading `model_type` column.

Key numbers to inspect: the `simple`/`L1_param`/**Fully-Aware** rows for each
scenario (`bias` ≈ 0, `coverage` ≈ 0.94–0.95), compared against the matching
`complex`/`L1_param` rows (which should reproduce the headline +0.012..+0.021
bias and sub-nominal-but-close coverage).

## Smoke-gate decision rule (STOP-and-report vs proceed)
The control's whole point is that, under a *correctly specified* DGP, the
pipeline is unbiased and nominal. So:

- **PASS / proceed**: under `model_type='simple'`, the L1 **Fully-Aware** rows
  show `coverage` ≈ 0.94–0.95 and `|bias|` within ~2×MCSE of 0
  (MCSE-of-the-mean = `emp_sd/sqrt(n_reps)`). This confirms the cross-fit/
  single-fit machinery is correct and the headline L1–L3 bias is *attributable
  to misspecification under `complex`*, not a pipeline defect — green light to
  write the Theorem 1/2 prose.
- **STOP-and-report (privately, before writing Theorem 1/2)**: if any
  `simple`/L1 Fully-Aware row has `coverage` outside ~0.93–0.96 OR `|bias|`
  exceeding ~2×MCSE of 0. That would be a **real pipeline finding** (the engine
  is biased even when the model is correct) and must be diagnosed before any
  theory claims. `aggregate.R` prints an explicit `OK`/`FLAG` verdict per
  scenario for exactly this gate.

Note for the smoke run itself: with only 20 reps the MCSE is large, so the smoke
is a *plumbing* check (does it run end-to-end and write the tagged CSV?), not the
statistical gate. The gate is evaluated on the full 1000-rep output.

## Dependencies & caveats
- Depends on the canonical engine in `codes/` (sourced, not copied). If those
  signatures change, re-check `run.R`'s calls to `make_population`,
  `draw_sample`, `run_estimators`, `deff_clust`.
- `make_population(..., truth_M=2e6L)` is recomputed once per task (deterministic
  from `POP_SEED`); the `simple` and `complex` populations differ only through
  `apply_L` at sample-draw time, but truth is model_type-specific because
  `make_population` carries `model_type` in `params` — the `simple` and `complex`
  Psi values are computed on the SAME latent C/u law (apply_L does not touch the
  outcome/treatment models), so the two Psi columns should match within MC error.
- Seeds reuse `config.R` (`POP_SEED=20260606L`, `SAMPLE_SEED_BASE=1000L`); the
  per-task RNG stream is keyed by `SAMPLE_SEED_BASE + task`, same convention as
  `run_sim.R`, so chunks/cells never collide.
- This run does **not** include `complex`/L2 (the headline already has it) or
  L3/L4 — it is the cheap matched control, not a re-do of the full ladder.
