# R03_isolation_2x2 — Isolating cross-fitting from de-weighting at a harmonized 0.05 floor

## Purpose
The locked simulation's **Fully-Aware vs Fully-Aware-CF** contrast confounds *three*
things at once. Fully-Aware is **single-fit + weighted-nuisance + EIF floor 1e-3**;
Fully-Aware-CF is **cross-fit + UNWEIGHTED-out-of-fold-nuisance + OOF/EIF floor 0.05**.
So when CF rescues L4 coverage (0.17 -> 0.99 in `results/sim_full_summary.csv`), we
cannot tell whether **cross-fitting** or **de-weighting the nuisances** is the active
ingredient — and the appendix is about to justify cross-fitting.

This run breaks the confound by crossing the two factors and **harmonizing the floor**.
Four arms, standard scenario, at the SAME propensity floor `g_floor = 0.05` applied in
BOTH the OOF/in-sample propensity AND the EIF reconstruction:

|              | weighted-OOF nuisance | unweighted-OOF nuisance        |
|--------------|-----------------------|--------------------------------|
| **single-fit** | `SF-W` (FA-like)    | `SF-U` (**the dangerous tail**) |
| **cross-fit**  | `CF-W`              | `CF-U` (engine FA-CF-like)      |

All four hand pre-computed `(Q, g1W)` to `tmle()` for the **same weighted pooled
targeting** the engine's CV/CF arms use; the ONLY differences per cell are
(a) in-sample vs out-of-fold nuisance fits and (b) `obsWeights` to `.sl()`.

Rungs: **L4_aggressive** (`SL.ranger.deep`, the demonstration) and **L1_param**
(`SL.glm`, a placebo — GLM is Donsker, so all four arms should cover).

## Files
- `estimators_isolation.R` — the helper. Defines `run_isolation(obs, learners, g_floor=0.05, V_cf=5L)`
  returning `list(results=data.frame(method,b,se,df), diagnostics=list(eif_sfw, drow))`.
  Reuses the engine building blocks exposed by sourcing `R/estimators.R`:
  `.sl`, `.eif_from_tmle` (called with `gbound=g_floor`), `.se_des`, `make_cf_folds`.
  **Edits no canonical file.**
- `run.R` — SLURM-array driver. Sources the canonical engine read-only + this helper.
  Builds the population ONCE per task (`make_population("standard","complex")`), seeds reps
  by `SAMPLE_SEED_BASE + global_rep_index` (identical to `run_sim.R`), runs the 4 arms,
  writes one per-task RDS + manifest. Has a `SMOKE=1` path.
- `aggregate.R` — combines per-task RDS into `results/R03_isolation_2x2_summary.csv`
  using the SAME coverage/se_ratio/MCSE math and column names as `simulation/aggregate.R`.
- `submit.slurm` — array job (mirrors the project template; 32 cpus / 64G / 04:00:00).

## Outputs (non-destructive)
- Per-task RDS + manifest: `sim_output/arc_runs/R03_isolation_2x2/R03_standard_<rung>_chunk###.rds`
  (on ARC; locally `DATA_ROOT` resolves the same way).
- Summary CSV: `results/R03_isolation_2x2_summary.csv` (+ `_combined.rds`).
- **Never** touches `results/sim_full_summary.csv` or `DATA_INTERMEDIATE`.

Summary columns match the engine:
`scenario,rung,method,n_reps,Psi,bias,emp_sd,mean_se,se_ratio,coverage,mcse_cov,deff_clust,icc_eif`
where `bias=mean(b)-Psi`, `emp_sd=sd(b)`, `se_ratio=mean(se)/emp_sd`,
`coverage=mean(|b-Psi| <= t.975(df)*se)`, `mcse_cov=sqrt(cov*(1-cov)/n_reps)`.

## Smoke-test (local, ~2-4 min — DO THIS FIRST)
The smoke runs **50 reps of L4_aggressive only** (one cell), then prints a per-arm
summary AND the smoke-gate verdict. `SMOKE=1` overrides reps/rung internally, so you
only need the paths + a small core count.

PowerShell (native Windows R; recommended on this dev box):
```powershell
$env:REPO_ROOT="<repo-root>"; $env:SIM_CODE="<repo-root>/codes"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R03_DIR="<repo-root>/simulation/enhancements/R03_isolation_2x2"
$env:R03_OUT="<repo-root>/sim_output/arc_runs/R03_isolation_2x2"
$env:SMOKE="1"; $env:SLURM_CPUS_PER_TASK="6"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "$env:R03_DIR/run.R"
```
Bash/Linux equivalent:
```bash
cd /path/to/survey-tmle2
SMOKE=1 SIM_CODE=$PWD/codes REPO_ROOT=$PWD DATA_ROOT=$PWD/sim_output \
  R03_DIR=$PWD/simulation/enhancements/R03_isolation_2x2 \
  R03_OUT=$PWD/sim_output/arc_runs/R03_isolation_2x2 \
  SLURM_CPUS_PER_TASK=6 Rscript simulation/enhancements/R03_isolation_2x2/run.R
```
Writes `..._chunk001_SMOKE.rds`. Assumes R 4.4/4.5 + the same packages as the engine
(SuperLearner, tmle, survey, surveyCV, ranger).

### Smoke-gate decision rule (PRE-STATED)
Expected: at L4, **SF-U and SF-W UNDER-cover** (single-fit on a non-Donsker learner),
while **CF-W and CF-U RECOVER** ~nominal coverage. That would prove cross-fitting is
the active ingredient.

- **PROCEED** to the full run if `SF-U` coverage at L4 is **well below nominal**
  (the script flags PROCEED when `SF-U coverage < 0.90`).
- **STOP-AND-REPORT** if `SF-U` (single-fit / UNWEIGHTED) **already covers ~nominal**
  (`>= 0.90`) at L4. That is the dangerous tail: it would mean **de-weighting**, not
  cross-fitting, rescues coverage — so do NOT proceed and do NOT write the appendix's
  de-weighting / cross-fitting justification before discussing with the team.
  (50 reps -> MCSE on a 0.95 coverage ~ 0.03, so treat a borderline 0.88-0.92 SF-U as
  "inconclusive, bump smoke reps before deciding".)

## Submit the full run (ARC)
```bash
cd "$SLURM_SUBMIT_DIR"          # the repo root on ARC
sbatch simulation/enhancements/R03_isolation_2x2/submit.slurm
```
Array `1-20` = standard x {L4_aggressive, L1_param} x 10 chunks of 100 reps
(1000 reps/rung). Tasks 1-10 = L4_aggressive; 11-20 = L1_param. Then aggregate:
```bash
SIM_CODE=$PWD/codes REPO_ROOT=$PWD DATA_ROOT=$PWD/sim_output \
  R03_OUT=$PWD/sim_output/arc_runs/R03_isolation_2x2 \
  Rscript simulation/enhancements/R03_isolation_2x2/aggregate.R
```
`aggregate.R` reprints the L4 decision verdict on the full 1000-rep result.

## Expected output / what to inspect
`results/R03_isolation_2x2_summary.csv`, 8 rows (4 arms x 2 rungs). Key cells:
- **L4_aggressive**: `coverage` for SF-W ~ SF-U **low** (think ~0.15-0.55), CF-W ~ CF-U
  **~0.95**; `se_ratio` < 1 for the SF arms (anti-conservative), ~1 for CF arms.
- **L1_param** (placebo): all four arms `coverage` ~ 0.93-0.95, `se_ratio` ~ 1 — confirms
  the 2x2 machinery is unbiased when the learner is Donsker (no spurious arm effect).
- `bias` ~ small for every arm at both rungs (TMLE double-robustness); the story is in
  `coverage`/`se_ratio`, not `bias`.

## Dependencies
- Canonical engine (read-only): `codes/{config,dgp,estimators,diagnostics,learners}.R`.
- R 4.4 (ARC) / 4.5 (dev) + SuperLearner, tmle, survey, surveyCV, ranger, gam, earth, glmnet.
- No new population/DGP override needed: this run reuses the locked standard-scenario DGP
  exactly. The only *new* logic is the 2x2 nuisance-fitting in `estimators_isolation.R`.

## Caveats
- The harmonized floor (0.05 everywhere) means these SF arms are NOT byte-identical to the
  locked Fully-Aware (which floors the EIF at 1e-3). That is the point — we are removing the
  floor asymmetry so the only thing that varies is cross-fit x weighting. SF-W will therefore
  differ slightly from the locked Fully-Aware numbers; do not expect an exact match.
- Cross-fit folds are randomized within the parallel RNG stream (`make_cf_folds` calls
  `sample`); seeding follows `run_sim.R` (`clusterSetRNGStream(iseed=SAMPLE_SEED_BASE+task)`),
  so the full run is reproducible but per-rep fold assignments differ from the locked CF arm.
- `df` is the clustered design df (~`#PSU - #strata`), used for the t critical value, exactly
  as the engine does.
- If `SF-U` coverage lands borderline (0.88-0.92) at 50 smoke reps, increase smoke reps before
  applying the gate — 50 reps has MCSE ~0.03 on coverage.
