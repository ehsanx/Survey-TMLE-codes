# R12_highdeff — High-DEFF "Design C" (stress the clustering the method is sold for)

## Purpose

The headline simulation (`standard` scenario) realises an EIF clustering
design effect (`deff_clust`, computed on the **Fully-Aware EIF**) of only
**~1.25–1.40**. Reviewers asked for a regime where `deff_clust` reaches
**~2.5–4** so the clustering-aware SE that the method is sold for is genuinely
stressed.

**Design C** re-parameterises the canonical `standard` scenario by raising the
two PSU random-effect variances that drive EIF clustering and enlarging the
within-PSU (stage-2) sample size:

| knob        | canonical `standard` | Design C | role |
|-------------|----------------------|----------|------|
| `sigma2_Y`  | 1.5                  | **4.5**  | OUTCOME-logit PSU random effect (drives EIF clustering) |
| `sigma2_A`  | 0.8                  | **2.2**  | TREATMENT-logit PSU random effect (clusters the clever covariate) |
| `base_n0`   | 25                   | **55**   | units sampled per PSU at stage 2 (more correlated units / PSU) |

`sigma2_C`, `beta_strat`, `alpha_g`, `p_treat_target`, `gamma0`, `gamma_C`,
`te_log_odds` are **unchanged** from `standard`, so the estimand and overlap
regime are the same — only the clustering intensity is dialled up.

**Expected story:** Fully-Aware and Fully-Aware-CF stay ~nominal (clustering in
the SE), while **Partially-Aware** (no clustering in the SE) degrades *more*
than at low DEFF — sharpening the "clustering matters" point. Non-Aware is worst
(unweighted + iid).

**Arms reported:** Partially-Aware vs Fully-Aware vs Fully-Aware-CF (plus
Non-Aware, and Fully-Aware-CV on the multi-learner L3 rung) at three rungs:
`L1_param`, `L3_adaptive`, `L4_aggressive`.

## Non-destructive design

- The driver `run.R` **sources the canonical engine read-only**
  (`R/config.R`, `dgp.R`, `estimators.R`, `diagnostics.R`, `learners.R`).
- All Design C behaviour is in **`design_c.R`** (this folder), which only
  **passes overridden arguments** to the canonical `make_population()` and
  `draw_sample()` — both already expose `sigma2_Y`, `sigma2_A`, `base_n0`,
  `base_m`, `alpha_strat`. **No new scenario is added to `R/dgp.R`; no
  canonical file is edited.**
- Outputs go to **private** locations only; the locked headline results
  (`results/sim_full_summary.csv`, `sim_output/intermediate/...`) are never
  touched:
  - per-task RDS  → `sim_output/arc_runs/R12_highdeff/`
  - summary CSV   → `results/R12_highdeff_summary.csv`
  - combined RDS  → `sim_output/results/arc_runs/R12_highdeff/r12_combined.rds`
  - manifests     → `sim_output/manifest/arc_runs/R12_highdeff/`

## Files

| file          | what it is |
|---------------|------------|
| `design_c.R`  | Design C wrapper: `make_population_designC()`, `draw_sample_designC()`, `designC_knobs()`. Every knob is env-overridable (`DC_SIGMA2_Y`, `DC_SIGMA2_A`, `DC_SIGMA2_C`, `DC_BASE_N0`, `DC_BASE_M`, `DC_ALPHA_STRAT`) for the retune loop. |
| `run.R`       | SLURM-array driver. One task = one (rung × rep-chunk). Builds the Design C population once (`POP_SEED`), runs the chunk across cores, writes per-task RDS + per-task summary + manifest, and **prints the realised `deff_clust`**. Has a `SMOKE` mode. |
| `aggregate.R` | combines per-task RDS → `results/R12_highdeff_summary.csv` (columns match `aggregate_sim.R`). |
| `submit.slurm`| ARC array job (mirrors `codes/sim.slurm` resources; array/time/mem adjusted). |

## Smoke test (run this FIRST, locally — 1–3 min)

The smoke mode runs **L4_aggressive only, 100 reps, one chunk** — exactly the
gate cell the spec calls for. (You can shrink reps with `SIM_N_REPS` for a faster
plumbing check, e.g. 8 reps in ~1 min.)

```bash
# from the repo root <repo-root>
SMOKE=1 SIM_N_REPS=100 SLURM_CPUS_PER_TASK=6 \
  SIM_CODE="<repo-root>/codes" \
  REPO_ROOT="<repo-root>" \
  DATA_ROOT="sim_output" \
  Rscript "<repo-root>/simulation/enhancements/R12_highdeff/run.R"
```

PowerShell equivalent:

```powershell
$env:SMOKE=1; $env:SIM_N_REPS=100; $env:SLURM_CPUS_PER_TASK=6
$env:SIM_CODE="<repo-root>/codes"; $env:REPO_ROOT="<repo-root>"
$env:DATA_ROOT="sim_output"
Rscript "<repo-root>/simulation/enhancements/R12_highdeff/run.R"
```

The smoke run prints a line like:

```
[R12_highdeff task 1] REALIZED deff_clust (FA eif): mean=2.7x  median=... range=[...]  | icc_eif mean=...
```

and writes `results/R12_highdeff_summary_SMOKE.csv` plus
`sim_output/arc_runs/R12_highdeff/r12_DesignC_L4_aggressive_chunk001_SMOKE.rds`.

### Smoke-gate decision rule (STOP-and-report vs proceed)

1. **Realised DEFF.** Look at the printed `REALIZED deff_clust (FA eif): mean`.
   - **mean ≥ ~2.5** → target met → **proceed** to the full submit.
   - **mean < ~2.5** → **STOP. Do NOT submit the full run.** Retune by raising
     the variances and resubmit the smoke, e.g.
     `export DC_SIGMA2_Y=5.0 DC_SIGMA2_A=2.5 DC_BASE_N0=60` (these probed at
     ~2.9 locally). Local GLM-probe ladder for reference:
     `(Y=4.5,A=2.2,n0=55)→2.59`, `(Y=5.0,A=2.5,n0=60)→2.89`,
     `(Y=4.0,A=2.0,n0=70)→2.67`. Re-run the smoke until the mean clears ~2.5,
     then bake the winning values into `submit.slurm` (uncomment the `DC_*`
     exports) before the full run.
   - **mean > ~5** (DEFF too high, sample may be unrealistically clustered) →
     report and consider dialling back so the regime stays plausible (2.5–4).

2. **CF honesty.** If `Fully-Aware-CF` shows `se_ratio > 1.5` (over-conservative
   under the heavy clustering), this is **NOT a STOP** — frame it in the text as
   honest conservatism (CF trades a little efficiency for valid coverage under
   the stress regime). Report the number; do not retune for it.

3. **Story sanity (expected, not a hard gate at 100 reps).** Fully-Aware /
   Fully-Aware-CF coverage should be near nominal; Partially-Aware coverage
   should be visibly *lower* than its low-DEFF (`standard`) value. If
   Partially-Aware does NOT degrade relative to the headline, the run does not
   make its point — report and reconsider before spending full-run compute.

4. **Design checks.** `sumw_over_N` ≈ 1 and `min_psu_str` ≥ 2 (printed in the
   per-rep check columns / manifest). If either fails, the design is broken —
   STOP and report.

## Submit (full run on ARC)

From the SLURM submit dir (scratch copy of the repo root):

```bash
sbatch simulation/enhancements/R12_highdeff/submit.slurm
```

- Full: `SIM_N_REPS=1000`, `SIM_CHUNK=100` → 3 rungs × 10 chunks = **30 tasks**
  (`--array=1-30`, already set).
- For **500 reps** (the lighter option in the spec): set `SIM_N_REPS=500` and
  `--array=1-15` in `submit.slurm`.
- Rung/chunk mapping (run.R grid order): tasks 1–10 = `L1_param`,
  11–20 = `L3_adaptive`, 21–30 = `L4_aggressive` (chunk = ((task-1) %% 10) + 1).

## Aggregate (after the array finishes)

```bash
SIM_CODE="$PWD/codes" SIM_OUT="$PWD/sim_output/arc_runs/R12_highdeff" \
  REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  Rscript simulation/enhancements/R12_highdeff/aggregate.R
```

Writes `results/R12_highdeff_summary.csv` (one row per rung × method,
columns: `scenario,rung,method,n_reps,Psi,bias,emp_sd,mean_se,se_ratio,coverage,
mcse_cov,deff_clust,icc_eif` — identical to `results/sim_full_summary.csv`) and
prints the **realised `deff_clust` per rung** (the headline number reviewers
asked for). Smoke RDS files are excluded unless `INCLUDE_SMOKE=1`.

## Expected output / what to inspect

- **`deff_clust` column ≈ 2.5–4** across rungs (the whole point — verify it
  actually landed there; it is computed on the realised Fully-Aware EIF, the
  same gate as the headline sim).
- **Fully-Aware** and **Fully-Aware-CF** `coverage` ~0.93–0.96.
- **Partially-Aware** `coverage` materially *below* its `standard`-scenario
  value (degrades more under high DEFF) — this is the sharpened "clustering
  matters" result.
- **Non-Aware** worst (unweighted + iid SE).
- `se_ratio` for Partially-Aware < 1 (SE too small because it ignores
  clustering); for Fully-Aware/CF ~1 (possibly >1 for CF — frame as honest
  conservatism, see gate rule 2).

## Resources

- `--cpus-per-task=32`, `--mem=96G`, `--time=05:00:00`. L4 deep-RF is the heavy
  rung and Design C draws larger samples (`base_n0` 55 vs 25, n≈3.3k vs ≈1.5k),
  so memory and walltime are raised vs the headline `codes/sim.slurm` (64G / 1h).
  After the first run, tighten with `seff <jobid>`.

## Dependencies

- Canonical engine: `codes/{config,dgp,estimators,diagnostics,learners}.R`
  (sourced read-only).
- R packages (same as the headline sim): `SuperLearner, tmle, survey, surveyCV,
  earth, gam, glmnet, ranger, parallel`. On ARC these come from
  `R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/` (R 4.4.0). Locally
  validated under R 4.5.1.
- Seeds reused from `config.R`: `POP_SEED=20260606`, `SAMPLE_SEED_BASE=1000`,
  `RNGkind("L'Ecuyer-CMRG")`. Reps seeded `SAMPLE_SEED_BASE + global rep index`;
  per-task streams via `clusterSetRNGStream(SAMPLE_SEED_BASE + task)` — matches
  `run_sim.R`, fully reproducible.

## Caveats / cuttability

- **This is flagged as the most cuttable run.** If compute is tight, the *text
  argument plus the realised-DEFF report* may suffice on their own: the smoke
  run alone (L4, 100 reps) already produces a defensible
  "deff_clust reached ~2.5–4 and Partially-Aware degraded" statement. The full
  30-task array is the belt-and-suspenders version (all three rungs, 1000 reps,
  publication-grade Monte-Carlo error).
- The tuned defaults were chosen with a **fast GLM-TMLE EIF probe**; the SL deep-RF
  EIF (L4) typically realises a slightly *higher* DEFF, so ~2.5 is a conservative
  floor. **Always confirm with the smoke run** before the full submit.
- Design C is a re-parameterised `standard` scenario, **not `R1`** — the `R1`
  geometry (2 PSUs/stratum) constrains the CF fold split differently and is a
  separate axis; mixing it in would confound "more clustering" with "fewer PSUs".
- All summary columns match `aggregate_sim.R`/`run_sim.R` so the existing
  figure/table code can consume `results/R12_highdeff_summary.csv` directly.
