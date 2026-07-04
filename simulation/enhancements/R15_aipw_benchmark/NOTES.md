# R15_aipw_benchmark — survey-weighted AIPW external benchmark (spec A7)

## What / why

Implements **spec item A7 — "External survey-causal benchmark (compared to what?)"** from
`Writing/comments/phase4-arc-sim-specs.md`: all five paper arms are the SAME estimator
(survey-TMLE) design-ablated; A7 adds an **independent, literature-standard competitor** —
survey-weighted AIPW with a survey **JKn replicate** SE — as a full arm across the ladder
(L1–L4) × both designs (`standard`, `R1`), and on one NHANES example (**E1**, best overlap).
This lets the reader judge (i) whether the paper's coverage story is TMLE-specific, and
(ii) the **two-fold-splitting SE cost** (the locked Fully-Aware-CF pays SE/SD ≈ 1.44 at L4).

NON-DESTRUCTIVE: sources `codes/{config,dgp,estimators,diagnostics,learners}.R` +
`codes/arc_runs/_checkpoint.R` read-only; all new logic lives in this folder; outputs go to
`$R15_OUT` / `$R15_E1_OUT` / `results/arc/`.

## Files

| file | role |
|---|---|
| `aipw_helpers.R` | `aipw_arms()` (AIPW-SF + AIPW-CF, JKn + Eq-8 SEs), `.aipw_jkn()`, `one_rep_aipw()` |
| `run.R` | simulation array driver (mirrors `R04_nuisance_rate/run.R`) |
| `submit.slurm` | array job: **FULL `--array=1-80`** (2 scen × 4 rungs × 10 chunks of 100 reps) |
| `nhanes_e1.R` | NHANES E1 driver (single task, B-split convention) |
| `submit_e1.slurm` | single-task job for `nhanes_e1.R` |
| `aggregate.R` | integrity guards (duplicate-rep `stop` + per-cell completeness WARN) + summary + the three decision views |

## Estimator design (locked)

* **Nuisances two ways**: **AIPW-SF** = WEIGHTED single-fit `.sl()` Q on (A,W) and g on W,
  in-sample predictions (the AIPW analogue of the Fully-Aware single-fit arm);
  **AIPW-CF** = UNWEIGHTED out-of-fold fits over `make_cf_folds()` — the same PSU-within-stratum
  protocol as the engine's primary CF arm. Both: g floored to **[0.05, 0.95]**, Q truncated
  to **[1e-3, 1−1e-3]** (engine clamps).
* **Weight normalization for the weighted fits (smoke-driven fix)**: the SF nuisances are fit
  with `weights = w/mean(w)`. The weighted MLE is invariant to the weight SCALE, but RAW
  survey weights (1e1–1e5) make glm's IRLS diverge spuriously (coefs → ±1e15, fitted → {0,1};
  observed on E1 and on sim reps 3/20 with raw weights, converged & sane when normalized).
  `tmle()` applies the SAME normalization internally (`obsWeights/sum(obsWeights)*n`), so this
  makes AIPW-SF the honest analogue of the locked weighted Fully-Aware arm rather than a
  numerically sabotaged strawman. The Hajek point estimate and BOTH survey variances use the
  RAW weights. Residual "collapsed-g" reps (every raw g outside the floor) would now indicate
  GENUINE weighted quasi-separation; the aggregator counts them per cell (expected ~0).
* **Point estimate** (each arm): Hajek-weighted AIPW,
  `D_i = A_i(Y_i−Q1_i)/g_i − (1−A_i)(Y_i−Q0_i)/(1−g_i) + Q1_i − Q0_i`,
  `psi_hat = Σ w_i D_i / Σ w_i`.
* **Variance TWO ways per arm**:
  * `se_jkn` — `svydesign → as.svrepdesign(type="JKn") → withReplicates()` of the Hajek
    ratio of D with the **nuisances FROZEN** (Q, g fit once; only weights perturbed across
    replicates). Freezing the nuisances is the **standard practice** for replicate-variance
    AIPW; the replicate SE is conditional on the fitted nuisances, exactly as the paper's
    Eq-8 SE plugs the fitted EIF into the design linearization — so the two SEs estimate the
    same (first-order) target and are directly comparable.
  * `se_lin` — the engine's `.se_des()` (Eq-8) on the **centered** pseudo-outcome
    `(D − psi_hat)`: apples-to-apples with the five arms' EIF-based SE.
  * `df` = the engine's design df from `.se_des` (#PSU − #strata; domain-aware). Both
    coverages use a t(df) reference.
* **NHANES E1 domain JKn**: replicate design built on the **FULL** MEC design, then
  `subset(rep_des, inpop)` → `withReplicates` (the replicate-weight analogue of the engine's
  subset-the-DESIGN rule). If the subset path errors, falls back to the algebraically
  identical domain-zeroed full-sample ratio `Σ wts·D·1{inpop} / Σ wts·1{inpop}`; the realized
  mode is recorded as `jkn_mode` (smoke observed: `subset_repdes`).

## Comparability (preserve!)

Reps are seeded `draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i)` with **i the GLOBAL rep
index** — the exact formula of the locked headline run (`run_sim.R`) — so for each
(scenario, rep) the analysis sample is **byte-identical** to the one the five locked arms saw.
AIPW-vs-TMLE contrasts are therefore paired (estimator differences, not sampling MC noise).
Worker RNG streams (folds + SL internal CV only) use
`clusterSetRNGStream(iseed = SAMPLE_SEED_BASE + 1000*cell_index + task)` — the 1000·cell offset
keeps streams distinct across cells/chunks (task ≤ 80 < 1000, so no collisions).

## Grid / cost

FULL: 2 scenarios × 4 rungs × 1000 reps in chunks of 100 → **80 tasks**, `--array=1-80`
(order: tasks 1–10 standard×L1, 11–20 R1×L1, …, 61–70 standard×L4, 71–80 R1×L4).
Per rep = 12 `.sl` calls (SF 2 + CF 5 folds × 2) + 2 svrep builds — **cheaper than one
five-arm `run_estimators()` rep** (~16 fits + 3 targeting passes). Expected per-chunk wall on
32 cores: L1 ≈ 2–5 min, L2 ≈ 5–15 min, L3 ≈ 10–30 min (3-learner ladder w/ 10-fold internal
CV), L4 ≈ 5–15 min (deep RF, single learner, V=2 internal CV). `--time=04:00:00` is the
L4-deep-RF binding budget with ≥4× headroom.

**Walltime recovery (if a chunk hits the limit anyway).** Resubmit **ONLY the failed task
ids at the UNCHANGED `SIM_CHUNK=100`** with a longer `--time`, e.g.
`sbatch --array=63,67 --time=08:00:00 codes/arc_runs/R15_aipw_benchmark/submit.slurm`.
This is truly idempotent: the checkpoint filename is keyed on the chunk NUMBER, whose rep
range is derived from `SIM_CHUNK` at runtime, so completed chunks re-skip and failed chunks
recompute the exact same reps. **Do NOT change `SIM_CHUNK` on a partially complete run**
(e.g. halving to 50 with `--array=1-160`): the chunk→rep mapping shifts, so stale
`CHUNK=100` files alias the new chunk numbers — missing reps are never computed and
surviving reps are double-counted (`aggregate.R`'s duplicate guard will `stop`, and its
completeness check WARNs on short cells, but the run itself is wasted). If you MUST change
`SIM_CHUNK`, FIRST delete/move ALL `r15_<scenario>_<rung>_chunk*.rds` (and their manifests)
for the affected cells, then resubmit the full new array.

Rough total ≤ ~300 core-hours. E1: B=20 splits in parallel ≈ one split's wall
(primary library on ~31k rows, ~15–45 min) → `--time=02:00:00`.

## Commands

**Simulation SMOKE (run locally first; ~2–4 min)** — one cheap cell (standard × L1_param),
20 reps, 2 cores, truth_M = 2e5:

PowerShell:
```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R15_DIR="<repo-root>/simulation/enhancements/R15_aipw_benchmark"
$env:R15_OUT="<repo-root>/sim_output/arc_runs/R15_aipw_benchmark"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R15_aipw_benchmark/run.R"
```

bash (ARC login node or Linux):
```bash
cd /path/to/survey-tmle2
SMOKE=1 REPO_ROOT=$PWD SIM_CODE=$PWD/codes DATA_ROOT=$PWD/sim_output \
R15_DIR=$PWD/codes/arc_runs/R15_aipw_benchmark \
R15_OUT=$PWD/sim_output/arc_runs/R15_aipw_benchmark \
SLURM_CPUS_PER_TASK=2 Rscript codes/arc_runs/R15_aipw_benchmark/run.R
```

**NHANES E1 wiring SMOKE (~1–2 min)** — B=2 splits, `SL.glm` only:

PowerShell:
```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:R15_DIR="<repo-root>/simulation/enhancements/R15_aipw_benchmark"
$env:NH_ANA="<repo-root>/Nhanes/analytic"
$env:NHANES_B="2"; $env:R15_LIB_OVERRIDE="SL.glm"; $env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R15_aipw_benchmark/nhanes_e1.R"
```

bash:
```bash
cd /path/to/survey-tmle2
SMOKE=1 REPO_ROOT=$PWD SIM_CODE=$PWD/codes \
R15_DIR=$PWD/codes/arc_runs/R15_aipw_benchmark NH_ANA=$PWD/Nhanes/analytic \
NHANES_B=2 R15_LIB_OVERRIDE=SL.glm SLURM_CPUS_PER_TASK=2 \
Rscript codes/arc_runs/R15_aipw_benchmark/nhanes_e1.R
```

**Submit (ARC, from the scratch repo root):**
```bash
sbatch codes/arc_runs/R15_aipw_benchmark/submit.slurm        # sim, --array=1-80
sbatch codes/arc_runs/R15_aipw_benchmark/submit_e1.slurm     # NHANES E1, single task
# TEST first (optional): sbatch --array=1,61 codes/arc_runs/R15_aipw_benchmark/submit.slurm
```

**Aggregate (after the array completes; local or ARC):**
```bash
R15_OUT=$PWD/sim_output/arc_runs/R15_aipw_benchmark SIM_CODE=$PWD/codes \
Rscript codes/arc_runs/R15_aipw_benchmark/aggregate.R
```
(PowerShell: set `$env:R15_OUT`, `$env:SIM_CODE`, then run the same Rscript. The E1 CSV is
written by `nhanes_e1.R` itself — no separate aggregation step.)

## Expected outputs

* `$R15_OUT/r15_<scenario>_<rung>_chunk###.rds` × 80 + `manifest/` RDS per task.
* `results/arc/R15_aipw_benchmark_summary.csv` — long format, `se_type ∈ {jkn, lin}` leading,
  then the project-convention columns (scenario, rung, method, n_reps, Psi, bias, emp_sd,
  mean_se, se_ratio, coverage, mcse_cov, deff_clust, icc_D).
* `$R15_OUT/R15_aipw_benchmark_combined.rds` — per-rep + summary + per-rep diagnostics.
* `results/arc/R15_aipw_nhanes_E1.csv` — per arm: b, b_split_sd, se_jkn, se_lin, df,
  lcl/ucl (JKn, t-df), lcl_lin/ucl_lin, B, library, jkn_mode (SF arm) + jkn_mode_cf
  (CF arm; recorded separately so a divergent subset-vs-fallback choice between the
  two arms would be visible).

**KEY NUMBERS to inspect (full run):**

1. **L1–L3 sanity**: AIPW-CF (and AIPW-SF at L1–L2) bias ≈ locked FA/CF bias (~0.01–0.02),
   coverage ~0.93–0.95 with BOTH SEs; `n_gsf_collapsed = 0` everywhere.
2. **JKn ≡ Eq-8 identity (sim only)**: `jkn_over_lin` must be **exactly 1** in every sim cell.
   With stratum-constant weights and equal PSU takes, deleting a PSU leaves the JKn-rescaled
   Hajek denominator unchanged, so the JKn ratio variance reduces *algebraically* to the
   linearization (verified to machine precision vs a toy design where it breaks). This is a
   built-in plumbing check, and it means the jkn-vs-lin comparison is only *informative* on
   NHANES (where weights vary within strata: E1 smoke ratio 1.0001).
3. **L4 headline ("compared to what?")**: does AIPW-SF break down like the locked single-fit
   FA (bias ~0.07–0.12, coverage ~0.18–0.22) while AIPW-CF holds near-nominal like the locked
   CF (bias ≤ ~0.01, coverage ~0.95+)? If yes, the failure-and-rescue is estimator-generic,
   not TMLE-specific.
4. **Fold-splitting SE cost at L4**: locked CF se_ratio ≈ **1.44** (standard 1.4389,
   R1 1.4403). Compare AIPW-CF se_ratio — if it is also ≈ 1.3–1.5, the SE cost is the price
   of PSU-level fold splitting per se; if ≈ 1.0, it is specific to the TMLE pooled-targeting
   construction.
5. **E1**: AIPW-SF/CF b should sit near the locked E1 arms (~0.035–0.045; the smoke with
   SL.glm gave SF 0.0436, CF 0.0448 vs locked CF 0.0363); JKn CI ≈ linearized CI
   (df = 184 − 90 = 94); `jkn_mode = subset_repdes`.

## SMOKE-GATE DECISION RULE

* **Sim smoke** prints `[SMOKE-GATE] PASS` iff: both arms present; ALL b/se_jkn/se_lin finite;
  per arm `mean(se_jkn)/mean(se_lin) ∈ [0.5, 2]`; per arm `|MEDIAN bias| ≤ 0.05` (median is
  the gate statistic so a single extreme AIPW rep cannot mask/flag the gate via the mean);
  no dropped reps; zero collapsed CF g-fits.
* **E1 smoke** prints `[SMOKE-GATE] PASS` iff: both arms finite (both SEs);
  `se_jkn/se_lin ∈ [0.5, 2]`; `jkn_mode ∈ {subset_repdes, domain_zeroed_full}`.
* **PASS → submit the FULL arrays. STOP → fix this folder's code (never the engine) and
  re-smoke.** If the full run later shows AIPW-SF *covering* at L4 (≥0.90), do not reframe
  silently — that would say the single-fit failure is TMLE-specific; flag for the authors.

**Smoke executed 2026-06-12 (R 4.5.1, Windows, 2 cores): BOTH PASS.**
* Sim (standard × L1_param, 20 reps, ~0.3 min): AIPW-SF bias 0.0236 (median 0.0177),
  emp_sd 0.0443, cov 0.95/0.95; AIPW-CF bias 0.0269 (median 0.0161), emp_sd 0.0457,
  cov 0.95/0.95; jkn_over_lin = 1 exactly (the structural identity above); df 50;
  collapsed-g reps SF=0 CF=0; deff ≈ 1.22; Psi = 0.20925.
* E1 wiring (B=2, SL.glm, ~0.7 min): AIPW-SF b 0.04355 (se_jkn 0.007896, se_lin 0.007895),
  AIPW-CF b 0.04478 (se 0.00795), df 94, jkn_mode = jkn_mode_cf = subset_repdes (the subset
  path worked in BOTH arms; no fallback needed); locked five-arm context printed alongside
  (locked CF b 0.03633).
* An earlier smoke (pre-fix) exposed the raw-weight IRLS divergence in the weighted SF fits
  (sim: 6/20 reps with every raw g outside the floor, b outliers −1.39/+3.01; E1: SF
  b = 0.442) → fixed by the mean-1 weight normalization (see Estimator design) and re-run
  clean. This is why the gate uses the median and counts collapsed-g reps.

## Caveats / assumptions

* **Frozen-nuisance replicate variance**: the JKn SE does NOT refit Q/g per replicate (that
  would be ~185× the cost and is not standard practice); it is conditional on the fitted
  nuisances, the same first-order status as the Eq-8 plug-in EIF SE. State this in any
  manuscript text that reports se_jkn.
* **Hajek AIPW has no targeting step**: its point estimate differs from the TMLE arms'
  (tmle's Qstar-based ψ); it does not respect the [0,1] outcome bounds. That is the point —
  it is the standard external competitor, not another TMLE ablation.
* AIPW-SF couples **weighted nuisances + in-sample prediction** (the FA analogue) and AIPW-CF
  couples **unweighted nuisances + cross-fitting** (the primary-CF analogue) — the same
  two-axis bundling as the paper's headline contrast; the R03 isolation run is the place
  where the axes are unbundled, not here.
* The harmonized g floor (0.05 both arms) differs from the locked FA arm's effective 1e-3
  EIF floor; AIPW-SF is therefore a *floored* single-fit competitor (the cleaner benchmark).
* In the sim, `se_jkn` carries no information beyond `se_lin` (the exact identity in KEY
  NUMBERS #2); any manuscript claim about replicate-vs-linearization agreement should cite
  the NHANES E1 numbers, where the two genuinely differ.
* `df` is taken from the Eq-8 design (`.se_des`), also used for the JKn t-CI; on these designs
  `degf(rep-design source)` is identical (sim: #PSU−#strata; E1: every PSU intersects the
  domain).
* E1 manifest records the split seeds (20260607+1:B); POP_SEED/TRUTH_SEED are sim-only and
  appear in the sim manifests.
* **Partially failed chunks still checkpoint as complete** (failed reps are dropped, not
  retried): `run.R` prints a loud `[task N] WARNING: x/y reps failed in this chunk ...`
  line before saving and records `n_failed` in the chunk RDS + manifest; to recompute,
  delete that chunk's RDS and resubmit the task id. `aggregate.R`'s completeness check
  WARNs on the resulting short cell.
* Local smoke leaves `SMOKE_*` files under `$R15_OUT`, `$R15_E1_OUT`, and
  `results/arc/SMOKE_R15_aipw_nhanes_E1.csv` — harmless, excluded by all aggregators.
