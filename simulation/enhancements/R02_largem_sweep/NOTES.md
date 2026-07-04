# R02_largem_sweep — Fixed-L4 larger-m sweep

## Purpose

Settle whether the **L4 cross-fit (Fully-Aware-CF) over-coverage** seen in the
locked results (`results/sim_full_summary.csv`: standard/L4_aggressive CF has
`se_ratio = 1.4389`, `coverage = 0.985`; R1/L4 CF `se_ratio = 1.4403`,
`coverage = 0.992`) is a **finite-m artifact** that shrinks as the number of
sampled PSUs `m` grows (supporting **Theorem 2 as conservative-consistency**),
or is **structural** (persists at large m).

The DGP, finite population, and SuperLearner library are held **FIXED** at the
`standard` scenario; the only thing that varies is the number of **sampled
PSUs per stratum** `m = base_m`, swept across:

| level | base_m | m_total = H·base_m (H=10) | df_design = H·(base_m−1) | multiple |
|------:|-------:|--------------------------:|-------------------------:|---------:|
| baseline | 6  | 60  | 50  | 1×    |
| 2×       | 12 | 120 | 110 | 2×    |
| ~3.3×    | 20 | 200 | 190 | ~3.3× |
| 5×       | 30 | 300 | 290 | 5×    |

`base_m` is `draw_sample()`'s **existing argument** (stage-1 SRS-WOR PSU count).
Increasing it raises the number of *independent* sampled PSUs — the asymptotic
dimension `m` in Theorem 2 — while holding population geometry (`H=10`,
`J_per_stratum=60`, `M_per_psu=200`) and units-per-PSU (`base_n0`, default 25)
**fixed**, so cluster size and the weight structure don't change. Constraint:
`base_m ≤ J_per_stratum = 60`, so `30` is the largest valid 5× level.

Each task runs **four rungs** at its `m`:
- `L4_aggressive` (`SL.ranger.deep`) — the over-covering deep-RF rung (primary).
- `RF_shallow` (`SL.ranger.t1`) — default-depth RF **contrast**: its `se_ratio`
  should approach 1 *faster* than deep-RF if the over-coverage is overfit-driven.
- `L1_param` (`SL.glm`) and `L3_adaptive` — **bias track** to test the O(1/m)
  Hájek-ratio explanation for the R01 smooth-rung bias (cheap; ride along).

1000 reps per (m, rung).

## What is NON-DESTRUCTIVE here

- Sources the canonical engine read-only: `codes/{config,dgp,estimators,
  diagnostics,learners}.R`. **No canonical file is edited.**
- All new logic lives in `helpers.R` (this folder): `one_rep_m()` (forwards
  `base_m` to `draw_sample`) and `summarise_sweep()` (matches `aggregate_sim.R`
  columns). The m knob is the engine's existing `base_m` arg — `R/dgp.R` is
  untouched.
- Per-task RDS → `sim_output/arc_runs/R02_largem_sweep/`
  (ARC: `$DATA_ROOT/arc_runs/R02_largem_sweep/`). Manifests in a `manifest/`
  subdir. Summary CSV → `results/R02_largem_sweep_summary.csv`. The locked
  `results/sim_full_summary.csv` and `DATA_ROOT/{intermediate,results}` are
  never touched.

## Files

- `run.R`       — array driver (one task = one m-level × one rep-chunk; all rungs).
- `helpers.R`   — `one_rep_m()` + `summarise_sweep()` (the only new logic).
- `aggregate.R` — combine per-task RDS → summary + diagnostics CSV.
- `submit.slurm`— SLURM array (1–40), 32 cpus, 64G, 06:00:00.

## How to SMOKE-TEST (do this first)

Two tiers — both write `SMOKE_*` RDS so they can't be mistaken for real task
files (the aggregator ignores `SMOKE_` files). Both build the population with a
reduced `truth_M=2e5` and cap workers at 2 cores to stay light on RAM.

**Tier 1 — pipeline validation (default; ~3-4 min, validated locally on R 4.5).**
Runs the cheap `L1_param` (GLM) rung at the largest m (base_m=30), 20 reps. This
confirms the m knob, the column schema, and the design df end-to-end WITHOUT the
heavy deep-RF cost. This is the gate the user should run first.

ARC login node (or any box with R 4.4/4.5 + SuperLearner/tmle/survey/ranger):

```bash
cd "$SLURM_SUBMIT_DIR"          # repo root / submit dir
module load gcc/9.4.0 r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/
export REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" SIM_CODE="$PWD/codes"
SMOKE=1 SLURM_ARRAY_TASK_ID=1 Rscript simulation/enhancements/R02_largem_sweep/run.R
```

Windows/dev (R 4.5):

```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:SIM_CODE="<repo-root>/codes"; $env:SLURM_ARRAY_TASK_ID="1"
Rscript <repo-root>/simulation/enhancements/R02_largem_sweep/run.R
```

Expected Tier-1 output: banner `[SMOKE] m=30, rung=L1_param, 20 reps ...`, then
`[task 1] rung=L1_param done (20 reps, ~4 min)`, then a saved
`SMOKE_sim_m30_chunk001.rds` (80 est rows = 20 reps × 4 methods — L1 is a single
learner so the Fully-Aware-CV arm does not fire) and a manifest. VALIDATED:
`df_design = 290` = H·(base_m−1) = 10·29, `n ≈ 7590`, `w_cv ≈ 1.41`, all
`b`/`se` finite, and `summarise_sweep()` emits the aggregate_sim.R column set.

**Tier 2 — heavy decision-cell eyeball (optional; minutes PER REP, only if you
want to see deep-RF at m=30 locally).** Set `SMOKE_RUNG=L4_aggressive` (or
`RF_shallow`) to run that rung at base_m=30 with 4 reps:

```bash
SMOKE=1 SMOKE_RUNG=L4_aggressive SLURM_ARRAY_TASK_ID=1 \
  Rscript simulation/enhancements/R02_largem_sweep/run.R
```

This is the cell that gets a 6h ARC walltime; locally a single deep-RF rep
exercises every arm (two weighted `tmle()` SL fits + the CF fold loop) and can
take several minutes, so 4 reps is the most you'd want on a laptop. On a 16 GB
laptop keep `cores=2` (the default smoke cap) — running 4+ concurrent deep-RF
workers can exhaust RAM and kill a PSOCK worker (`unserialize` error). ARC's 64G
× 32 cpus has no such issue.

### Smoke-gate decision rule (STOP-and-report vs proceed)

Tier 1 is a PIPELINE gate; Tier 2 (or the full run) is the SCIENCE gate.

- **STOP-and-report immediately** (do NOT submit) if Tier 1 errors, produces any
  non-finite `b`/`se`, or `df_design ≠ 290` at m=30 — that would mean the
  `base_m` knob isn't changing the design dimension as intended and the whole
  premise of the sweep is broken.
- **PROCEED to submit the full array** if Tier 1 passes (finite estimates,
  df_design=290, correct columns). The full 1000-rep run is what actually
  estimates the L4 se_ratio-vs-m curve.
- **Science decision (apply after the full run, or after a Tier-2 eyeball):** if
  the deep-RF (`L4_aggressive`) Fully-Aware-CF `se_ratio` at the largest m
  (m_total=300) is **meaningfully below the baseline 1.44** (e.g. ≲ 1.25) and
  `RF_shallow` is closer to 1 than deep-RF → finite-m / conservative-consistency
  (Theorem 2), keep the **"se_ratio → 1"** narrative. If the deep-RF CF
  `se_ratio` is **still > 1.3 at m=300** → pivot the manuscript claim to
  **"bounded valid conservatism, decreasing in m"** rather than "se_ratio → 1".
  A Tier-2 eyeball with 4 reps cannot resolve 1.30 vs 1.44 — only the 1000-rep
  run can; Tier 2 is just a direction/sanity peek.

## How to SUBMIT (full run)

From the repo root on ARC:

```bash
sbatch simulation/enhancements/R02_largem_sweep/submit.slurm
```

Array 1–40 = 4 m-levels × 10 rep-chunks (1000 reps / 100 per chunk). All four
rungs run inside every task. Walltime is 06:00:00 because deep-RF (L4) at
m_total=300 (base_m=30) fits 500 full-depth trees per fold per rep and is the
heaviest cell; 32 cores parallelize the 100-rep chunk.

## Aggregate (after the array finishes)

```bash
export REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" SIM_CODE="$PWD/codes"
Rscript simulation/enhancements/R02_largem_sweep/aggregate.R
```

Writes:
- `results/R02_largem_sweep_summary.csv` — one row per (m, rung, method)
  with columns matching `aggregate_sim.R` plus the swept knobs:
  `scenario,rung,base_m,m_total,method,n_reps,Psi,bias,emp_sd,mean_se,
  se_ratio,coverage,mcse_cov,mean_df,deff_clust,icc_eif`.
- `results/R02_largem_sweep_diag.csv` — per (m, rung) targeting/positivity
  (`eps_fa,g_fa_*,eps_cf,g_cf_*,cf_V_eff,deff,icc_eif`).
- `sim_output/arc_runs/R02_largem_sweep/R02_combined.rds` — full per-rep detail.

The aggregator also prints three decision views: (1) full summary, (2)
**Fully-Aware-CF se_ratio & coverage vs m** for L4 + RF_shallow (the headline),
(3) **Fully-Aware bias vs m at L1/L3** (the O(1/m) Hájek check).

## Expected output / what to inspect

- **Headline (settle the question):** in the CF view, plot `se_ratio` vs
  `m_total` for `L4_aggressive`. Monotone decrease toward 1 ⇒ finite-m,
  conservative-consistent (Theorem 2). A plateau well above 1 ⇒ structural →
  pivot to "bounded valid conservatism, decreasing in m". `coverage` should
  stay ≥ 0.95 throughout (conservatism is valid, never anti-conservative).
- **Contrast:** `RF_shallow` CF `se_ratio` should decay toward 1 **faster** than
  deep-RF, consistent with the overfit interpretation of the deep-RF excess.
- **Bias track:** `Fully-Aware` `bias` at `L1_param`/`L3_adaptive` should shrink
  ~∝ 1/m if the R01 smooth-rung bias is a Hájek-ratio finite-m effect.
- **Sanity:** `mean_df` ≈ H·(base_m−1) = {50,110,190,290}; `coverage` MCSE
  `mcse_cov ≈ 0.007` at 1000 reps.

## Dependencies

- R 4.4 (ARC) / 4.5 (dev); packages: `SuperLearner, tmle, survey, surveyCV,
  ranger, earth, gam, glmnet` (same set as `run_sim.R`). No new packages.
- Reuses `config.R` seeds (POP_SEED=20260606, SAMPLE_SEED_BASE=1000,
  L'Ecuyer-CMRG). Per-task RNG stream offset by `1000·m_idx + task` so m-levels
  and chunks never collide.

## Caveats

- The m knob is **sampled** PSUs per stratum, not population PSUs. The
  population (and hence the truth `Psi`) is identical across all m — only the
  design's PSU count and df change. This is exactly the Theorem-2 asymptotic
  regime (m → ∞ with fixed cluster size), which is what we want to probe.
- `base_m ≤ 60` (= `J_per_stratum`). To go beyond 5× you'd have to rebuild the
  population with more population PSUs (a different study); not done here.
- The R1 scenario is intrinsically `m_h = 2` PSUs/stratum and is **not** swept
  (its design df is structurally thin); this run is `standard` only, where the
  L4 over-coverage is cleanest to trace.
- Deep-RF at m=30 is the heavy cell; if any task hits the 6h walltime, split
  `SIM_CHUNK` to 50 and double the array (`--array=1-80`, 20 chunks/m).
