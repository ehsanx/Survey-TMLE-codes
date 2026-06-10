# R10_fpc — FPC / first-stage-fraction sensitivity (estimand/regime defense)

Self-contained, **non-destructive** ARC enhancement run. Sources the canonical
engine in `codes/*.R` and reuses its building blocks; **edits nothing** in
`codes/*.R`. All new logic lives in this folder (`fpc_helpers.R`). Outputs land in
`sim_output/arc_runs/R10_fpc/` and `results/arc/` — the locked results
(`results/sim_full_summary.csv`, the `intermediate/` tree) are never touched.

## Purpose

Defend two related design/estimand assumptions raised in review:

1. **The `m/M_pop -> 0` (negligible first-stage fraction) assumption**, and
2. **Design B / R1's realized first-stage sampling fraction `f ≈ 0.25`**
   (R1 samples `base_m = 2` of `J_per_stratum = 8` population PSUs per stratum).

Two analyses:

### PART A — finite-population-correction (FPC) branch
On the **R1 design only**, and **only for rungs `L1_param` and `L2_smooth`**,
recompute the **Fully-Aware-CF** design SE **with** an FPC
(`fpc = J_per_stratum`, the population PSU count per stratum) and compare coverage
**with vs without** FPC. The Fully-Aware (single weighted fit) arm is carried
through both SE formulas too, for context.

**We deliberately do NOT apply the FPC at `L4_aggressive`.** L4-CF already sits at
≈0.985 coverage (over-coverage). The FPC shrinks the SE, which at L4 would tip the
arm into **under-coverage** — the opposite of a robustness argument. The correct
place to demonstrate the FPC is where coverage is near-nominal (L1/L2), where a
modest SE narrowing should keep coverage ≥0.94. This is stated in `run.R` and
`fpc_helpers.R` headers.

The FPC SE wrapper `.se_des_fpc()` mirrors the canonical
`codes/estimators.R::.se_des(..., clustered = TRUE)` exactly — same Hájek
`svymean` of the same EIF — except the stage-1 `svydesign` carries `fpc = ~fpc1`.
With a single-stage `ids = ~cluster` design, `survey` applies the stage-1 FPC
`(1 - m_h / M_h)` to the between-PSU variance.

### PART B — first-stage-fraction sweep (`f -> 0` insensitivity)
Grow the finite population's PSUs-per-stratum `J_per_stratum ∈ {8, 16, 32, 64}`
while holding `base_m = 2` fixed (R1's default), so the realized first-stage
fraction `f = base_m / J = 2/J ∈ {0.25, 0.125, 0.0625, 0.03125}` shrinks toward 0.
`M_per_psu` is held at the R1 default (300), so the population `N` grows with `J`.
For each `J` we run the **canonical (no-FPC)** Fully-Aware-CF arm (`L1_param`) and
check that **super-population coverage stays flat as `f -> 0`** — i.e. the
`m/M_pop -> 0` assumption is innocuous because the super-population estimand and
its EIF do not depend on the fraction. We grow `J` (not shrink `m`) so that `m_h`
stays ≥2 and the design df / CF fold split are unaffected (R1 is intrinsically
`m_h = 2`).

## Files

| file | role |
|------|------|
| `run.R`          | driver. Sources `codes/*.R` + `fpc_helpers.R`; selects PART A/B by SLURM array index. |
| `fpc_helpers.R`  | **all new logic**: `.se_des_fpc()` (FPC SE), `cf_eif_for_rep()` (re-derives the CF arm EIF so both SE formulas hit the SAME EIF), `summarize_arm()`. |
| `aggregate.R`    | combines per-task RDS into the two final CSVs. Run after the array finishes. |
| `submit.slurm`   | ARC array job (mirrors `codes/sim.slurm` resources; array `1-5`). |
| `_selftest.R`    | dev-only single-core sanity check of the new logic on a tiny pop. Not part of the run; safe to delete. |

## Task grid (built in `run.R`)

```
task 1      PART A : R1, rungs {L1_param, L2_smooth}, fpc vs no-fpc      (engine x2/rep + local CF refit)
task 2      PART B : J=8   (f=0.250)  L1_param CF, no-fpc
task 3      PART B : J=16  (f=0.125)  L1_param CF, no-fpc
task 4      PART B : J=32  (f=0.0625) L1_param CF, no-fpc
task 5      PART B : J=64  (f=0.03125) L1_param CF, no-fpc   (N = 50*64*300 = 960k)
```

`FPC_N_REPS = 1000` per cell.

## How to SMOKE-TEST (locally, before submitting)

Tiny version: PART A only, **20 reps**, R1 / {L1_param, L2_smooth}, fpc vs no-fpc,
4 cores. Finishes in a few minutes on a quiet box.

```bash
# from the repo root (Linux/ARC login node or Windows Git-bash)
SMOKE=1 \
  SIM_CODE="$PWD/codes" \
  FPC_RUN_DIR="$PWD/codes/arc_runs/R10_fpc" \
  SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=4 \
  Rscript codes/arc_runs/R10_fpc/run.R
```

PowerShell:

```powershell
$env:SMOKE="1"; $env:SIM_CODE="$PWD\codes"
$env:FPC_RUN_DIR="$PWD\codes\arc_runs\R10_fpc"
$env:SLURM_ARRAY_TASK_ID="1"; $env:SLURM_CPUS_PER_TASK="4"
Rscript codes\arc_runs\R10_fpc\run.R
```

Smoke writes `*_partA_SMOKE.rds` and `results/arc/R10_fpc_partA_summary_SMOKE.csv`
and prints the FA/CF × fpc/no-fpc summary table to stdout.

Optional even-faster dev check of ONLY the new functions (tiny pop, single core):
`Rscript codes/arc_runs/R10_fpc/_selftest.R` (prints `SELFTEST OK`).

### Smoke-gate decision rule (STOP-and-report vs proceed)

After the smoke prints the PART A summary, **proceed to the full submit only if
all of these hold**:

1. The run completes and writes the RDS + CSV with no R error.
2. For both L1_param and L2_smooth, the **CF `fpc=yes` SE is ≤ the `fpc=no` SE**
   (mean_se column) and the shrink is **modest** (ratio roughly 0.90–1.00 — at
   `f=0.25` the FPC factor on the between-PSU variance is `1-0.25=0.75`, so the SE
   shrinks by at most ~`sqrt(0.75)≈0.87` if the variance were purely first-stage;
   in practice less). A ratio < 0.80 or > 1.00 means the `fpc=` wiring is wrong.
3. The local no-fpc CF SE roughly tracks the engine's CF SE
   (`se_cf_engine_check` column ≈ the `Fully-Aware-CF / fpc=no` SE; small
   differences are OK — fold RNG differs between the engine and the local refit).
4. Coverage is in a sane range (not 0, not 1, ~0.85–0.99 at 20 reps with wide
   MCSE). At 20 reps coverage is noisy, so this is only a smell test.

**STOP-and-report instead of submitting if:** the smoke errors; the FPC SE
**exceeds** the no-FPC SE (fpc misapplied); the FPC ratio is implausibly small
(< 0.80, suggesting `fpc` is being read as a 2nd-stage correction or a per-row
population that double-counts); or `se_cf_engine_check` is wildly different from
the local CF SE (CF replay diverged from the engine). In any of these cases the
finite-population claim would be mis-stated — report the discrepancy rather than
burning the array job.

## How to SUBMIT (ARC)

From the repo root on ARC (== `$SLURM_SUBMIT_DIR` == scratch):

```bash
sbatch codes/arc_runs/R10_fpc/submit.slurm
```

5 array tasks, 32 cpus / 64G each, 2 h wall (mirrors `codes/sim.slurm`; time
bumped because PART A runs the engine twice per rep + a CF refit, and PART B
grows `N` to 960k at `J=64`). Logs in `sim_output/logs/R10_fpc_%A_%a.{out,err}`.

## Aggregate (after the array finishes)

```bash
Rscript codes/arc_runs/R10_fpc/aggregate.R
```

Writes:
- `results/arc/R10_fpc_partA_summary.csv` — FA & CF × {fpc, no-fpc} for L1/L2 on R1.
- `results/arc/R10_fpc_fraction_sweep.csv` — CF coverage vs `J`/`f`/`N_pop`.

## Expected output

**Per-task RDS** in `sim_output/arc_runs/R10_fpc/`:
`R10_fpc_partA.rds`, `R10_fpc_partB_J008.rds … _J064.rds`. Each holds per-rep rows,
`Psi`, `truth`, `params`.

**Summary CSVs** in `results/arc/` (columns follow the run_sim/aggregate
convention — `method, fpc, n_reps, Psi, bias, emp_sd, mean_se, se_ratio,
coverage, mcse_cov`, plus `J, frac, N_pop` for the sweep).

**Manifests** in the manifest dir: `manifest_R10_fpc_partA.rds`,
`manifest_R10_fpc_partB_J###.rds` (seeds, package versions, git rev, sessionInfo
bits — same reviewer-proofing as `run_sim.R`).

What to inspect:
- **PART A**: for L1_param and L2_smooth, the CF (and FA) `mean_se` is modestly
  smaller WITH fpc, and `coverage` **stays ≥ 0.94**. The narrative: the FPC is a
  legitimate, available correction that tightens the SE a little at `f=0.25`
  without breaking coverage — so the headline (no-FPC) numbers are if anything
  conservative.
- **PART B**: `coverage` is **flat (≈ the no-fpc R1/L1 level, ~0.91–0.95)** across
  `J = 8 → 64` as `frac` falls `0.25 → 0.031`. Flat coverage = the
  `m/M_pop -> 0` assumption is innocuous for the super-population estimand.

## Aggregation -> figure/table

The two CSVs feed a small supplement table/figure: PART A as a 2×2×2 (rung × arm ×
fpc) coverage/SE table; PART B as a coverage-vs-`f` line (flat line is the result).
Both reuse the locked column names so the existing figure code can ingest them.

## Dependencies

R 4.4/4.5 with the project library: `SuperLearner, tmle, survey, surveyCV, earth,
gam, glmnet, ranger` (same as the locked sim). No new packages. On ARC the
`module load gcc/9.4.0; module load r/4.4.0` + `R_LIBS=$HOME/...` env from
`codes/sim.slurm` is mirrored in `submit.slurm`.

## Caveats / assumptions

- **Single-stage FPC.** `.se_des_fpc()` applies a stage-1 (PSU-level) FPC only,
  matching the survey-design SE the paper uses (`ids = ~cluster`). The second
  stage (units within PSU) FPC is intentionally omitted — the between-PSU term
  dominates the clustered design variance, and adding a 2nd-stage FPC would
  require a two-stage `svydesign` not used elsewhere in the pipeline. The result
  is a *conservative* (slightly larger-than-fully-corrected) SE, which only
  strengthens the "coverage stays ≥0.94" claim.
- **`fpc1 = J_per_stratum` is the POPULATION PSU count per stratum**, constant
  within a stratum and ≥ the number sampled (`m_h = 2`), so `survey` accepts it.
  Verify in the smoke that `fpc=yes` SE ≤ `fpc=no` SE.
- **CF EIF replay.** `cf_eif_for_rep()` re-fits the CF arm to get its EIF (the
  engine returns `eif_fa` but not the CF EIF). The fold split uses fresh RNG, so
  the local no-fpc CF SE may differ slightly from the engine's; the comparison
  that matters (fpc vs no-fpc) uses the SAME local EIF, so it is internally
  consistent. `se_cf_engine_check` is carried per-rep for auditing.
- **`L4` excluded from the FPC branch by design** (see Purpose) — do not "complete
  the grid" by adding it; that would manufacture under-coverage.
- The PART B super-pop truth `Psi` is recomputed per `J` (the population changes),
  but is essentially `J`-invariant by construction; small `Psi` differences across
  `J` are Monte-Carlo truth noise (`truth_M = 2e6`), not a regime effect.
- Outputs are namespaced under `arc_runs/R10_fpc/` and `results/arc/`; rerunning
  overwrites only this run's files.
```
