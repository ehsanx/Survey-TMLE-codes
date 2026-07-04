# R17_e4_floor — E4 floor sensitivity + share-at-floor for N6 (spec A9b)

**Spec.** Item **A9b — "E4 floor sensitivity + overlap-restricted estimand (N6)"** in
`Writing/comments/phase4-arc-sim-specs.md`: for E4 (6% exposure; single-fit min ĝ → ~0.003
on the N1-corrected data) the locked 0.05 propensity floor **binds**; re-run the primary
arm at floor ∈ {0.05, 0.025, 0.01}, report the point estimate, the **share of units at the
floor**, and the **weighted exposure mass below 0.05**; plus the share-at-floor per example
(E2–E4) for the N6 estimand-reframe wording (truncated / overlap-restricted target where
the floor binds).

## What it runs

ONLY the **primary** estimator — Fully-Aware-CF at the primary library
`LIB = c("SL.glm","SL.earth","SL.glmnet")` (== `Nhanes/R/03_run_estimators.R`):
PSU-within-stratum `make_cf_folds`, **unweighted** per-fold `.sl` fits, OOF propensity
clipped at the **cell's floor**, OOF Q clipped at the locked `1e-3`, weighted **pooled**
`tmle()` targeting on the domain, EIF reconstructed with `gbound = the cell's floor`,
domain design SE via the full-design `subset()` path (`.se_des(..., nest=TRUE, inpop)`).
B = `NHANES_B` (default 20) split-repeats per cell, split seeds `20260607 + 1:B`
(the locked `nhanes_arc.R` convention).

| task | example | floor | role |
|---|---|---|---|
| 1 | E2 | 0.05 | N6 share row |
| 2 | E3 | 0.05 | N6 share row |
| 3 | E4 | 0.05 | N6 share row + sensitivity anchor (locked floor) |
| 4 | E4 | 0.025 | sensitivity |
| 5 | E4 | 0.01 | sensitivity |

**Per-split extra records** (computed from the **RAW pre-clipping** OOF ĝ over the domain
rows, with domain weights `w`): `share_clip_w = Σw·1(ĝ_raw<floor)/Σw`,
`share_clip_unw = mean(ĝ_raw<floor)`, `mass05_w = Σw·1(ĝ_raw<0.05)/Σw` (**fixed 0.05
reference** regardless of cell floor), `expmass05_w = Σw[A=1, ĝ_raw<0.05]/Σw[A=1]` (share
of the **exposed** below 0.05), `g_raw_min/max`, plus `b, se, df` per split. Lower tail
only (rare exposure ⇒ the floor binds from below); the upper tail is monitored via
`g_raw_max` and `share_clip_hi_w`.

## Files

- `run.R` — array driver (1 task = 1 cell). Sources `codes/{estimators,diagnostics}.R` +
  `codes/arc_runs/_checkpoint.R` read-only; mirrors `R05_harmonized_floor/run.R`
  structurally (NHANES loading + one-hot encoding byte-identical to `nhanes_arc.R`,
  `nest=TRUE`, inpop domain SE, B-split convention, per-cell checkpoint).
- `estimators_r17.R` — ALL new logic: `run_cf_floor()`, the primary CF arm with a
  tunable `g_floor` + the raw-ĝ floor/share records. Reuses the exposed engine blocks
  `.sl`, `make_cf_folds`, `.eif_from_tmle(gbound=)`, `.se_des` — nothing reimplemented.
- `aggregate_R17.R` — per-cell RDS → the two CSVs + combined RDS + decision views.
- `submit.slurm` — SLURM array (mirrors R05's; only array/driver/R17 paths changed).

## Key design facts (read before interpreting)

- **Floor really binds in targeting — by source inspection, not by a runtime check.**
  `tmle` 2.1.1's internal default `gbound = NULL` resolves to `5/√n/log n` ≈ 0.0051
  (E3/E4) / 0.0028 (E2) with ATE bounds `[lb, 1]` — **below every cell floor**
  (0.0051 < 0.01 < 0.025 < 0.05) — so `.bound` is a no-op on the pre-clipped `g1W`
  and the floor we set is the floor used. NOTE the SMOKE-gate assert
  `min(cf$g$g1W) ≥ floor` canNOT verify this: `tmle` returns the user-supplied `g1W`
  unmodified in `cf$g$g1W` (its internal bounding applies to an unexposed
  `g1W.total`), so the assert is tautologically true for any pre-clipped input. It is
  kept only as a cheap sanity invariant on our own clipping; the floor-correctness
  claim rests on the source-inspection numbers above.
- **Exact floor attribution.** `set.seed(20260607+b)` is called per split and
  `make_cf_folds` is the **first** RNG consumer in `run_cf_floor`, so for a given split
  the folds, the per-fold SL fits, and the **raw OOF ĝ are identical across the three E4
  floor cells**; the floor enters only via post-hoc clipping + targeting. Differences in
  `b` across floors are the pure floor effect (no split noise between floor cells).
- **Comparability to the locked run.** The analytic frames are FIXED data — byte-identical
  inputs to the locked headline ladder and the ITEM-0 E4 re-run (no `draw_sample`; the
  sim packages' byte-identical-samples property translates here to byte-identical input
  frames + the same split-seed base). Caveat: realized CF folds differ from the locked
  *all-arm* run at the same seed, because there FA/NA/CV consume RNG before the CF fold
  draw; the E4/f050 cell is therefore a fresh primary-arm estimate, not a byte-replay of
  the locked CF row (expect it within the locked `b_split_sd` of the locked E4 CF value).
- **Provenance.** `Nhanes/analytic/E*_imputed.rds` carry the **N1-corrected** data (post
  RHD180 fix, commit `dcce692`; pre-fix copies in `Nhanes/analytic_backup_N1/`). The E4
  numbers here therefore extend the ITEM-0 refreshed E4 results, not the pre-fix ones.

## Cost estimate

Per cell: B=20 splits × 1 arm × (V_eff≈2 folds × 2 nuisances × 3-learner SL with 10-fold
internal CV). The locked 5-arm ladder cell at this scale ran ~1–2 h on 32 cores; one arm
≈ 20–45 min/cell (E2, n_domain≈29k, heaviest; E3/E4 ≈11k light). FULL array = 5 cells ×
≤1 h ≈ **≤5 node-hours**; `--time=04:00:00` is generous. SMOKE (GLM-only, B=2, E4) ≈ 1–3
min locally.

## SMOKE (run locally BEFORE submitting; wiring check, NOT the primary library)

PowerShell:
```powershell
$env:SMOKE="1"
$env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R17_DIR="<repo-root>/Nhanes/sensitivity/R17_e4_floor"
$env:R17_OUT="<repo-root>/Nhanes/nhanes_output/arc_runs/R17"
$env:NH_ANA="<repo-root>/Nhanes/analytic"
$env:R17_LIB_OVERRIDE="SL.glm"
$env:NHANES_B="2"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/Nhanes/sensitivity/R17_e4_floor/run.R"
```

bash (ARC login node, from the repo root):
```bash
export SMOKE=1 REPO_ROOT=$PWD SIM_CODE=$PWD/codes DATA_ROOT=$PWD/sim_output \
  R17_DIR=$PWD/Nhanes/arc_runs/R17_e4_floor \
  R17_OUT=$PWD/Nhanes/nhanes_output/arc_runs/R17 \
  NH_ANA=$PWD/Nhanes/analytic \
  R17_LIB_OVERRIDE=SL.glm NHANES_B=2 SLURM_CPUS_PER_TASK=2
Rscript "$R17_DIR/run.R"
```

SMOKE forces the **E4/f050** cell, B=2, ≤2 cores, GLM-only library; outputs get a
`SMOKE_` filename prefix (never collide with real cells; safe to leave in place).
Target < 8 min (measured ~1–2 min).

### SMOKE-GATE DECISION RULE

The run prints `[SMOKE-GATE] PASS` or `[SMOKE-GATE] STOP: <reason>` after these checks:

**PROCEED to the full ARC submit iff ALL hold** (the gate asserts them):
- finite `b`, finite `se > 0`, `df > 0` for every split; no split errored;
- `share_clip_w`, `share_clip_unw`, `mass05_w`, `expmass05_w` all in `[0,1]`;
- `g_raw` range sensible (`0 < min < max < 1`) **and the floor binds on E4**
  (`g_raw_min < 0.05` — with GLM-only ĝ this should be comfortably true; post-N1
  single-fit min ĝ ≈ 0.003);
- `min(cf$g$g1W) ≥ floor` (cheap sanity invariant on our own clipping ONLY — `tmle`
  echoes the supplied `g1W` unmodified, so this line cannot detect a passthrough
  failure; the f025/f010 cells are protected by the source-inspection argument in
  "Key design facts", not by this assert);
- `eps_cf < 20` (pooled weighted targeting did not diverge).

**STOP and fix (do NOT submit) if** the verdict line is STOP — most likely causes:
a covariate-encoding mismatch (assert fires), a worker missing an exported object,
or a clipping bug in `estimators_r17.R` (the `g_used` sanity invariant fails — note
it would NOT catch tmle re-bounding the supplied ĝ, since `cf$g$g1W` echoes the
input; that risk is excluded by the source-inspection numbers, not by the gate).
Fix in `estimators_r17.R`, never in `codes/`.

Note: SMOKE's `b` is a GLM-only wiring number, NOT comparable to the primary table;
do not read it as science.

## Submit + aggregate (ARC)

```bash
sbatch Nhanes/arc_runs/R17_e4_floor/submit.slurm          # FULL --array=1-5
# test variant: sbatch --array=3 ... (single E4/f050 cell)
# after the array finishes:
Rscript Nhanes/arc_runs/R17_e4_floor/aggregate_R17.R
```
Re-submitting the full array is idempotent (per-cell `arc_skip_if_done`).

## Expected outputs + KEY NUMBERS

- 5 cell files `Nhanes/nhanes_output/arc_runs/R17/nh17_{E2,E3,E4}_{f050,f025,f010}.rds`
  (each: `summary`, `per_split` (B rows), `diagnostics` incl. the split-1 raw-ĝ vector
  for overlap/N6 figures) + `manifest/`.
- `results/arc/R17_floor_sensitivity.csv` — the 3 E4 rows across floors:
  `floor, b, b_split_sd, se, df, lcl, ucl, share_clip_w, mass05_w, expmass05_w`.
- `results/arc/R17_floor_share.csv` — E2/E3/E4 at floor 0.05 (the N6 share-at-floor
  table): `example, n_domain, A_prev, share_clip_w, share_clip_unw, mass05_w,
  expmass05_w, g_raw_min, g_raw_max, b, se, lcl, ucl`.
- `Nhanes/nhanes_output/arc_runs/R17/R17_combined.rds`.

**Key numbers to inspect (the aggregator prints both decision views):**
1. **E4 stability across floors**: `max |b(f_i) − b(f_j)|` vs the mean design `se`
   (ratio ≪ 1 ⇒ the E4 estimate is floor-stable and the N6 reframe is precautionary
   wording, not a correction; ratio ≳ 1 ⇒ the floor materially moves E4 — report this
   to the authors BEFORE using the E4 row, since the manuscript's truncated-estimand
   wording must then carry the table). Same-split raw ĝ is identical across floor
   cells, so this difference is the pure floor effect.
2. **N6 share table**: expect `share_clip_w`/`expmass05_w` to be material on E4
   (rare exposure; locked diagnostics had a non-trivial near-bound mass) and ≈0 on
   E2/E3. `mass05_w` is the cross-floor-comparable reference (fixed 0.05 cut).
3. `share_clip_w` must DECREASE monotonically as the floor drops 0.05→0.025→0.01 on
   E4 (it counts `ĝ_raw < floor`); `mass05_w` must be (split-noise) constant across
   the E4 cells — a built-in consistency check on the raw-ĝ wiring.

## Caveats / assumptions

- **Single imputation** (m=1) inherited from `02_impute.R`; SEs understate imputation
  uncertainty (stated paper caveat; the MI sensitivity lives in R06).
- `b_split_sd` is across B cross-fit/SL splits; `se` is the mean design SE (locked
  `nhanes_arc.R` convention). `lcl/ucl` use `qt(.975, max(1, df))`.
- Share-at-floor counts the **lower** tail only (the spec's quantities); the upper tail
  (`ĝ_raw > 1−floor`) is recorded per split as `share_clip_hi_w` and via `g_raw_max`
  but is expected ≈0 (rare exposure).
- The E4/f050 cell is a fresh primary-arm estimate on the N1-corrected data, not a
  byte-replay of the locked E4 CF row (different RNG path; see "Key design facts").
- `g_raw_min/max` in the CSVs are means of per-split min/max; the per-split values are
  in `per_split`, the full split-1 ĝ vector in each cell's `diagnostics$g_raw`.
- E2's `mass05_w`/`expmass05_w` use the same fixed 0.05 cut as its cell floor, so for
  f050 cells `mass05_w == share_clip_w` by construction (both reported for uniformity
  with the E4 f025/f010 cells, where they differ).
