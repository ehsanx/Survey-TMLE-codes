# R05_harmonized_floor — NHANES harmonized-floor + weighted-OOF-CF ladder

De-confound the real-data (NHANES) ladder. The locked NHANES pipeline
(`Nhanes/R/nhanes_arc.R` → `codes/run_estimators`) compares **Fully-Aware**
(single weighted `tmle()` fit; EIF re-bounded at `gbound = 1e-3`) against
**Fully-Aware-CF** (cross-fit, per-fold nuisances **unweighted**, OOF propensity
floored at `0.05`). Those two arms differ on **three** axes at once:
(single-fit vs cross-fit) × (weighted vs unweighted nuisance) × (`1e-3` vs `0.05`
floor). This run isolates them on the **E3** and **E4** ladders (the examples with
the full L1–L4) by adding two arms that change **one axis at a time**, without
editing any canonical file:

| arm | what it is | isolates |
|---|---|---|
| `Fully-Aware` | locked FA, EIF floor `1e-3` (re-emitted unchanged) | baseline |
| `Fully-Aware-h05` | **same** FA `tmle()` fit, EIF re-bounded at **0.05** | the **floor** axis |
| `Fully-Aware-CF` | locked CF, per-fold **unweighted** OOF, floor `0.05` (re-emitted) | cross-fit baseline |
| `Fully-Aware-CF-wOOF` | cross-fit, per-fold nuisances **weighted** (`obsWeights`), floor `0.05` | the **de-weighting** axis |

**Expected (the framing the manuscript should lock):**
1. Harmonizing the FA floor to `0.05` **leaves the FA point estimate unchanged**
   (`b_FA == b_FA_h05` exactly — FA `psi` comes from the targeted `Qstar`, not the
   EIF), so the qualitative single-fit-vs-CF **point divergence survives** the
   floor change. Only the FA SE/CI width moves with the floor.
2. The weighted-OOF-CF arm shows that **cross-fitting** (not de-weighting) is what
   stabilizes the ladder: if `b_CF-wOOF` tracks `b_CF` (not the single-fit `b_FA`),
   the stabilization is attributable to cross-fitting per se.

## Files
- `run.R` — driver. Sources `codes/{estimators,diagnostics,learners}.R` read-only,
  sources this folder's `estimators_nhanes_iso.R`, mirrors `nhanes_arc.R`'s
  split-repeat (B) + diagnostics structure. Grid `{E3,E4} × {L1..L4}` = 8.
- `estimators_nhanes_iso.R` — helper: `run_estimators_iso()`, a thin
  reimplementation of `run_estimators()`'s FA + CF arms that **reuses** the exposed
  building blocks `.eif_from_tmle` (tunable `gbound`), `make_cf_folds`, `.sl`,
  `.se_des`. Adds the two new arms. `inpop`/`nest`/`W_cols` domain handling is
  byte-identical to `run_estimators()` (same `sub <- which(inpop)`, same `se_des`
  closure), so design df and the domain variance match the locked pipeline.
- `aggregate_R05.R` — combines per-cell RDS → `results/arc/R05_*.csv` + combined RDS.
- `submit.slurm` — SLURM array (mirrors `Nhanes/nhanes.slurm`; only array/driver/
  `NH_OUT` changed).

## New logic vs the canonical engine
All in `estimators_nhanes_iso.R`. (1) `Fully-Aware-h05`: calls the **already-exposed**
`.eif_from_tmle(fa, ..., gbound = fa_gbound)` with `fa_gbound = 0.05` on the **same**
`fa` `tmle()` object — no re-fit. (2) `Fully-Aware-CF-wOOF`: re-uses the **same** fold
object from `make_cf_folds`, but the per-fold `.sl(..., weights = w[tr], ...)` is
**weighted** (the locked CF arm passes `weights = NULL`); OOF propensity floored at
`g_oof_bound = 0.05` exactly as the locked CF arm. The locked `Fully-Aware` (1e-3) and
`Fully-Aware-CF` (unweighted-OOF) arms are re-emitted **unchanged** as comparators.
Diagnostics extend the canonical `drow` with `eps_cfw`, `g_cfw_*`, `fa_gbound`,
`g_oof_bound`, and floor-sensitive near-bound masses; the EIF passed to `deff_clust`
is the canonical `1e-3` FA EIF (so DEFF matches the locked numbers), with a parallel
harmonized-floor DEFF for sensitivity.

## Non-destructive outputs
- per-cell RDS → `Nhanes/nhanes_output/arc_runs/R05/nh_<E>_<rung>.rds`
- manifests   → `Nhanes/nhanes_output/arc_runs/R05/manifest/`
- summaries   → `results/arc/R05_summary.csv`, `R05_point_divergence.csv`,
  `R05_diagnostics.csv`, `R05_combined.rds`

Never touches the locked `Nhanes/nhanes_output/{intermediate,results}` or any
`codes/*.R` / `Nhanes/R/*.R`. Column names match the locked NHANES tables
(`example,label,rung,method,b,b_split_sd,se,df,lcl,ucl`) plus the new arms.

## Smoke test (validate locally before submitting)
SMOKE mode forces the **decision-gate cell E4 / L4_aggressive, B=2** (the deepest
non-Donsker rung, where the single-fit-vs-CF divergence is largest), ~2–4 min:

```powershell
# Windows (R 4.4/4.5 with SuperLearner, tmle, survey, surveyCV, earth, gam, glmnet, ranger)
$env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="$env:REPO_ROOT/codes"
$env:R05_DIR="$env:REPO_ROOT/Nhanes/arc_runs/R05_harmonized_floor"
$env:NH_ANA="$env:REPO_ROOT/Nhanes/analytic"
$env:NH_OUT="$env:REPO_ROOT/Nhanes/nhanes_output/arc_runs/R05_smoketest"
$env:NH_MANIFEST="$env:NH_OUT/manifest"
$env:SLURM_CPUS_PER_TASK="2"; $env:SMOKE="1"
Rscript "$env:R05_DIR/run.R"
```

```bash
# Linux/ARC login node
export REPO_ROOT=$PWD SIM_CODE=$PWD/codes \
  R05_DIR=$PWD/Nhanes/arc_runs/R05_harmonized_floor \
  NH_ANA=$PWD/Nhanes/analytic \
  NH_OUT=$PWD/Nhanes/nhanes_output/arc_runs/R05_smoketest \
  NH_MANIFEST=$PWD/Nhanes/nhanes_output/arc_runs/R05_smoketest/manifest \
  SLURM_CPUS_PER_TASK=2 SMOKE=1
Rscript "$R05_DIR/run.R"
```

A wiring-only check that finishes in <1 min (GLM rung, E4/L1, no SMOKE): set the
same env **without** `SMOKE`, add `NHANES_B=2` and `SLURM_ARRAY_TASK_ID=2`
(task 2 = E4/L1_param), and run `run.R`.

Then aggregate the smoke cell (writes to a throwaway results dir):
```
NH_OUT=.../arc_runs/R05_smoketest R05_RESULTS=.../results/arc_smoketest \
  Rscript Nhanes/arc_runs/R05_harmonized_floor/aggregate_R05.R
```

### Smoke-gate decision rule (READ BEFORE the full submit)
Inspect the `[DECISION]` line and the `R05_point_divergence.csv` for the smoke
cell. The framing is **locked on the point-estimate divergence, which is
floor-robust** — NOT on the near-bound mass (`g_fa_near_bound`), which is
floor-sensitive and must not drive the decision.

**PROCEED to the full run iff:**
- `|b_FA − b_FA_h05| ≈ 0` (≤ ~1e-6). This is a built-in invariant: re-flooring the
  FA EIF cannot move the FA point estimate. A non-zero value means the helper is
  wrong (it must reuse the SAME `tmle` fit) — **STOP and fix `estimators_nhanes_iso.R`.**
- `|b_FA − b_CF|` at L4 is materially **> 0** (the headline single-fit-vs-CF
  divergence persists under the harmonized `0.05` floor). In the locked table,
  E4/L4 has `b_FA = −0.338` vs `b_CF = −0.023` → divergence ≈ 0.32, so expect a
  large value here.

**STOP-and-report (do NOT launch the full array) if:**
- `|b_FA − b_FA_h05|` is **not** ~0 → helper bug (see above).
- Harmonizing the floor **collapses** the FA-vs-CF point divergence at L4 (i.e.
  `|b_FA − b_CF|` becomes small, comparable to its B-split SD). That would mean the
  paper's divergence was a **floor artifact**, not a cross-fitting effect — a
  finding that changes the manuscript framing and must be reported before spending
  the full array.
- Any arm errors, returns `NA`/non-finite `b`/`se`, or the run aborts.

Note: the **weighted-OOF-CF** arm (`Fully-Aware-CF-wOOF`) can be numerically
**unstable on the parametric rungs** (weighted per-fold WLS under NHANES's heavy
weights + rare exposure in E4 can extrapolate; the L1 wiring check produced
`b ≈ −0.5`). That instability is itself part of the message (weighted-OOF is the
fragile axis); it is **not** a smoke-gate failure. The gate keys only on the two
FA-vs-CF point quantities above.

## Submit (ARC)
From the repo root (= `SLURM_SUBMIT_DIR` = scratch checkout):
```
sbatch Nhanes/arc_runs/R05_harmonized_floor/submit.slurm
```
- `--array=1-8` runs all `{E3,E4} × {L1..L4}`. To run only the data-adaptive rungs
  L3,L4 (where divergence is largest): edit to `--array=5-8`.
- After the array finishes:
  `Rscript Nhanes/arc_runs/R05_harmonized_floor/aggregate_R05.R`

## Expected output (full run)
- 8 cell files `Nhanes/nhanes_output/arc_runs/R05/nh_{E3,E4}_{L1..L4}.rds`, each a
  list with `summary` (4 arms × `example,label,rung,method,b,b_split_sd,se,df,lcl,ucl`),
  `per_split` (B × 4 rows), and `diagnostics`.
- `results/arc/R05_summary.csv` — 32 rows (8 cells × 4 arms).
- `results/arc/R05_point_divergence.csv` — 8 rows: `b_FA, b_FA_h05, b_CF, b_CF_wOOF,
  d_FA_CF, d_FA_FAh05, d_CFwOOF_CF, d_CFwOOF_FA, cfwoof_closer_to`.
- `results/arc/R05_diagnostics.csv` — overlap (`g_fa/g_cf/g_cfw` ranges), DEFF at
  both floors, near-bound mass, targeting epsilons.

**Key numbers to inspect:** for each E×rung, `d_FA_FAh05` (must be 0),
`d_FA_CF` (the headline; should stay large at L3/L4 under the harmonized floor),
and `cfwoof_closer_to` (expect `"CF"` if cross-fitting — not de-weighting — drives
the stabilization).

## Aggregation
`aggregate_R05.R` reads `NH_OUT` (the R05 subfolder), orders arms
`FA → FA-h05 → FA-CF → FA-CF-wOOF` and rungs `L1→L4`, and writes the three CSVs +
combined RDS to `results/arc/` (override with `R05_RESULTS`). It computes the
floor-robust divergence quantities used by the decision rule.

## Dependencies
- Engine (read-only): `codes/estimators.R` (exposes `.eif_from_tmle`,
  `make_cf_folds`, `.sl`, `.se_des`, `run_estimators`), `codes/diagnostics.R`
  (`deff_clust`), `codes/learners.R` (`SL_LADDER`, custom RF wrappers).
- Data: `Nhanes/analytic/{E3,E4}_imputed.rds` (from `Nhanes/R/01_*`/`02_impute.R`),
  each with `attr(,"covs")`, `attr(,"example")`, `inpop`, and the design columns
  `SDMVSTRA, SDMVPSU, WTMEC_POOLED, A, Y`. **E3/E4 must be imputed first.**
- R 4.4 (ARC) / 4.5 (dev) + `SuperLearner, tmle, survey, surveyCV, earth, gam,
  glmnet, ranger` (and `mice` only for the upstream impute step).
- Reproducibility: split seeds `20260607L + 1:B` (same base as `nhanes_arc.R`);
  a manifest per cell records package versions, git rev, `fa_gbound`, `g_oof_bound`.

## Caveats
- **Single imputation** (m=1) is inherited from `02_impute.R`; this run does not add
  Rubin's rules, so SEs understate imputation uncertainty (a stated paper caveat).
- The **weighted-OOF-CF** arm is included as a *diagnostic contrast*, not a
  recommended estimator; it can be unstable on parametric/low-overlap rungs (see the
  L1 note above). Read its `b_split_sd` alongside `b`.
- `b_split_sd` is across B cross-fit/SL splits; the reported `se` is the
  design-based SE from split-1's `df` averaged over splits (mirrors `nhanes_arc.R`).
- E3 has **no `Fully-Aware-CV` arm** in the locked table because it is a single-
  learner rung at L1/L4; this run does not emit FA-CV at all (it is not part of the
  de-confounding contrast). That is intentional — the four arms here are exactly the
  FA/CF axis isolation.
- Smoke writes under `…/arc_runs/R05_smoketest/` and `results/arc_smoketest/`
  (throwaway) so it never collides with the real `…/arc_runs/R05/` outputs.
```
