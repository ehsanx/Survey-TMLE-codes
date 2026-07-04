# R13_null_typeI — Null / type-I error (ψ = 0) + small-effect power

**Spec.** `Writing/comments/phase4-arc-sim-specs.md` item **A5** — *"Null / type-I error
(ψ=0) + small-effect power"*: confirm type-I control at a true null (the E3 "apparent
association dissolves" story is never simulated at ψ=0) and add a small-effect power
point. Knobs: `make_population(te_log_odds = 0)` for the null, `te_log_odds = 0.3` for
power. Both designs (`standard`, `R1`), all rungs L1–L4, all five arms,
`run_estimators()` **unchanged**.

**What this run produces.** Rejection rate of H0: ψ=0 (target 0.05) per
te × scenario × rung × method, plus the standard coverage-of-true-Ψ / bias / emp_sd /
se_ratio columns computed with `aggregate_sim.R`'s exact formulas — feeding the new
Web-D type-I/power table and a §5/§7 sentence, with audit checks on the CF rejection
rate ≈ 0.05.

## Design (locked)

- **Cells:** `te ∈ {0, 0.3}` × `scenario ∈ {standard, R1}` × all 4 `SL_LADDER` rungs
  = **16 cells**. FULL = 1000 reps/cell, `SIM_CHUNK=100` → 10 chunks/cell →
  **`--array=1-160`**. Cells are te-major (tasks 1–80 te=0, 81–160 te=0.3); within te,
  `expand.grid(scenario, rung)` scenario-fastest; chunks fastest within a cell.
- **Population per cell:** `make_population(scenario, model_type="complex",
  te_log_odds=TE, truth_M=2e6L)` (SMOKE: `2e5L`), deterministic from `POP_SEED`.
- **Null truth check (te=0):** verify `|truth$psi| <= 4*truth$se_mc`; printed and
  recorded as `null_ok` in the manifest (**warn, do not stop**). At te=0 the GH blip is
  identically 0 (m1 == m0 bitwise when θ=0), so ψ = se_mc = 0 exactly; the check uses
  `max(4*se_mc, 1e-12)` to stay well-defined in that degenerate case.
- **Estimators:** `run_estimators()` UNCHANGED — Fully-Aware, Partially-Aware,
  Non-Aware, Fully-Aware-CV (multi-learner rungs L2/L3 only), Fully-Aware-CF.
  Every results row is tagged with `te=TE`.
- **Comparability property (preserved):** reps use
  `draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i)` with **i the GLOBAL rep
  index**. `draw_sample` is independent of te, and in `dgp.R` θ enters ONLY the final
  Y draw (after the C and A RNG, and the treatment model has no θ), so for a given
  (scenario, rep) the **sampled units, C, A, weights and design are byte-identical to
  the locked headline run (te=1.5)** — only the Y column differs. Differences vs the
  headline table are therefore attributable to the effect size, not to sampling noise.
- **RNG streams:** `clusterSetRNGStream(iseed = SAMPLE_SEED_BASE + 1000L*cell_index +
  task)` — unique WITHIN this run (the 1000-per-cell offset keeps streams distinct
  across cells/chunks) and disjoint from the OLDER arc runs' `SAMPLE_SEED_BASE + task`
  range (R01/R03/R04/R07/R11/R12). The sibling batch runs (R14/R15/R16/R18/R19, and
  R02 via `1000L*m_idx + task`) share the same formula, so cross-run iseed
  coincidences DO occur — harmless: worker streams only drive SL CV folds / ranger /
  `make_cf_folds`, while the samples self-seed inside
  `draw_sample(SAMPLE_SEED_BASE + global rep index)` regardless of the stream.
- **Output files:** `r13_te{000|030}_<scenario>_<rung>_chunk###.rds` under
  `R13_OUT` (+ per-task manifests in `R13_OUT/manifest/`). SMOKE outputs get a
  `SMOKE_` prefix. Per-chunk checkpoint via `arc_skip_if_done()` (full-array
  resubmission is idempotent); checkpoint disabled in SMOKE.

## Cost estimate

Per-task cost is identical to the headline `codes/run_sim.R` task (same engine, same
per-rep work; the te knob is free). Headline: L4/R1 ≈ 440 s/rep single-core → ~23–45
min per 100-rep chunk at 32 cores; L1–L3 chunks are minutes. This run is exactly **2×
the headline grid** (16 vs 8 cells): 160 tasks × 32 cpus, worst-case wall ≈ 45 min for
the 40 L4 tasks, minutes for the rest → roughly **2.5–4k cpu-hours** total, well inside
one overnight window with `--time=04:00:00` (conservative cap; cheap cells exit early).

## Smoke (run locally BEFORE submitting; ~2–3 min)

One cheap cell (te=0 × standard × L1_param), 20 reps, 2 cores, `truth_M=2e5L`.

PowerShell:
```powershell
$env:SMOKE="1"
$env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R13_DIR="<repo-root>/simulation/enhancements/R13_null_typeI"
$env:R13_OUT="<repo-root>/sim_output/arc_runs/R13_null_typeI"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R13_null_typeI/run.R"
```

bash:
```bash
SMOKE=1 \
  REPO_ROOT="<repo-root>" \
  SIM_CODE="<repo-root>/R" \
  DATA_ROOT="<repo-root>/sim_output" \
  R13_DIR="<repo-root>/simulation/enhancements/R13_null_typeI" \
  R13_OUT="<repo-root>/sim_output/arc_runs/R13_null_typeI" \
  SLURM_CPUS_PER_TASK=2 \
  Rscript codes/arc_runs/R13_null_typeI/run.R
```

Writes `SMOKE_r13_te000_standard_L1_param_chunk001.rds` (+ a `SMOKE_manifest_*` in
`manifest/`) and prints an inline per-method table plus a `[SMOKE-GATE]` verdict line.

### SMOKE-GATE DECISION RULE

- **PASS → submit the full array** iff ALL of:
  1. every per-rep `b` and `se` is finite;
  2. 4–5 method rows present (exactly **4** at L1_param — the CV arm only exists on
     multi-learner rungs);
  3. `|reject_rate − 0.05| ≤ 0.25` for every method (lenient: 20 reps ⇒ rejection
     counts of 0–6 pass for a nominal-0.05 arm);
  4. the null truth check passes (`null_ok=TRUE` printed by the NULL CHECK line).
- **STOP-and-report** if any check fails — in particular a non-zero ψ at te=0
  (DGP/truth bug) or non-finite SEs (engine/separation issue), BEFORE burning the
  160-task run.

Local validation (2026-06-12, R 4.5.1, 41 s wall on 2 cores; numbers identical
across two runs ⇒ stream seeding is deterministic): **[SMOKE-GATE] PASS** — all b/se
finite; 4 method rows; reject_rate FA/CF/PA = 0.00, NA = 0.10 (max |reject−0.05| =
0.05 ≤ 0.25); NULL CHECK `|psi|=0 ≤ max(4·se_mc,1e-12)` OK. CF mean_b = 0.0180,
emp_sd = 0.0251, mean_se = 0.0375 — consistent with the headline L1/standard scale.
aggregate.R was also exercised end-to-end on the smoke RDS (sandboxed under TEMP):
null-truth audit, both decision views, CSV + combined RDS all produced.

## Submit + aggregate (ARC)

```bash
# from the repo root on scratch (= SLURM_SUBMIT_DIR)
sbatch codes/arc_runs/R13_null_typeI/submit.slurm     # FULL --array=1-160
```

```bash
# after all tasks finish (counts: 160 r13_*.rds expected)
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R13_OUT="$PWD/sim_output/arc_runs/R13_null_typeI" \
  Rscript codes/arc_runs/R13_null_typeI/aggregate.R
```

PowerShell (local aggregation after copying outputs back):
```powershell
$env:SIM_CODE="<repo-root>/R"
$env:REPO_ROOT="<repo-root>"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R13_OUT="<repo-root>/sim_output/arc_runs/R13_null_typeI"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R13_null_typeI/aggregate.R"
```

## Expected outputs + KEY NUMBERS to inspect

- **Per task:** `sim_output/arc_runs/R13_null_typeI/r13_te{000|030}_<scenario>_<rung>_chunk###.rds`
  (run id, te, knobs, learners, chunk, reps, Psi, truth, params, `null_ok`, per-rep
  results + diagnostics) + `manifest/manifest_*.rds` (params, truth, seeds incl.
  POP_SEED/SAMPLE_SEED_BASE/TRUTH_SEED + the stream iseed, R/package versions,
  timestamp, git sha, sysname, `null_ok`).
- **Aggregated:** `results/arc/R13_null_typeI_summary.csv` + 
  `R13_OUT/R13_null_typeI_combined.rds`. Columns: `te, scenario, rung, method, n_reps,
  Psi, bias, emp_sd, mean_se, se_ratio, coverage, mcse_cov, reject_rate, mcse_reject,
  deff_clust, icc_eif`.
- **KEY NUMBERS:**
  - **Decision view 1 (type-I at te=0; the audit object):** `reject_rate` for the
    **Fully-Aware-CF** rows — target **≈ 0.05** (flag band
    0.05 ± 2·sqrt(.05·.95/1000) ≈ [0.036, 0.064]). Expected pattern elsewhere:
    PA/NA **over-reject** (unclustered/iid SEs; headline coverage 0.86–0.46 ⇒ type-I
    up to ~0.5 for NA), single-fit FA over-rejects at **L4** (Donsker failure; headline
    L4 coverage 0.17–0.22 at te=1.5, so expect severe over-rejection if the L4
    single-fit bias persists at the null). **Reading a CF FLAG:** under the locked
    `model_type="complex"` DGP the CF arm carries the known nuisance-misspecification
    bias at **L1–L3** (~ +0.018 at the null, matching the headline bias rows), so the
    te=0 rejection of H0: ψ=0 conflates bias with SE calibration — CF rows can leave
    the [0.036, 0.064] band **upward** from that bias (L1–L3), or **downward** from
    the conservative CF SEs (`se_ratio` ≈ 1.4 ⇒ effective critical value ≈ 2.8
    empirical-SD units). A FLAG is therefore NOT by itself a type-I failure: read it
    **jointly with the `bias` and `se_ratio` columns** (R01 simple-control is the
    matched evidence separating nuisance misspecification from variance
    miscalibration). The cleanest audit-target cell is **CF at L4** (bias ≈
    0.002–0.006): a flag there that bias + se_ratio cannot explain ⇒ STOP-and-report
    before writing Web-D prose.
  - **Decision view 2 (power at te=0.3):** FA/CF/CV rejection rates. At te=0.3 the RD
    truth Ψ is ≈ 0.04–0.05 (vs 0.21 at te=1.5) with emp_sd ≈ 0.035–0.045, so expect
    **modest power** (~0.15–0.40); the point is a calibrated power VALUE, not high
    power. CF power may sit slightly below FA (its SE is wider by design at L4).
  - The `null-truth audit` table printed by aggregate.R: every te=0 population row
    must show `null_ok=TRUE` (ψ = 0 exactly by GH construction).

## Caveats / assumptions

- **Rejection test:** `|b|/se > qt(.975, max(1, df))` — same t-reference as the
  coverage convention; at te=0, coverage of Ψ=0 equals 1 − reject_rate by construction
  (a built-in consistency check between the new and standard columns).
- **Flag band uses the nominal-rate MCSE** `sqrt(.05*.95/n)`: the spec's "0.05 ±
  2*MCSE" is implemented with the MCSE under the nominal rate, since the empirical
  `mcse_reject` collapses to 0 when a cell rejects never/always (it is still reported
  as a CSV column).
- **PA/NA "power" at te=0.3 is not honest power** (their SEs are invalid; their type-I
  at te=0 is the relevant display) — the power table prints all arms for completeness
  with a note.
- **te=0.3 truth:** ψ(0.3) is computed by the same 2e6-draw GH/MC integral as the
  headline; the small |Ψ| ≈ 0.04 makes the relative MC error larger in relative terms
  but se_mc ≈ 4e-5 is still ~1000× smaller than emp_sd — negligible for coverage.
- **L1/L4 rungs have no CV arm** (single-learner libraries) — 4 methods there, 5 at
  L2/L3; aggregate.R groups by method so the tables are simply shorter for those rungs.
- **Non-destructive:** sources `codes/*.R` read-only; writes only under
  `sim_output/arc_runs/R13_null_typeI/` and `results/arc/R13_null_typeI_summary.csv`.
  The locked `results/sim_full_summary.csv` and `sim_output/intermediate/` are never
  touched.
- Dependencies are the canonical set (`SuperLearner, tmle, survey, surveyCV, earth,
  gam, glmnet, ranger`); ARC env = `gcc/9.4.0 + r/4.4.0 + R_LIBS=$HOME/R/...`.
