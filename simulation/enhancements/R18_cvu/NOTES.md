# R18_cvu — Unweighted-nuisance internal-CV simulation arm (CV-u)

**Spec.** Item **A4** in `Writing/comments/phase4-arc-sim-specs.md`
("Unweighted-nuisance internal-CV simulation arm (CV-u)").

**Why.** The APPLICATION's CV arm fits its nuisances **unweighted** (the
real-data path of `run_estimators()`, taken when `inpop` is supplied — see the
arm-4 note in `codes/estimators.R`), but the locked SIMULATION's CV foil fits
them **weighted** (CV-w). A reviewer can therefore object that the simulation's
"internal CV ≠ cross-fitting" demonstration was made on a *different object*
than the one the application reports. This run adds the unweighted CV variant
(**Fully-Aware-CVu**) to the simulation, on the SAME samples as the locked
headline run, so the contrast is shown for the application's actual object.

**Design (locked).** `run_cvu()` in `cvu_arm.R` mirrors the engine's
Fully-Aware-CV block (`codes/estimators.R` lines ~160–182) **exactly** except
that `obsWeights` for BOTH nuisance SuperLearner fits is `rep(1, n)`. Downstream
is byte-identical to the engine CV arm: same `surveyCV::folds.svy` cluster-aware
internal-CV folds, same Q truncation `[1e-3, 1-1e-3]`, same `g_oof_bound = 0.05`
floor, same **weighted** `tmle()` targeting (`obsWeights = w`, `Q`/`g1W`
supplied), EIF via `.eif_from_tmle(..., gbound = 1e-3)` (engine default — keeps
rows comparable to the locked CV-w/CF rows), clustered design SE + df via
`.se_des`. Only this ONE arm is fit per rep (cheap).

In this DGP sampling is ignorable given C (selection depends only on
stratum/PSU; A, Y depend on the stratum only through C), so the unweighted
nuisance fits are *consistent* — any CVu under-coverage at L3 is attributable to
the **lack of sample-splitting**, not to weighting bias. That is exactly the A4
claim.

**Comparability (do not break).** Each rep `i` draws its sample via
`draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i)` with `i` the GLOBAL rep
index — the samples are **byte-identical to the locked headline run** for the
same scenario, so CVu rows are directly comparable to the locked
Fully-Aware-CV / Fully-Aware-CF / Fully-Aware rows on the same draws.

## Files

| file | role |
|---|---|
| `run.R` | driver; sources the canonical engine read-only; one (scenario×rung×chunk) cell per array task |
| `cvu_arm.R` | ALL new logic: `run_cvu()` (the CV-u arm; reuses `.eif_from_tmle`, `.se_des`) |
| `aggregate.R` | per-task RDS → `results/arc/R18_cvu_summary.csv` + decision views vs the locked table |
| `submit.slurm` | SLURM array (mirrors R04 resources; `--time=02:00:00`, `--array=1-40`) |

## Grid / cost

- Cells = 2 scenarios (`standard`, `R1`) × 2 rungs (`L2_smooth`, `L3_adaptive`)
  = 4 (the CV arm exists only on multi-learner rungs; L1/L4 are single-learner).
- FULL = 1000 reps/cell, `SIM_CHUNK=100` → 4 × 10 = **40 array tasks**
  (`--array=1-40`). Task order (cell-major): 1–10 standard×L2, 11–20 R1×L2,
  21–30 standard×L3, 31–40 R1×L3. TEST submission: `--array=1,11,21,31`.
- Per rep: 2 SuperLearner fits (3-learner library, ~5-fold internal CV) + a
  cheap pre-fit-Q/g `tmle()` targeting ≈ 20–60 s single-core. 100 reps over 32
  cores ≈ **5–15 min/task** → the 2 h wall is generous. Whole run ≈ 40 tasks ×
  ≲0.25 h × 32 cpus ≈ **≤ 320 core-hours** (upper bound).
- Re-submitting the full array is idempotent (`arc_skip_if_done` per chunk).
- RNG streams: `clusterSetRNGStream(iseed = SAMPLE_SEED_BASE + 1000*cell_index
  + task)` — unique within the run and offset away from sibling runs that use
  `SAMPLE_SEED_BASE + task`. Streams affect only learner-internal RNG
  (SuperLearner CV splits, `folds.svy`, ranger); the SAMPLES are seeded
  independently per rep (see Comparability).

## Smoke test (run locally before submitting; target < 8 min)

One cheap cell (`standard × L2_smooth`), 10 reps, 2 cores, `truth_M = 2e5`.

PowerShell:
```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R18_DIR="<repo-root>/simulation/enhancements/R18_cvu"
$env:R18_OUT="<repo-root>/sim_output/arc_runs/R18_cvu"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R18_cvu/run.R"
```

bash (repo root):
```bash
SMOKE=1 \
  REPO_ROOT="<repo-root>" \
  DATA_ROOT="<repo-root>/sim_output" \
  SIM_CODE="<repo-root>/R" \
  R18_DIR="<repo-root>/simulation/enhancements/R18_cvu" \
  R18_OUT="<repo-root>/sim_output/arc_runs/R18_cvu" \
  SLURM_CPUS_PER_TASK=2 \
  Rscript codes/arc_runs/R18_cvu/run.R
```

Writes `sim_output/arc_runs/R18_cvu/SMOKE_r18_standard_L2_smooth_chunk001.rds`
(+ `manifest/SMOKE_manifest_...`), prints an inline summary row and an explicit
`[SMOKE-GATE] PASS/STOP` line.

### SMOKE-GATE DECISION RULE
- **PASS** (submit the full `--array=1-40`) iff ALL of:
  1. every per-rep `b` and `se` is finite;
  2. the reported df match the design df (standard: 60 PSUs − 10 strata = **50**;
     tolerance ±2);
  3. coverage over the 10 reps ∈ **[0.6, 1]** (10-rep binomial noise around the
     expected ~0.94 at L2).
- **STOP** and investigate otherwise (e.g. `se` = NA → `.se_des`/design issue;
  df ≪ 50 → wrong clustering columns; coverage < 0.6 → a bias bug in the arm,
  since L2 CV-w/CF/FA all cover ~0.94–0.95 in the locked table).

## Submit (full run on ARC)

```bash
# from the SLURM_SUBMIT_DIR (the repo root on scratch)
sbatch codes/arc_runs/R18_cvu/submit.slurm
```

### Aggregate after the run
```bash
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R18_OUT="$PWD/sim_output/arc_runs/R18_cvu" \
  Rscript codes/arc_runs/R18_cvu/aggregate.R
```
(PowerShell: set the same four env vars, then
`& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" codes/arc_runs/R18_cvu/aggregate.R`.)

## Expected outputs + KEY NUMBERS

- **Per task:** `sim_output/arc_runs/R18_cvu/r18_<scenario>_<rung>_chunk###.rds`
  (per-rep `b/se/df`, DEFF/ICC on the CVu EIF, design checks) +
  `manifest/manifest_r18_..._chunk###.rds` (params, truth, seeds, package
  versions, git sha).
- **Aggregated:** `results/arc/R18_cvu_summary.csv` (+
  `R18_cvu_combined.rds` under the OUT dir) and TWO printed decision views:
  the long CVu-vs-locked table and the wide coverage table
  (`cov_CVu, cov_CVw, cov_CF, cov_FA, CVu_minus_CVw, CVu_minus_CF`).
- **Locked comparators** (from `results/sim_full_summary.csv`, 1000 reps):

  | scenario | rung | CV-w | CF | FA |
  |---|---|---|---|---|
  | standard | L2_smooth   | 0.940 | 0.952 | 0.942 |
  | standard | L3_adaptive | 0.847 | 0.951 | 0.895 |
  | R1       | L2_smooth   | 0.943 | 0.941 | 0.940 |
  | R1       | L3_adaptive | 0.877 | 0.939 | 0.899 |

- **The decision** is CVu coverage vs CV-w vs CF at L2–L3:
  - **Expected (supports A4):** CVu ≈ CV-w at both rungs — ~0.94 at L2,
    **under-covering ~0.85–0.89 at L3** — and clearly below CF (~0.94–0.95) at
    L3. Then "internal CV ≠ cross-fitting" holds for the unweighted object the
    application reports, and the paper can cite CVu directly.
  - **Surprise (STOP-and-report):** CVu ~nominal at L3 (≥ CF) would mean the
    weighted nuisance fit, not the missing sample-split, drove the locked CV
    under-coverage — do NOT fold A4 into the paper without review.
  - Also inspect `se_ratio` (CVu should sit near CV-w's ~0.74–0.80 at L3 if the
    story is anti-conservative SEs from in-sample overfitting).

## Caveats / assumptions

- **Not rep-paired on fold splits.** Samples are seed-identical to the locked
  run, but the internal-CV fold split and learner RNG draw from a different
  stream (offset documented above), so CVu vs CV-w is a distribution-level
  comparison on identical data draws, not a per-rep paired one. This matches
  how all other arms are compared.
- CVu is undefined on single-learner rungs (L1_param, L4_aggressive) — same
  restriction as the engine's CV arm (`length(learners) > 1`).
- The locked CV-w rows came from `run_estimators()` where the CV arm shares the
  rep with 4 other arms; the worker RNG state at the CV fit therefore differed
  from this run's. Identical point estimates per rep are NOT expected; identical
  *samples* are.
- `gbound = 1e-3` in the EIF reconstruction (engine default) is intentionally
  NOT the 0.05 OOF floor — it mirrors the locked CV-w/CF rows exactly.
- Outputs are non-destructive: everything lands under
  `sim_output/arc_runs/R18_cvu/` and `results/arc/`; the locked
  `results/sim_full_summary.csv` is only ever READ.
