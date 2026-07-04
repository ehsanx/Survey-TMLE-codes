# R16_harmonized_sim — Harmonized-truncation headline simulation (spec A9a)

**Spec.** `Writing/comments/phase4-arc-sim-specs.md` item **A9a — Harmonized-truncation
headline simulation**: *"Remove the floor/cross-fitting confound from Figure 1 / Table 1:
apply the SAME propensity floor + Q-truncation across all five arms at every rung. ...
Output: re-run Figure 1 / Table 1 numbers with harmonized truncation. Targets: the
headline coverage table (audit §1) — re-derive if the harmonized version is adopted as
primary."*

**What.** The locked headline grid (`codes/run_sim.R` → `run_estimators`) lets the
truncation rule vary WITH the arm: the single-fit arms let `tmle()` bound its internal
SL propensity and reconstruct the EIF at `gbound=1e-3`, while the CV/CF arms floor their
internal-CV/out-of-fold propensity at `g_oof_bound=0.05` before targeting. The
single-fit-vs-CF coverage contrast therefore confounds **cross-fitting** with the
**floor**. R16 re-runs the FULL headline grid (2 scenarios × 4 rungs, 1000 reps) with
ONE rule for all FIVE paper arms:

- g predictions floored to `[g_floor, 1−g_floor]` (default **0.05**) BEFORE targeting,
- Q predictions truncated to `[q_lo, 1−q_lo]` (default **1e-3**),
- targeting always via `tmle(Y, A, W, family="binomial", obsWeights=<arm>, Q=, g1W=)`
  with PRE-FIT nuisances (the pattern validated by `R03_isolation_2x2`),
- EIF via `.eif_from_tmle(fit, Y, A, w_arm, gbound = g_floor)`,
- SEs via the engine's `.se_des()`.

| arm | nuisance fit | targeting weight | SE |
|---|---|---|---|
| `Fully-Aware-h` | weighted single-fit (`.sl(weights=w)`), in-sample | `w` | clustered design (Eq 8) |
| `Partially-Aware-h` | SAME fit/EIF as arm 1 | `w` | weighted, NO clustering |
| `Non-Aware-h` | unweighted single-fit | `rep(1,n)` | `sqrt(var(eif)/n)`, df=n−1 |
| `Fully-Aware-CV-h` | WEIGHTED internal-CV (surveyCV folds; sim-path foil); multi-learner rungs only | `w` | clustered design |
| `Fully-Aware-CF-h` | UNWEIGHTED per-fold OOF (`make_cf_folds`) | `w` (pooled) | clustered design |

`Fully-Aware-CF-h` is **identical to the locked CF arm** except the EIF `gbound`
1e-3 → `g_floor` (a no-op at the default 0.05 since the OOF g handed to `tmle` is
already in [0.05, 0.95]) — so any CF difference vs the locked table is pure MC/SL
randomness, a built-in sanity anchor. The arms that actually CHANGE are the single-fit
trio, whose nuisances are now pre-fit, floored at 0.05, and Q-truncated before targeting.

**Why.** The headline question: does the L4 single-fit undercoverage (locked coverage
0.220 standard / 0.216 R1) **survive** harmonization? Expected **YES** per
`R03_isolation_2x2` (SF-W at the harmonized 0.05 floor still collapsed at L4 while CF
held ~nominal). If yes, the Figure 1 / Table 1 contrast is NOT a truncation artifact
and the harmonized table can be reported (or adopted as primary per the spec's Targets
note). If no → STOP-AND-REPORT (see decision rule).

**Comparability (preserved).** The population is built ONCE per task from `POP_SEED`,
and each rep draws via `draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i)` with `i`
the GLOBAL rep index — **samples are byte-identical to the locked headline run for the
same scenario**, so harmonized-vs-locked differences are attributable to the truncation
rule alone, not sampling noise. Worker RNG streams use
`iseed = SAMPLE_SEED_BASE + 1000L*cell_index + task` (cell stride 1000 > max task 80,
so no (cell, chunk) collision); the stream only governs SL-internal CV splits / ranger
randomness, never the sample realizations.

---

## Files

| file | role |
|---|---|
| `estimators_harmonized.R` | all NEW logic: `run_harmonized()` = the five paper arms at one truncation rule; reuses `.sl`/`.eif_from_tmle`/`.se_des`/`make_cf_folds` from `codes/estimators.R` (read-only) |
| `run.R` | driver; one (scenario × rung × chunk) cell per array task; checkpointed; SMOKE gate |
| `aggregate.R` | per-task RDS → `results/arc/R16_harmonized_sim_summary.csv` + harmonized-vs-locked decision view |
| `submit.slurm` | SLURM array (mirrors R04 resources; `--array=1-80`, 04:00:00) |

`_inspect_smoke.R` is a throwaway local check (safe to delete).

Outputs per task: `r16_<scenario>_<rung>_chunk###.rds` (SMOKE runs prefix `SMOKE_`) +
`manifest/<SMOKE_>manifest_r16_<scenario>_<rung>_chunk###.rds` (params, truth, seeds incl.
POP_SEED/SAMPLE_SEED_BASE/TRUTH_SEED + the per-task RNG-stream iseed, pop audit,
R/package versions, timestamp, git sha, hostname).

---

## Smoke test (run locally before submitting)

Tiny single cell (standard × L1_param), 20 reps, 2 cores, `truth_M=2e5`. Finishes in
~2–5 min. PowerShell:

```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R16_DIR="<repo-root>/simulation/enhancements/R16_harmonized_sim"
$env:R16_OUT="<repo-root>/sim_output/arc_runs/R16_harmonized_sim"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R16_harmonized_sim/run.R"
```

bash:

```bash
# from the repo root
SMOKE=1 REPO_ROOT="$PWD" SIM_CODE="$PWD/codes" DATA_ROOT="$PWD/sim_output" \
  R16_DIR="$PWD/codes/arc_runs/R16_harmonized_sim" \
  R16_OUT="$PWD/sim_output/arc_runs/R16_harmonized_sim" \
  SLURM_CPUS_PER_TASK=2 \
  Rscript codes/arc_runs/R16_harmonized_sim/run.R
```

Writes `sim_output/arc_runs/R16_harmonized_sim/SMOKE_r16_standard_L1_param_chunk001.rds`
(+ a `SMOKE_manifest_...` under `manifest/`), prints the inline per-arm summary table and
the explicit `[SMOKE-GATE] PASS/STOP` verdict line.

### SMOKE-GATE DECISION RULE

- **PASS / submit the full run** iff ALL of:
  1. exactly the four arms `Fully-Aware-h`, `Partially-Aware-h`, `Non-Aware-h`,
     `Fully-Aware-CF-h` are present and `Fully-Aware-CV-h` is ABSENT (L1 is a
     single-learner rung — the CV arm must be skipped, like the engine);
  2. the UNWEIGHTED-nuisance arms (`Non-Aware-h`, `Fully-Aware-CF-h`) are
     strictly all-finite (R03 reference: SF-U/CF-U 1000/1000 finite), AND the
     WEIGHTED-nuisance FA-h/PA-h divergence share is **≤ 0.25** with every NA
     exactly matching a `diverged` flag. (The spec's original "all finite" gate
     assumed a clean GLM rung; the locked R03 full run already showed the
     weighted single-fit-nuisance + pooled-targeting path — which FA-h IS —
     diverging in 172/1000 L1 reps, `results/arc/R03_isolation_2x2_summary.csv`.
     The guard marks those reps NA; the aggregate counts them as `n_diverged`.)
  3. the FA-h and PA-h `psi` values are IDENTICAL per rep (they share ONE fit;
     only the SE differs) — `run.R` also hard-`stopifnot`s this invariant;
  4. sanity (not gated, but inspect): CF-h coverage at L1 should be high
     (locked CF L1 standard = 0.942; with 20 reps expect ≥ ~0.85, mcse ~0.08);
     FA-h on surviving reps should look like R03's SF-W L1 row (coverage ~0.95,
     `se_ratio` ~1.6 — conservative, never anti-conservative).
- **STOP** otherwise; the verdict line names the failed condition. Fix the code
  (never the engine) and re-run the smoke.

---

## Submit (full run on ARC)

```bash
# from the SLURM_SUBMIT_DIR (the repo root on scratch)
sbatch codes/arc_runs/R16_harmonized_sim/submit.slurm
```

- Grid = 2 scenarios × 4 rungs × ceil(1000/100)=10 chunks = **80 array tasks**
  (`--array=1-80`; TEST variant `--array=1,11,21,31,41,51,61,71` = chunk 1 of each cell).
- 32 cpus / 64 G / 04:00:00 per task (account `YOUR_ALLOCATION`), mirroring the R04
  exemplar. Re-submission is idempotent (per-chunk `arc_skip_if_done` checkpoint).
- Keep `--array`, `SIM_N_REPS=1000`, `SIM_CHUNK=100` consistent; `run.R` stops on a
  bad index.

### Aggregate after the run

```bash
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R16_OUT="$PWD/sim_output/arc_runs/R16_harmonized_sim" \
  Rscript codes/arc_runs/R16_harmonized_sim/aggregate.R
```

(Run it where `results/sim_full_summary.csv` exists — repo root locally or on scratch —
so the decision view can print the harmonized-vs-locked comparison.)

---

## Cost estimate

Per rep, `run_harmonized` does ~2× the engine's SL work on single-learner rungs
(4 single-fit `.sl` fits at V=2 + 10 CF fold-fits + 3 cheap targeting steps vs the
engine's 2 internal tmle SL fits + CF) and ~1.2× on multi-learner rungs. The locked
full run completed comfortably within 4 h/task at 32 cores with the same chunk size;
L4 (deep RF, single learner, V=2 internal CV) is the heaviest rung here but ranger on
n≈1.5–2k rows is seconds per fit. Expected wall time: **≤ ~1 h per L1/L4 task,
~1.5–3 h per L2/L3 task** (earth/gam ensemble at V=10) → the 04:00:00 cap is safe.
Total ≈ 80 tasks × ≲2 h × 32 cpus ≈ **≤ 5k core-hours** upper bound (likely ~half).

---

## Expected outputs + KEY NUMBERS to inspect

- `results/arc/R16_harmonized_sim_summary.csv` — columns
  `g_floor, q_lo, scenario, rung, method, n_reps, n_diverged, Psi, bias, emp_sd,
  mean_se, se_ratio, coverage, mcse_cov, deff_clust, icc_eif`.
- `R16_OUT/R16_harmonized_sim_combined.rds` — per-rep rows + summary.
- The aggregate prints two decision views:
  1. the harmonized summary table;
  2. **harmonized-vs-locked** (`results/sim_full_summary.csv`) coverage + se_ratio per
     (scenario, rung, base arm), then the HEADLINE verdict.
- **Key numbers:**
  - `Fully-Aware-h` coverage at **L4_aggressive** vs locked FA (0.220 standard /
    0.216 R1). HEADLINE: undercoverage **survives** harmonization if it stays `< 0.90`
    (expected: still collapsed, ~0.2–0.4 per R03's SF-W ≈ harmonized single-fit).
  - `Fully-Aware-CF-h` coverage at L4 vs locked CF (0.985 / 0.992) — should match
    within MC noise (built-in anchor; if it deviates > ~3 mcse, suspect a bug).
  - L1–L3 rows: the UNWEIGHTED-nuisance arms (NA-h, CF-h) should track the locked
    table closely. The WEIGHTED-nuisance arms (FA-h/PA-h, and CV-h at L2/L3) carry
    the R03-documented signature instead: at L1 expect `n_diverged` ≈ 150–200/1000
    (R03 SF-W: 172) with surviving-rep coverage ~0.95 and `se_ratio` ~1.6
    (conservative). This non-convergence of the pre-fit weighted GLM path is
    itself a harmonization finding — the locked FA (tmle-internal nuisances)
    never diverges — and must be reported alongside the table.
  - `n_diverged` at L4 should be ~0 (R03 SF-W at L4: 0/1000 — the deep RF does
    not separate; the divergence is a weighted-GLM-on-complex-covariates issue).

### FULL-RUN DECISION RULE (printed by aggregate.R)

- **CONFIRMS** (expected): `Fully-Aware-h` L4 coverage `< 0.90` in BOTH scenarios while
  `Fully-Aware-CF-h` holds ~nominal → the headline contrast is NOT a truncation
  artifact; report the harmonized table as the robustness companion (or adopt as
  primary per A9a's Targets note and re-derive the audit §1 numbers).
- **STOP-AND-REPORT**: `Fully-Aware-h` L4 coverage `≥ 0.90` in ANY scenario → the
  locked contrast is (partly) a floor artifact; reconcile with R03 before adopting
  either table as primary.

---

## Caveats / assumptions

- **What "harmonized" changes for the single-fit arms:** the locked FA/PA/NA arms let
  `tmle()` fit nuisances internally (its own `gbound` bounding) and reconstruct the EIF
  at 1e-3; R16 pre-fits them with `.sl`, floors g at 0.05 and truncates Q at 1e-3 BEFORE
  targeting. So harmonized FA-h ≠ locked FA numerically even at L1 — small shifts on
  smooth rungs are expected and are part of the answer, not a bug.
- `tmle()` applies its own internal `gbound` to a supplied `g1W`; at the default
  `g_floor=0.05` this never binds. If overriding `R16_GFLOOR` below ~0.025, tmle's
  internal bound (not ours) becomes the binding floor — keep `g_floor ≥ 0.025`.
  Overriding `R16_GFLOOR`/`R16_QLO` also requires a **fresh `R16_OUT`** (or deleting
  the old chunks): the per-chunk checkpoint cannot distinguish knob values, so a
  re-submit into the same dir would silently pool mixed-floor chunks (aggregate.R
  stops on a mixed-knob pool).
- The CV-h arm keeps the engine's WEIGHTED internal-CV nuisances (the simulation-path
  foil; the real-data path is unweighted) — emitted only on multi-learner rungs
  (L2_smooth, L3_adaptive), exactly like the engine, so locked CV rows exist for the
  same cells.
- Divergence guard (house pattern from R03/R05): reps whose fluctuation diverges
  (`max|eps| > 20`) are marked non-convergent (`b/se = NA`, counted as `n_diverged`)
  rather than polluting the 1000-rep means. KNOWN to trigger for the weighted-GLM
  single-fit path at L1 (~17% of reps; locked R03 SF-W = 172/1000, CF-W = 166/1000)
  — the weighted GLM on the Kang–Schafer complex covariates separates and the
  pooled weighted fluctuation blows up (eps → ~1e14). Verified in the local smoke
  (2/20 reps, eps ≈ 1e14, both flagged + NA'd; unweighted arms 20/20 finite).
  If the harmonized table is adopted as primary, the L1 FA-h `n_diverged` must be
  disclosed (e.g. a table footnote), since the locked FA arm does not have this
  failure mode.
- Non-destructive: everything lands under `R16_OUT` and `results/arc/`; the locked
  `results/sim_full_summary.csv` and `sim_output/intermediate/` are never touched.
- Local smoke validated on R 4.5.1 (Windows); ARC runs r/4.4.0 with the shared
  `R_LIBS` — same package set as the locked run (`SuperLearner, tmle, survey,
  surveyCV, earth, gam, glmnet, ranger`).
