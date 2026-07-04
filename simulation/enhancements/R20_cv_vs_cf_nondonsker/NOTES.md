# R20_cv_vs_cf_nondonsker — internal CV vs cross-fitting at a non-Donsker library

**Type:** simulation | **Cost:** 20 array tasks × 32 cpus, ≤6 h each. L5 contains the deep RF, so per-rep cost ≈ the headline L4 rung plus the extra internal-CV arm. The standard-design cells (n≈1518) are lighter than R1.

## What & why (the reviewer gap this closes)

The manuscript ladder deliberately makes its two **non-Donsker** rungs **single-learner**
(L1 = `SL.glm`, L4 = `SL.ranger.deep`) — see `codes/learners.R` header: the aggressive rung uses a
*single* data-adaptive learner "so that SuperLearner's internal CV cannot hide overfitting behind a
smooth ensemble member." A side effect: the **Fully-Aware-CV** arm is gated `if (length(learners) > 1)`
in `codes/estimators.R`, so it is **structurally absent at L1 and L4**. Consequently the head-to-head a
reviewer asks for — *does internal cross-validation substitute for cross-fitting at a genuinely
non-Donsker library?* — is never shown. The CV arm only exists at L2/L3 (≈Donsker), where it already
undercovers at L3 (0.847 standard / 0.877 R1).

R20 runs **all five arms** (including the weighted internal-CV arm) at a **multi-learner, non-Donsker**
library:

```
L5_nondonsker = c("SL.glm", "SL.earth", "SL.ranger.deep")
```

= L3's shape with the **aggressive deep RF** swapped in for the tame `SL.ranger.t1`. Now the CV arm
fires (length > 1) at a library with genuine non-Donsker capacity, so we can compare
**Fully-Aware-CV vs Fully-Aware-CF vs Fully-Aware (single-fit)** coverage and se_ratio at one rung.

**Two legitimate outcomes, both reviewer-relevant** (the run decides which):
1. CV still **undercovers** (internal CV keeps weight on the overfitter / picks it in-sample) while CF
   holds → clean: *internal CV is not a substitute for cross-fitting.*
2. CV **recovers** coverage (its CV down-weights the deep RF toward the smooth members) while CF also
   holds → nuance: internal CV can mask the problem by retreating to smooth learners — but then you are
   not actually using the flexible learner, and CF lets you use it safely with a guarantee.

Either way it directly answers "what about CV at a non-Donsker library?" — aggregate.R prints the verdict.

## Engine reuse / non-destructive

SOURCES `codes/{config,dgp,estimators,diagnostics,learners}.R` read-only. `run_estimators()` is
UNCHANGED; the length-3 library makes it emit all five arms incl. the internal-CV arm. `L5` is defined
in `run.R` (NOT added to `SL_LADDER`), so no engine file is edited. Samples are seed-identical to the
locked headline run (`draw_sample(SAMPLE_SEED_BASE + global rep index)`) — only the library differs.

## Smoke test

```
# PowerShell (plumbing check: cheap glm+earth, 8 reps, ~2-4 min; validates the CV gate fires + all 5 arms emit)
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"; $env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"; $env:R20_DIR="<repo-root>/simulation/enhancements/R20_cv_vs_cf_nondonsker"
$env:R20_OUT="<repo-root>/sim_output/arc_runs/R20_cv_vs_cf_nondonsker"; $env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "$env:R20_DIR/run.R"
# bash: SMOKE=1 REPO_ROOT=$PWD SIM_CODE=$PWD/codes DATA_ROOT=$PWD/sim_output \
#   R20_DIR=$PWD/codes/arc_runs/R20_cv_vs_cf_nondonsker R20_OUT=$PWD/sim_output/arc_runs/R20_cv_vs_cf_nondonsker \
#   SLURM_CPUS_PER_TASK=2 Rscript codes/arc_runs/R20_cv_vs_cf_nondonsker/run.R
# Set R20_SMOKE_DEEP=1 to eyeball the REAL deep-RF L5 (slow, minutes/rep).
```
**SMOKE-GATE:** PASS iff all b/se finite, **5 arms present, the Fully-Aware-CV arm is present** (the
whole point — confirms the multi-learner CV gate fires), and the Fully-Aware-CF arm is present.

## Submit + aggregate

```
sbatch codes/arc_runs/R20_cv_vs_cf_nondonsker/submit.slurm          # FULL --array=1-20
# after the array completes:
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R20_OUT="$PWD/sim_output/arc_runs/R20_cv_vs_cf_nondonsker" \
  Rscript codes/arc_runs/R20_cv_vs_cf_nondonsker/aggregate.R
```

## Expected output

Per-task: `sim_output/arc_runs/R20_cv_vs_cf_nondonsker/r20_<scenario>_chunk###.rds` (+ manifest/).
Summary: `results/arc/R20_cv_vs_cf_nondonsker_summary.csv` (10 rows = 2 scenarios × 5 arms, standard
columns). aggregate.R prints the CV-vs-CF-vs-single-fit decision view per scenario + a verdict, and a
ladder-context view (CV at L3, CF at L4) from the locked `sim_full_summary.csv` if present.

**KEY NUMBERS to inspect:** Fully-Aware-CF coverage (expect ≈0.95, holds); Fully-Aware single-fit
coverage (expect a collapse, ≈0.1–0.3, matching L4); Fully-Aware-CV coverage (**the new datum** — does
it track CF or the single-fit?). se_ratio < 1 ⇒ anti-conservative.

## Caveats

- Smoke validates **plumbing** (cheap glm+earth); the deep-RF science path is exercised only on ARC
  (and is identical to the validated L4/R02/R16 deep-RF path). Set `R20_SMOKE_DEEP=1` to eyeball it.
- The CV arm here is the **weighted** internal-CV (the simulation foil, as in the locked CV arm). The
  application-matched **unweighted** CV-u is R18; if you want CV-u at L5 too, run R18 with its library
  overridden to L5 (a one-line env change documented in R18/NOTES.md) — not included here to keep the
  headline comparison clean.
- Local R 4.5.1 vs ARC R 4.4.0 — same standing cross-version note as every prior run.
