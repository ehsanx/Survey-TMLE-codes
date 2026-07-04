# R22_aipw_ladder_complete — NOTES

> **⚠ ARC-ONLY, AUTHOR SUBMITS.** >10 min/cell (the L5 deep-RF AIPW arms). Only
> `SMOKE=1` (8 reps, glm+earth) runs locally. **Run guide: `README.md`** (copy /
> run / download). **Spec / sign-off: `Writing/comments/R22_AIPW_COMPLETE_SPEC.md`.**

## Scope — the MINIMAL "α" completion (author decision 2026-06-14)

Computes ONLY the genuinely-missing AIPW cells of the unified L1–L7 figure; no
TMLE, no duplicate AIPW-SF/CF, no xgboost/HAL:

| rung | library | arms computed |
|---|---|---|
| `L2_smooth` | glm+gam+earth | **AIPW-CV** |
| `L3_adaptive` | glm+earth+ranger.t1 | **AIPW-CV** |
| `L5_nondonsker` | glm+earth+ranger.deep (= R20) | **AIPW-SF, AIPW-CV, AIPW-CF** |

`L1/L4 internal-CV` = **N/A by construction** (a lone learner has nothing to
weight → internal CV ≡ single-fit). After R22 the figure is gap-free except
those two principled N/A cells. **No locked number changes; audit only extends.**

Every output cell IS a deliverable (we recompute nothing locked) → no within-run
cross-check is possible. Trust = **code identity** with the already-validated AIPW
machinery (`R15/aipw_helpers.R::aipw_arms`, `R21/r21_helpers.R::aipw_cv_arm`),
applied to the L2/L3/L5 libraries on the SAME seeds. `aggregate.R` prints the new
cells next to their LOCKED neighbours (TMLE-CV@L2/L3, AIPW@L4, AIPW@L6) for a
plausibility eyeball.

## Grid / seeding

`{standard, R1} × {L2, L3, L5}` = 6 cells × 1000 reps (10×100-rep chunks) =
**60 tasks** (`--array=1-60`). Cell order (scenario fastest within rung): 1–10
std-L2 · 11–20 R1-L2 · 21–30 std-L3 · 31–40 R1-L3 · **41–50 std-L5 · 51–60 R1-L5**
(L5 = the expensive block; probe task 41). Same seeds as every locked run
(`draw_sample(sample_seed = SAMPLE_SEED_BASE + i)` → byte-identical samples);
stream `iseed = SAMPLE_SEED_BASE + 1000*cell_index + task`. Sources `_checkpoint.R`;
idempotent resubmission.

## SMOKE gate (local, done)

`SMOKE=1 Rscript codes/arc_runs/R22_aipw_ladder_complete/run.R` → standard,
glm+earth, AIPW "all", 8 reps; asserts AIPW-SF/CV/CF emit finite b/se.
**Result 2026-06-14: `[SMOKE-GATE] PASS`.**

## Provenance / standing rules

Mirrors R21's structure; sources `config/dgp/estimators/diagnostics/learners.R` +
R15's `aipw_helpers.R` + R21's `r21_helpers.R` **read-only** — no engine/R15/R21
file edited. On completion: extend `verify_inline_numbers.R` for any new headline
number; freeze into `survey-tmle-reproduce`. **Do NOT push/commit/submit — author
signs off.**
