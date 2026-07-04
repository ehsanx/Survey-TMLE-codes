# R21_deployable_cvcf — NOTES

> **⚠ ARC-ONLY, AUTHOR SUBMITS.** This is a >10 min/cell job (5-member
> CV-SuperLearner with tuned xgboost + deep RF, fit three ways on the TMLE side
> PLUS a full independent AIPW estimator). Local long-runs are duds and are
> forbidden by the standing rule. **Do NOT run the FULL sim locally.** Only the
> `SMOKE=1` plumbing gate (cheap glm+earth+glmnet, 8 reps) may run locally to
> validate wiring. The author submits to ARC.

Spec: `Writing/JASA/R21_L6_ARC_SPEC.md` (authoritative). Pre-authorized by
`Writing/JASA/REFRAME_PATHS.md:13` (more-realistic-aggressive learners).
Strengthens BOTH venues; decides nothing on the JASA-vs-Biometrics fork.

## What R21 settles (two gaps, one run)

1. **The lone-forest optic.** R20's L5 = `c(SL.glm, SL.earth, SL.ranger.deep)`
   is a 3-member library whose one aggressive member is a near-lone
   interpolating forest (the appendix calls it "the deployable counterpart of
   the lone-forest L4"). **L6 is the library one would ACTUALLY deploy**: a
   richer CV-SuperLearner with **tuned xgboost** carrying the aggressive load
   and the deep RF demoted to one of 5 members.
2. **The missing AIPW-CV cell (the C3 cell).** R15 has only AIPW-SF and AIPW-CF;
   the TMLE engine has a cluster-aware Fully-Aware-CV arm but its AIPW
   counterpart was never built. **R21 adds AIPW-CV** so "cluster-CV ≠ CF" is
   tested **estimator-agnostically** on a deployable class. Per the recon this
   is the single unobserved cell that **gates A1's framing**: if AIPW-CV stays
   nominal where TMLE-CV under-covers, the CV pathology localizes to TMLE and
   A1 weakens (the aggregator prints this verdict explicitly).

## Library L6 (defined in run.R, NOT in SL_LADDER)

```
L6 = c('SL.glm','SL.earth','SL.glmnet','SL.xgboost.tuned','SL.ranger.deep')
```
- `SL.xgboost.tuned` (defined in `r21_helpers.R`): `ntrees=500, shrinkage=0.05,
  minobspernode=5`, **`max_depth=6` fixed by default** (cheap + reproducible).
  Aggressive-but-DEPLOYABLE — not interpolating-by-construction (cf. the lone
  deep RF). Set `R21_XGB_TUNE=1` to genuinely tune depth over {3,6,8}, or
  `R21_XGB_MAXDEPTH=N` to fix a different depth.
- `SL.ranger.deep` (kept but **demoted** to one of 5): the rate-breaking member
  stays in the library, so L6 is unambiguously deployable yet still contains an
  interpolating learner. (Open author call at submission: keep-but-demote
  [recommended, implemented] vs replace with tuned-xgb as sole aggressive load.)
- `SL.hal9001` (OPTIONAL, behind `R21_WITH_HAL=1`, **OFF by default**): slow /
  can fail to install on some ARC images; never make the headline depend on it.

## Arm matrix — 6 per cell

| | single-fit | cluster-CV | CF |
|---|---|---|---|
| **TMLE** | Fully-Aware | Fully-Aware-CV | Fully-Aware-CF |
| **AIPW** | AIPW-SF | **AIPW-CV (new)** | AIPW-CF |

- TMLE arms: emitted by the **unchanged** `run_estimators()` (the CV arm fires
  because `length(L6) > 1`).
- AIPW arms: `aipw_arms_all()` in `r21_helpers.R` = R15's `aipw_arms()` (SF+CF,
  unchanged) + the new `aipw_cv_arm()`.

### AIPW-CV implementation + ⚠ WEIGHTING-CONVENTION GUARD

`aipw_cv_arm()` mirrors `estimators.R:178-198` on the AIPW side:
- **cluster-aware folds:** one weighted `SuperLearner(Q,g)` call whose
  `cvControl = list(validRows = split(seq_len(n), surveyCV::folds.svy(d,
  nfolds=5, clusterID='cluster')))` — whole-PSU internal-CV folds, exactly the
  engine recipe (with the same `tryCatch -> default folds` guard).
- **in-sample predict** for Q1/Q0/g (engine does `predict(qg, Xa1)` etc.).
- **clamps copied VERBATIM from R15** (`aipw_helpers.R:148-154`):
  `g -> [0.05, 0.95]`, `Q -> [1e-3, 1-1e-3]`.
- then the **Hájek AIPW point + JKn replicate SE (`.aipw_jkn`) + Eq-8 linearized
  SE (`.se_des`)** — the `one_arm()` body from `aipw_helpers.R:157-172`, raw
  weights throughout.

**The guard (critical):** the engine's CV arm is **WEIGHTED in the simulation
path** "as a deliberate foil" (`estimators.R:186-187`:
`w_cv <- if (is.null(inpop)) w else rep(1,n)`). So `aipw_cv_arm()` fits its CV
nuisances **WEIGHTED**, using R15's **mean-1 weight normalization** (`w/mean(w)`)
exactly as AIPW-SF does (raw survey weights spuriously diverge glm's IRLS;
SuperLearner/tmle apply the same internal normalization). The Hájek point and
**both** survey variances keep the **RAW** weights, identical to AIPW-SF. Net:
across the three AIPW arms only the protocol axis (single-fit vs cluster-CV-folds
vs cross-fit) differs — **not** the weighting — so the AIPW SF/CV/CF comparison
is on the CV-vs-CF axis, not confounded by the weighting protocol. (AIPW-CF
stays UNWEIGHTED out-of-fold, matching the primary TMLE-CF arm; that asymmetry
is intrinsic to the CF protocol and is identical in the TMLE engine, so the
TMLE-vs-AIPW comparison is apples-to-apples within each design column.)

## Hypotheses (vs locked targets)

Locked R20/L5 (TMLE): standard CF **0.948** (se_ratio 1.066) / single-fit
**0.905** / cluster-CV **0.850**; R1 CF 0.934 / SF 0.909 / CV 0.874.
Locked R15/L4 (AIPW): AIPW-SF **0.621** (se_ratio ~0.43) / AIPW-CF **0.96**.

- TMLE-single-fit ≈ 0.90; **TMLE-cluster-CV < 0.90** (CV does not rescue);
  TMLE-CF ≈ 0.93–0.95 (conservative, se_ratio > 1).
- AIPW-single-fit under-covers; **AIPW-cluster-CV ALSO under-covers** (the new
  estimator-agnostic point); AIPW-CF ≈ 0.96.
- **Confirming → cluster-CV ≠ CF is estimator-agnostic on a deployed class.**
  If internal CV instead recovers coverage on L6, the headline FLIPS (reportable
  either way, not a hidden failure: CV retreats to smooth members, so you lose
  the flexible learner; CF lets you keep it with a guarantee). If AIPW-CV is
  nominal but TMLE-CV under-covers, the CV pathology localizes to TMLE and A1
  weakens — also flagged by the aggregator.

## Grid / seeding

- Scenarios `{standard, R1}`; 1 rung (L6); `n_reps=1000`, `SIM_CHUNK=100` →
  10 chunks/cell → **20 tasks** (`--array=1-20`). tasks 1-10 standard, 11-20 R1.
- m baseline via `draw_sample` scenario defaults: standard `base_m=6` (m≈60,
  n≈1500), R1 `base_m=2` (m≈100, n≈2000).
- **OPTIONAL m-sweep** behind `R21_MSWEEP=1` (`base_m {6,20}`) → 4 cells →
  `--array=1-40`. **Gated OFF** so the core 20-task run ships first.
- Seeding reuses R20/R15 EXACTLY: per-rep `draw_sample(sample_seed =
  SAMPLE_SEED_BASE + i)` (samples **byte-identical** to the locked headline /
  R15 / R20 runs); parallel stream `iseed = SAMPLE_SEED_BASE + 1000*cell_index +
  task`. Sources `_checkpoint.R`; calls `arc_skip_if_done` early (idempotent
  resubmission of the same array+chunking).

## SMOKE gate

```
SMOKE=1 Rscript codes/arc_runs/R21_deployable_cvcf/run.R
```
Standard scenario, cheap `{glm,earth,glmnet}`, 8 reps. **Asserts all 6 arms —
incl. BOTH `*-CV` arms (TMLE-CV and the new AIPW-CV) — emit finite b/se** and
all six method labels are present. Prints `[SMOKE-GATE] PASS` or `STOP` with the
per-check reasons. (Plumbing only; the deep-RF / xgboost science path is
validated by R15/R20.)

## Env flags

| flag | default | effect |
|---|---|---|
| `R21_WITH_HAL` | off | append `SL.hal9001` to L6 (slow; optional) |
| `R21_MSWEEP` | off | base_m sweep {6,20} → doubles the array to 40 |
| `R21_XGB_TUNE` | off | CV-tune xgboost depth over {3,6,8} (default: fixed 6) |
| `R21_XGB_MAXDEPTH` | 6 | fixed xgboost depth when not tuning |

## Run / aggregate (ARC; AUTHOR submits)

```bash
# repo-checkout path:
sbatch codes/arc_runs/R21_deployable_cvcf/submit.slurm        # --array=1-20, 12h
# after it completes:
Rscript codes/arc_runs/R21_deployable_cvcf/aggregate.R        # -> results/arc/R21_deployable_cvcf_summary.csv
```
WALLTIME: 12 h is a conservative cap (R21 > R20's 6 h: 5-member library, fit 3
ways on the TMLE side + a full extra AIPW estimator). PROBE one task first
(`sbatch --array=1 --time=12:00:00 ...`), read `sacct -j <id> --format=Elapsed`,
then run the rest at ≈3× that. Walltime-recovery: resubmit ONLY the failed task
ids at the **unchanged** SIM_CHUNK with a longer `--time` (the checkpoint keys on
chunk number, not rep range — never change SIM_CHUNK against existing files).

The self-contained ARC bundle entry (`<ARC-bundle-dir>`) wires
all paths automatically: `bash submit_all.sh R21` then `bash aggregate_all.sh R21`.

## Aggregate output

`aggregate.R` writes `results/arc/R21_deployable_cvcf_summary.csv` and prints the
headline **3×2 DECISION VIEW** `{single-fit, cluster-CV, CF} × {TMLE, AIPW}` with
**BOTH coverage AND se_ratio** (se_ratio mandatory: the honest claim is
ONE-SIDED — "single-fit/cluster-CV under-cover; CF does NOT under-cover" — and CF
is *conservative* (se_ratio > 1), not exactly nominal). Per-cell verdict lines
cover the four outcomes: estimator-agnostic CV≠CF (strengthens A1), headline-flip
(CV recovers), CV-pathology-localizes-to-TMLE (A1 weakens), or mixed.

## Provenance / standing rules

Mirrors R20's package structure; sources `config/dgp/estimators/diagnostics/
learners.R` + R15's `aipw_helpers.R` **read-only** — **no engine file edited**.
L6, `SL.xgboost.tuned`, and the AIPW-CV arm are defined locally. On completion:
freeze into the internal `survey-tmle-reproduce` archive alongside the locked
summaries; re-run `verify_inline_numbers.R` and extend it for any new headline
numbers before they enter the manuscript. **Do NOT push/commit/submit — author
signs off.**
