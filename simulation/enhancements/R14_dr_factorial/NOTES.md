# R14_dr_factorial — Double-robustness factorial (incl. the informative-sampling unweighted-nuisance threat)

**Spec item.** A6 in `Writing/comments/phase4-arc-sim-specs.md` ("Double-robustness factorial
(incl. informative-sampling unweighted-nuisance case)"). The paper claims double robustness but
only ever misspecifies BOTH nuisances jointly (the headline `model_type="complex"` transforms feed
both Q and g). This run crosses the misspecification **per nuisance** and adds the sharpest
theory-facing threat from the external review (Gemini's mechanism): an **unweighted** misspecified
out-of-fold library converges to the **sample** L2 projection, which under **informative** sampling
differs from the population projection — so the unweighted-OOF CF default could lose double
robustness exactly when sampling is informative.

## Design (locked)

- **Covariate device.** Samples are drawn with `model_type = "simple"`, so `obs$L1..L4` carry the
  RAW `C1..C4` (`apply_L` is the identity for `"simple"` — verified in `codes/dgp.R` lines 83–85).
  `dr_helpers.R` builds two design frames per sample: **C** (correct) = the raw columns; **W**
  (wrong) = `apply_L(L1..L4, "complex")` (the Kang–Schafer transforms, the headline run's canonical
  misspecification device). Per cell, the Q-learner sees `q_spec`'s frame (+A) and the g-learner
  sees `g_spec`'s frame, **independently**.
- **Population.** `make_population("standard", model_type = "simple", truth_M = 2e6L)`. The
  population data and the truth Ψ are **identical across model_type** (R01 established this;
  `model_type` only changes what the learners SEE via `apply_L` at sampling time).
- **Arms per rep** (3, in `dr_helpers.R`), all with the harmonized 0.05 propensity floor (nuisance
  preds AND `.eif_from_tmle(gbound = 0.05)`), Q clamped at 1e-3, weighted pooled `tmle(Q=, g1W=)`
  targeting (the R03 pattern), clustered design SE (Eq 8):
  - **FA-w** — weighted single-fit nuisances via `.sl(weights = w)` on the spec'd frames;
  - **CF-u** — UNWEIGHTED per-fold OOF fits (the paper default) over `make_cf_folds`;
  - **CF-w** — WEIGHTED per-fold OOF fits on the **same** folds (paired, as in R07's `cf_arms.R`),
    so the CF-u − CF-w bias gap is attributable solely to de-weighting the OOF fits.
  R03's divergence guard is carried over: reps where the targeting fluctuation exceeds
  `max|epsilon| > 20` are marked non-convergent (`b/se = NA`; counted as `n_diverged`).
- **Cells.** `q_spec ∈ {C,W} × g_spec ∈ {C,W} × sampling ∈ {info (oversample=TRUE), noninfo
  (oversample=FALSE)} × rung ∈ {L1_param, L2_smooth}` = **16 cells**. FULL = 1000 reps in chunks of
  100 → `--array=1-160`. All cells are cheap (GLM/smooth; no deep RF) → `--time=02:00:00`.
- **Comparability.** `draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, ...)` with `i` the
  GLOBAL rep index, and `draw_sample` is internally `.scoped_seed`-ed — so for the **info** cells
  the sampled rows are **byte-identical to the locked headline run's rep `i`** (same `standard`
  population; `model_type` only changes the L columns, not the selection). The noninfo cells differ
  only through `oversample = FALSE`.
- **RNG streams.** `clusterSetRNGStream(iseed = SAMPLE_SEED_BASE + 1000L*cell_index + task)`:
  `cell_index ∈ 1..16` gives offsets ≥ 1000 apart while `task ≤ 160`, so no two (cell, chunk)
  tasks ever share a stream. The stream only drives fold assignment + SL CV splits; the sample
  draw is governed by `sample_seed`.

## Interpretation notes (for the write-up)

- At **L1 (GLM)** the wrong spec is sharply Kang–Schafer-misspecified — the decisive DR read.
  At **L2 (smooth glm+gam+earth)** the library **partially recovers** the transforms (gam/earth can
  partially undo monotone warps), so W-cell biases shrink relative to L1 — that contrast is part of
  the story, not a bug.
- With `oversample = FALSE` all stratum sampling fractions are equal, so weights are constant
  (`w_cv ~ 0`): the design is **non-informative by construction**, and CF-u vs CF-w can differ only
  by noise — the clean control for the threat view.

## Files

| file | role |
|---|---|
| `run.R` | driver; one (cell × chunk) per array task; checkpointed; SMOKE gate |
| `dr_helpers.R` | all NEW logic: the two design frames + the three DR arms |
| `aggregate.R` | per-(q,g,sampling,rung,method) summary + the two decision views |
| `submit.slurm` | SLURM array (mirrors R04 resources; 02:00:00) |

`_inspect_smoke.R` / `_debug_rep20.R` are throwaway local checks (safe to delete).

Engine files (`codes/*.R`) are sourced read-only; nothing outside this folder,
`sim_output/arc_runs/R14_dr_factorial/`, and `results/arc/` is written.

## Cost estimate

- L1 cells: GLM-only `.sl` (internal V=2) — ~5–15 s/rep ⇒ a 100-rep chunk on 32 cores ≈ 2–5 min.
- L2 cells: 3-learner smooth library, 22 `.sl` calls/rep (2 single-fit + 20 fold fits) — ~1–3
  min/rep ⇒ a 100-rep chunk on 32 cores ≈ 10–25 min.
- FULL run: 160 tasks, each ≪ 2 h; whole array typically clears in 1–2 wall-hours given queue
  parallelism. Total ≈ 40–80 core-hours — small by R02/R04 standards.

## Smoke test (run locally before submitting)

One cheap both-correct cell (q=C, g=C, info, L1_param), 20 reps, 2 cores, `truth_M = 2e5` —
~2–3 min. Outputs get a `SMOKE_` filename prefix and are excluded from aggregation.

PowerShell (Windows dev box, R 4.5.1):
```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:SIM_CODE="<repo-root>/R"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:R14_DIR="<repo-root>/simulation/enhancements/R14_dr_factorial"
$env:R14_OUT="<repo-root>/sim_output/arc_runs/R14_dr_factorial"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R14_dr_factorial/run.R"
```

bash (ARC login node or any POSIX shell, from the repo root):
```bash
SMOKE=1 \
  REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" SIM_CODE="$PWD/codes" \
  R14_DIR="$PWD/codes/arc_runs/R14_dr_factorial" \
  R14_OUT="$PWD/sim_output/arc_runs/R14_dr_factorial" \
  SLURM_CPUS_PER_TASK=2 \
  Rscript codes/arc_runs/R14_dr_factorial/run.R
```

### SMOKE-GATE DECISION RULE

The smoke prints a per-arm table and an explicit `[SMOKE-GATE] PASS/STOP: reason` line.
- **PASS ⇒ submit the full array** iff ALL of:
  1. all three arms (FA-w, CF-u, CF-w) produced finite summary stats; **CF-u (the unweighted
     paper default) has `n_diverged = 0`** — any CF-u divergence is a real bug; the WEIGHTED arms
     only need a convergent share ≥ 70% (weighted GLM nuisance fits are *known* to diverge or
     collapse under these heavy-tailed weights — the locked R03 standard/L1 rows show SF-W
     172/1000 and CF-W 166/1000 diverged, and base-R `glm()` itself fails to converge on those
     draws — so zero-divergence would be a stricter bar than the validated engine behavior);
  2. both-correct bias small: `|bias| <= 2*emp_sd/sqrt(n_convergent)` for every arm (this is the
     q=C, g=C cell — any real bias here means the frames or the targeting are wired wrong);
  3. CF folds sane: a single `V_eff` in [2,5] (expect **5**: the standard design has 6
     PSUs/stratum and `V_cf = 5`).
- **STOP** otherwise: fix the run folder (never the engine) and re-smoke. A bias gate failure with
  a *correct* spec is a wiring bug, not a finding.
- Validated locally 2026-06-12 (R 4.5.1, 20 reps, ~1.5 min): FA-w 18/20 convergent
  (bias +0.0031), CF-u 20/20 (bias −0.0040), CF-w 19/20 (bias −0.0062), all
  `|bias| <= 2*emp_sd/sqrt(n)`, `V_eff = 5` ⇒ `[SMOKE-GATE] PASS`.

## Submit (full run on ARC)

```bash
# from the SLURM_SUBMIT_DIR (the repo root on scratch)
sbatch codes/arc_runs/R14_dr_factorial/submit.slurm
```
- FULL `--array=1-160` (16 cells × 10 chunks); TEST `--array=1-160:10` (chunk 1 of every cell).
- Keep `--array`, `SIM_N_REPS=1000`, `SIM_CHUNK=100` consistent; `run.R` stops on a bad index.
- Re-submission is idempotent: finished chunks skip in seconds (`arc_skip_if_done`).

### Aggregate after the run
```bash
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R14_OUT="$PWD/sim_output/arc_runs/R14_dr_factorial" \
  Rscript codes/arc_runs/R14_dr_factorial/aggregate.R
```

## Expected outputs + KEY NUMBERS to inspect

- **Per task:** `sim_output/arc_runs/R14_dr_factorial/r14_q<C|W>_g<C|W>_<info|noninfo>_<rung>_chunk###.rds`
  (per-rep `method,b,se,df,diverged` + design checks + diagnostics) and
  `manifest/manifest_*.rds` (params, truth, seeds incl. `iseed`/`cell_index`, package versions,
  git sha, pop audit).
- **Aggregated:** `results/arc/R14_dr_factorial_summary.csv` with leading knob columns
  `q_spec,g_spec,sampling,rung,method` then the standard
  `n_reps,n_diverged,Psi,bias,emp_sd,mean_se,se_ratio,coverage,mcse_cov,deff_clust,icc_eif,w_cv,df_design`;
  combined RDS under `R14_OUT`.
- **KEY NUMBERS:**
  1. **VIEW 1 (DR 2×2, L1_param × noninfo):** bias in (C,C), (C,W), (W,C) all ≈ 0 (within
     ~3×`emp_sd/sqrt(n)` ≈ 0.001) with coverage ≈ 0.95; (W,W) bias clearly nonzero (order of the
     headline complex-DGP L1 bias, ~+0.01–0.02) with degraded coverage. That is the DR table for
     Web-D.
  2. **VIEW 2 (THREAT, g_spec=W cells):** primary `gap` = mean per-rep **paired** difference
     `b(CF-u) − b(CF-w)` over reps where both CF arms converged (paired MCSE; `n_pairs` and CF-w
     `n_diverged` are printed on each line; the unpaired cell-mean gap follows as a secondary
     line). The Gemini mechanism predicts `|gap_info| > 0` (materially, > 2×paired mcse) while
     `|gap_noninfo| ≈ 0`. The aggregate prints a `FLAG` per (rung, q_spec) cell and an overall
     verdict line. Either outcome is publishable: FLAG ⇒ state the unweighted-OOF limitation under
     informative sampling; no FLAG ⇒ robustness result for the paper default.
  3. `w_cv` ≈ 1.4 in info cells (the headline standard-design weight spread) and ≈ 0 in noninfo
     cells (design-informativeness sanity check). `n_diverged` MUST be 0 for CF-u everywhere;
     for the WEIGHTED arms expect ~10–20% diverged in the **info** L1 cells (cf. the locked R03
     standard/L1 rows: SF-W 172/1000, CF-W 166/1000) and ~0 in the noninfo cells (equal weights
     cannot destabilize the fits).

## Caveats / assumptions

- The "wrong" frame is the **code-canonical** Kang–Schafer form (`apply_L` "complex" in
  `codes/dgp.R`, the numerically sane variant) — same device as the headline run, so the W-cells
  here are directly comparable to the locked complex-DGP results.
- `tmle()` is called with `Q=` and `g1W=` pre-supplied, so its `W` argument is inert for nuisance
  fitting; we pass the raw frame for determinism (engine pattern, `estimators.R` CV/CF arms).
- The threat view's **primary** gap is the mean of **per-rep paired differences**
  `b(CF-u) − b(CF-w)`, joined on rep within each cell and restricted to reps where **both** CF arms
  converged, with paired MCSE = `sd(diff)/sqrt(n_pairs)` (the arms share the sample and the folds,
  so pairing is exact). The unpaired cell-mean gap is printed as a **secondary** reference line
  only: it has two defects — (i) its independence MCSE over-states the noise (merely conservative),
  but (ii) it also conditions each arm on its **own** convergent subset (CF-w drops ~17% of
  informative L1 draws per the locked R03 analogue while CF-u keeps ~all), and since weighted-fit
  divergence is concentrated in exactly the informative cells the FLAG keys on, that selection
  asymmetry could **manufacture or mask** an informative-only gap — it is *not* conservative.
  The paired primary removes (ii) by conditioning both arms on the same reps.
- Sampling here is informative through the stratum device (`oversample=TRUE` oversamples high-`b_h`
  strata, shifting the sample C-distribution); it is NOT the R07 beyond-C selection (`uS`). A6's
  threat is about the sample-vs-population projection of a *misspecified* model, for which the
  stratum device suffices; R07 covers the beyond-C channel.
- **Weighted-fit instability (diagnosed, expected).** Under the info design the survey weights are
  heavy-tailed (range ≈ 22–500); on some draws base-R `glm()` itself fails to converge for the
  weighted nuisance fit and collapses fitted probabilities to ~0 (reproduced locally on rep 20's
  g-fit: `glm.fit: algorithm did not converge`, fitted ≡ 2.2e-16). This is the SAME phenomenon
  behind the locked R03 divergence counts. The R03 guard (max|epsilon| > 20 → `b/se = NA`,
  counted in `n_diverged`) catches the diverging fluctuations; the residual collapse reps that do
  NOT diverge inflate `mean_se`/`se_ratio` for FA-w/CF-w exactly as in the locked R03 rows
  (SF-W `se_ratio` 1.64, CF-W 1.87) — read VIEW 2 on **bias**, which is robust to this, and the
  per-task RDS keeps per-rep `g_min/g_max` so the collapse share can be audited later.
- `deff_clust` / `icc_eif` in the summary are means of the FA-w-EIF diagnostics over **all** reps
  **including diverged ones** (the EIF is computed before the divergence guard nulls `b/se`). In
  informative cells — where the weighted fits diverge at the R03-documented rates — they should
  **not** be quoted as clean design-effect estimates; treat them as diagnostics only.
- The DR statement is about **bias**; FA-w single-fit coverage can deviate for the usual
  single-fit reasons even in unbiased cells — read coverage against the CF arms.
