# R07_informative — Informative-sampling-beyond-C stress test

**Purpose.** Defend the de-weighted out-of-fold (OOF) nuisance choice in the
`Fully-Aware-CF` arm. That arm fits per-fold nuisances **unweighted**, justified
by the *sampling-ignorability* assumption **S ⊥ (A, Y) | C** (selection into the
sample is independent of treatment/outcome given the confounders C). A referee
will ask: *what if sampling is informative beyond C?* This run builds a DGP in
which sampling is **informative beyond C**, sweeps an informativeness knob
`rho ∈ {0, 0.3, 0.6, 0.9}`, and compares two cross-fitted arms that differ in
**exactly one line** — whether the OOF fits use the survey weights:

| arm | OOF nuisance fits | consistent when |
|---|---|---|
| `Fully-Aware-CF-unwt` | **unweighted** (paper default) | S ⊥ (A,Y) \| C holds |
| `Fully-Aware-CF-wt`   | **weighted**   (defended alt.) | also under informative-beyond-C |

Both arms share everything downstream: g flooring at `g_oof_bound=0.05`, the
pooled **weighted** `tmle()` targeting, and the Eq-8 clustered design SE. So any
**bias difference** is attributable *solely* to de-weighting the OOF fits. We
report the **bias and coverage curves** for both arms vs `rho`.

This is the most novel DGP work in the resubmission package; the selection
mechanism is documented in full below.

---

## The selection mechanism (math)

The canonical population (`codes/dgp.R::make_population`) draws **three
independent** PSU-level random intercepts, constant within a PSU:

- `uC` — PSU effect in the **confounders** C (shared by C1, C2).
- `uY` — PSU effect in the **outcome** logit (integrated out of the estimand).
- `uA` — PSU effect in the **treatment** logit.

Crucially `uY, uA ⊥ uC`, so the confounders C carry **no information** about the
outcome/treatment PSU effects. The canonical `draw_sample` selects PSUs by
**SRS within strata**, with the only informativeness coming through `b_h` (a
fixed stratum grid that *does* enter C) — hence canonical sampling is
*non-informative given C*: weighting is not needed for consistency, which is
what licenses the unweighted OOF fit.

We break this. For each PSU `g` we build a latent **selection driver**,
standardized across PSUs to mean 0 / unit variance:

```
r_g    = z( uY_g + uA_g )                     # "outcome/treatment" part, z = standardize
eps_g  ⟂ {uY, uA, uC, b_h}, Var=1             # idiosyncratic noise, residualized vs ALL PSU effects
uS_g   = rho * r_g + sqrt(1 - rho^2) * eps_g  # Var(uS_g) = 1,  Corr(uS_g, r_g) = rho  (exactly)
```

The `eps` noise is the **OLS residual of N(0,1) noise regressed on
`[1, uY, uA, uC, b_h]`** — i.e. projected onto the orthogonal complement of the
span of *every* PSU-level driver of the data law, then standardized. This is the
key negative-control construction (built with `.lm.fit`, NOT `qr.solve`, which
returns a pivoted basic solution rather than the least-squares fit). Two
properties follow, verified on the realized population:

- **`rho = 0` is a clean negative control.** `cor(uS, uY) = cor(uS, uA) =
  cor(uS, uC) ≈ 0` to machine precision, so selection is on a variable
  independent of `(A, Y, C)`. The unweighted in-sample `E[Y]` does **not** drift
  (≈ +0.001), and *both* arms are unbiased — de-weighting is provably harmless.
  (Residualizing only against `r = z(uY+uA)`, or leaving `uS` correlated with
  the confounder effect `uC`, both reintroduce spurious rho=0 bias — that is why
  the full 5-vector projection is used.)
- **Informativeness is genuinely BEYOND C.** As `rho → 1`, `uS` rotates toward
  `r = z(uY+uA)`; since `uY, uA ⟂ uC` in the population, `cor(uS, uC)` stays ≈ 0
  for **all** rho while `cor(uS, uY), cor(uS, uA)` grow. So C never sees the
  selection — the unweighted-OOF (C-only) fit cannot recover, exactly the regime
  the referee asks about. Verified sweep (lambda = 0.8):

  | rho | cor(uS,uY) | cor(uS,uA) | cor(uS,uC) | unwt E[Y] drift |
  |----:|-----------:|-----------:|-----------:|----------------:|
  | 0.0 |  0.00 |  0.00 |  0.00 | +0.001 |
  | 0.3 |  0.24 |  0.19 | −0.02 | +0.04  |
  | 0.6 |  0.48 |  0.37 | −0.04 | +0.09  |
  | 0.9 |  0.72 |  0.56 | −0.06 | +0.13  |

**Stage-1 PSU selection** (within each stratum `h`): size measure

```
s_g = exp(lambda * uS_g)                       # lambda = informative-selection strength
```

Select `m_h` PSUs **without replacement** by **sequential / systematic PPS**
(Madow) with first-order inclusion probabilities

```
pi1_g = min(1,  m_h * s_g / sum_{g' in h} s_{g'}),   normalized so sum_g pi1_g = m_h
```

(units whose proportional `pi` exceeds 1 are taken with certainty and the
remaining sample size is reallocated — standard Tillé inclusion-probability
capping). **Stage-2** is the canonical SRS-WOR of `n0` units within a chosen PSU
(unchanged), `pi2 = n0 / M`. The **Horvitz–Thompson weight** is

```
w_i = 1 / (pi1_g * pi2),     so   E[ sum_{i in sample} w_i ] = N_pop,
```

i.e. weights are design-unbiased for population totals, so the **Hájek mean**
the estimators target stays consistent for `Psi` **by construction** — the
weighted arms remain valid; only the *unweighted-OOF* fit is put at risk.

**Two knobs, two roles:**

- `lambda` (fixed at **0.8**) sets the **strength** of informative selection =
  how unequal the weights are for a given `uS`. Tuned so `mean(sum w / N) ≈ 1`
  and the selection visibly perturbs the in-sample (A,Y) law as rho grows
  (`w_cv` rises from ≈0.16 at rho=0 to ≈0.9 at rho=0.9 — a realistic survey
  weight spread). Not pathological (no near-zero `pi1`).
- `rho` (the **swept** knob) sets *what* the selection is on. At `rho = 0` it is
  pure noise independent of `{uY, uA, uC}` → ignorable given C (indeed
  unconditionally) → **both arms unbiased**. As `rho → 1`, `uS → r_g` →
  selection is driven by the very `uY/uA` effects that C cannot see →
  **S ⊥ (A,Y) | C is violated** → the de-weighted OOF fit loses consistency
  while the weighted OOF fit stays closer to truth.

### On the weight spread vs rho
The `eps` part of `uS` is the residual after projecting noise off the 5-vector
`{1,uY,uA,uC,b_h}`, so at `rho=0` its variance (hence the `exp(lambda*uS)` size
spread) is modest → `w_cv ≈ 0.16`. As rho grows, `uS` loads on the full-variance
`r=z(uY+uA)` and the weights spread out (`w_cv ↑`). Both effects are *intended*:
at every rho the weights are **unequal** (selection bites), but only beyond
`rho=0` does that selection couple to (A,Y) — isolating the **cost of
de-weighting** from the mere presence of unequal weights. The rho=0 cell is a
genuine negative control (clean weights on a C-independent nuisance), not a
no-weighting cell.

---

## Files in this folder

| file | role |
|---|---|
| `dgp_infsamp.R` | **new DGP logic.** `attach_uS(pop, rho, seed)` builds the per-PSU driver; `draw_sample_infsamp(pop, sample_seed, rho, lambda, uS_tab)` does PPS-WOR stage-1 + canonical stage-2 and returns the same `obs` contract as `draw_sample`. Reuses `apply_L` from `codes/dgp.R`. |
| `cf_arms.R` | the **two CF arms.** `run_cf_pair(obs, learners)` runs `Fully-Aware-CF-unwt` and `Fully-Aware-CF-wt`, reusing the building blocks exposed by `codes/estimators.R`: `make_cf_folds`, `.sl`, `.eif_from_tmle`, `.se_des`. The only difference between arms is `weights = NULL` vs `weights = w[tr]` in the OOF `.sl()` calls. |
| `run.R` | the **driver.** Sources the canonical engine + the two helpers. One SLURM array task = one (rho × rung × rep-chunk). Builds population + truth + `uS(rho)` once, runs its rep-chunk in parallel, writes one RDS + a manifest. |
| `aggregate.R` | combines per-task RDS into `results/arc/R07_informative_summary.csv` (matching `aggregate_sim.R` columns) and prints the **decision-rule readout**. |
| `submit.slurm` | SLURM array job (mirrors `codes/sim.slurm` resources exactly). |
| `_dgp_check.R` | local de-risking: cheap cross-fitted **GLM-TMLE** ATE-bias sweep vs rho (no SuperLearner) — confirms both arms unbiased at rho=0 and the bias/selection moves with rho. Safe to delete. |
| `_qcheck.R` | instant DGP sanity: prints `cor(uS, uY/uA/uC)` and unweighted `E[Y]` drift vs rho (no estimator fits). Safe to delete. |

**Non-destructive:** nothing in `codes/*.R` is edited. Outputs go to
`sim_output/arc_runs/R07_informative/` (per-task RDS + manifests) and
`results/arc/R07_informative_summary.csv` — none of the locked results are
touched.

---

## How to smoke-test (locally, ~2–4 min)

Runs 20 reps for `rho ∈ {0, 0.9}` at `L2_smooth` only (one population build).

**PowerShell (Windows, R 4.5):**
```powershell
$env:SMOKE="1"
$env:SIM_CODE="<repo-root>\codes"
$env:DATA_ROOT="<data-root>\sim_output"
$env:R07_DIR="<repo-root>\codes\arc_runs\R07_informative"
$env:R07_OUT="<data-root>\sim_output\arc_runs\R07_informative"
$env:SLURM_ARRAY_TASK_ID="1"; $env:SLURM_CPUS_PER_TASK="4"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "$env:R07_DIR\run.R"
```
Then run the next rho cell by setting `SLURM_ARRAY_TASK_ID=2` (SMOKE grid is
rho∈{0,0.9} × {L2_smooth} × 1 chunk = 2 tasks), and aggregate:
```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "$env:R07_DIR\aggregate.R"
```

**Bash/ARC login node (R 4.4):**
```bash
export SMOKE=1 SIM_CODE=$PWD/codes DATA_ROOT=$PWD/sim_output \
       R07_DIR=$PWD/codes/arc_runs/R07_informative \
       R07_OUT=$PWD/sim_output/arc_runs/R07_informative \
       SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=4
Rscript codes/arc_runs/R07_informative/run.R
```

**Even cheaper DGP-only sanity** (no SuperLearner, ~3–5 min, validates the
mechanism + decision rule with a GLM-TMLE):
```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>\codes\arc_runs\R07_informative\_dgp_check.R"
```

---

## How to submit (full run)

From the repo root (= SLURM submit dir / scratch):
```bash
sbatch codes/arc_runs/R07_informative/submit.slurm
```
Grid = rho{0,0.3,0.6,0.9} (4) × rung{L2_smooth,L3_adaptive} (2) × 10 rep-chunks
= **80 array tasks** (`--array=1-80`), 1000 reps per (rho × rung).
After all tasks finish:
```bash
SIM_CODE=$PWD/codes R07_OUT=$PWD/sim_output/arc_runs/R07_informative \
  Rscript codes/arc_runs/R07_informative/aggregate.R
```

---

## Expected output

- **Per task:** `sim_output/arc_runs/R07_informative/r07_rho<NNN>_<rung>_chunk<CCC>.rds`
  (e.g. `r07_rho090_L3_adaptive_chunk004.rds`) with `$results` (per-rep
  `method,b,se,df,deff,icc_eif,...` for both CF arms) and `$diagnostics`
  (`eps_cf,g_cf_min,g_cf_max,...`); a manifest under `.../manifest/`.
- **Aggregated:** `results/arc/R07_informative_summary.csv` with columns
  `rho, rung, method, n_reps, Psi, bias, emp_sd, mean_se, se_ratio, coverage,
  mcse_cov, deff_clust, icc_eif, w_cv, df_design` — the **bias and coverage
  curves vs rho** for both arms. Plus `R07_informative_combined.rds`
  (per-rep + summary).

**Key numbers to inspect** (the bias/coverage curve per rung):
- `bias` for `Fully-Aware-CF-unwt` should be ≈0 at `rho=0` and **grow in
  magnitude** as `rho → 0.9`; its `coverage` should **fall** below nominal.
- `bias` for `Fully-Aware-CF-wt` should stay **closer to 0** across rho and hold
  `coverage` near 0.95 — quantifying the cost of de-weighting and the regime
  where the unweighted default is safe.

**Validated reference (cheap GLM-TMLE de-risk, 20 reps, `_dgp_check.R`,
lambda=0.8).** A *correctly-specified GLM* nuisance is fairly weight-robust, so
the de-weighting cost is modest at this tier — but the negative control and the
"selection bites" property are unambiguous, which is what the smoke gate checks:

| rho | UNWT bias | WT bias | unwt E[Y] drift | w_cv |
|----:|----------:|--------:|----------------:|-----:|
| 0.0 | −0.007 | +0.001 | −0.003 | 0.20 |
| 0.3 | +0.005 | +0.004 | +0.038 | 0.31 |
| 0.6 | +0.007 | +0.007 | +0.067 | 0.53 |
| 0.9 | +0.010 | +0.010 | +0.098 | 0.78 |

The **full run uses SuperLearner** (`L2_smooth`, `L3_adaptive`) where flexible
OOF learners interact with the weights — that is where the *separation* between
`-unwt` and `-wt` is expected to widen at high rho (the GLM tier above
deliberately understates it). The smoke's job is to confirm the DGP is sound
(clean rho=0, monotone drift), not to reproduce the final separation.

---

## Decision rule (smoke-gate: STOP-and-report vs proceed)

`aggregate.R` prints a **selection-bites readout**: for each rung it reports
`|bias|` at `rho=0` and `rho=0.9` for both arms.

**PROCEED to the full run iff**, at the smoke/DGP-check stage:

1. **`rho=0` is a clean negative control** — `cor(uS, uY+uA) ≈ 0` (printed by
   `run.R`; the `_qcheck.R` helper additionally confirms `cor(uS, uY) =
   cor(uS, uA) = cor(uS, uC) ≈ 0` and unweighted `E[Y]` drift ≈ 0.001) and
   **both** arms' `|bias|` at `rho=0` are small (within ~2× the `Fully-Aware-CF`
   bias seen in the locked `sim_full_summary.csv`, ≲ 0.02 for L2/L3 standard),
   and `mean(sum w / N) ≈ 1`. In particular the **WEIGHTED** arm must be unbiased
   at rho=0 — if it is not, the weights are pathological or `uS` is not clean:
   **STOP** (this exact failure was caught and fixed in development; see the
   `qr.solve` vs `.lm.fit` note in `dgp_infsamp.R`).
2. **The selection genuinely bites** — the **WEIGHTED** arm's `|bias|` (or, in
   the cheap GLM check, the unweighted in-sample `E[Y]` drift) **moves with
   rho**. If the weighted-arm bias and the in-sample (A,Y) law are *flat* across
   rho, the DGP is mis-implemented (selection not coupling to uY/uA): **STOP and
   report** — a flat unweighted curve would then be uninformative, not a defense.
3. **Honest direction** — the **UNWEIGHTED** arm's `|bias|` is **non-decreasing**
   in rho (de-weighting costs something as informativeness grows). If the
   unweighted bias *shrinks* with rho, something is inverted in the size measure
   or weights: **STOP and investigate** before burning the full 80-task budget.

If the smoke shows (1) clean rho=0, (2) the weighted/in-sample law clearly
shifting with rho, and (3) unweighted bias growing with rho, **proceed**. If the
weighted curve is flat (selection not biting) or the unweighted curve is flat or
inverted, **STOP-and-report** rather than submit — the result would not support
the de-weighting defense.

---

## Dependencies

- Canonical engine: `codes/config.R`, `dgp.R`, `estimators.R`, `diagnostics.R`,
  `learners.R` (sourced, never edited).
- R packages: `SuperLearner, tmle, survey, surveyCV, earth, gam, ranger`
  (same set as `run_sim.R`; R 4.4 on ARC / 4.5 locally). The `L2_smooth` rung
  needs `gam, earth`; `L3_adaptive` adds `ranger`.
- Seeds reused from `config.R`: `POP_SEED=20260606`, `SAMPLE_SEED_BASE=1000`,
  `RNGkind("L'Ecuyer-CMRG")`. `uS(rho)` is deterministic in `(POP_SEED, rho)`.

## Caveats

- `uS` is built from `pop$pop$uY, uA` — these columns are emitted by the
  canonical `make_population` (verified). If a future engine change drops them,
  `attach_uS` will error loudly (no silent fallback).
- `lambda=0.8` is a tuned default for the **standard** scenario. Re-tune if you
  switch scenarios (R1 has only 2 PSUs/stratum → PPS reduces to near-SRS and the
  informative channel weakens; the spec fixes this run to **standard**).
- We deliberately **balance** stage-2 `n0` and `m_h` across strata (no `b_h`
  oversampling) so the *only* informativeness is the new `uS` channel — keeping
  the "beyond C" interpretation clean. This differs from canonical `draw_sample`
  (which oversamples high-`b_h` strata); that is intentional and documented.
- Systematic PPS can occasionally return one too few/many PSUs from rounding; the
  sampler tops-up/trims deterministically to exactly `m_h` (guarded in code).
