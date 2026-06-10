# R04_nuisance_rate — Nuisance-rate diagnostic ladder

**Purpose.** Convert the one "not-fixable" theory item — *no weighted-clustered
oracle inequality for the survey-weighted SuperLearner L2 rate* — into a
**DEMONSTRATED** fact. For each ladder rung L1–L4 and each scenario
(`standard`, `R1`), per rep we fit the **cross-fit** nuisances exactly as the
`Fully-Aware-CF` arm does (unweighted out-of-fold fits, OOF propensity floored at
`g_oof_bound=0.05`), then measure the **weighted L2 errors** of the OOF nuisance
estimates against the **population truth**:

```
eQ  = || Qhat - Q0 ||_{L2(w)}      (combined over arms a in {0,1})
eg  = || ghat - g0 ||_{L2(w)}
prod         = eQ * eg                         <- the TMLE remainder rate object
prod_sqrtm   = prod * sqrt(m)   with m = #sampled PSUs   <- the DEMONSTRATION
```

**Expected result.** `prod_sqrtm` **shrinks across L1 -> L3** (the product rate
`||Qhat-Q0||·||ghat-g0|| = o_P(m^{-1/2})` is met, so root-m-consistent
inference is justified) and is **deliberately violated at L4** (interpolating
deep RF, `SL.ranger.deep`, min.node.size=1, mtry=p), where it stays large /
grows — the empirical signature of the failed rate that the locked coverage
table already attributes to L4. This makes the "no oracle inequality" caveat a
*shown* fact rather than an assumed one.

This is instrumentation only: it reuses the canonical engine read-only and adds
~1 helper file. No `codes/*.R` is edited.

---

## How Q0 / g0 are extracted from the population (the crux)

Read from `codes/dgp.R::make_population()`. The realized finite population
`population$pop` carries, **per unit**, the true nuisance values used to generate
the data:

| population column | meaning | used as |
|---|---|---|
| `pscore = plogis(alpha0 + alpha_g*Csum + uA_i)` | true propensity, **incl. realized treatment PSU effect `uA_i`** | secondary g0 (`.g0_real`) |
| `Q1_real/Q0_real = plogis(gamma0 + theta*a + gamma_C*Csum + uY_i)` | outcome mean **incl. realized outcome PSU effect `uY_i`** | secondary Q0 (`.Q*_real`) |
| `Q1_int/Q0_int = gh_expit(gamma0 + theta*a + gamma_C*Csum, sigma2_Y)` | **uY-integrated** E[Y\|A=a,C=c] | **PRIMARY Q0 target** (`.Q*_int`, already on `obs`) |

**Which truth is the right target?** The cross-fit SuperLearner sees only the
covariates `C` (never the latent random effects `uY`, `uA`), so it is consistent
for the **integrated** nuisances, not the realized-`u` ones:

- **Primary Q0** = the **uY-integrated** conditional mean `Q*_int` — this is the
  estimand the TMLE actually targets (see `dgp.R` RC-1 comment), and the norm in
  which the product-rate condition is stated. `draw_sample()` already propagates
  `.Q1_int/.Q0_int` onto `obs`, so no re-keying is needed for the primary Q.
- **Primary g0** = the **uA-integrated** marginal propensity `E[A|C]`, recomputed
  with the canonical `gh_expit(alpha0 + alpha_g*Csum, sigma2_A)` (helper
  `.g0_marginal`). The g-learner also sees only `C`, so this — not `pscore` — is
  its consistency target.

**Assumption / re-keying.** To get `Csum` (needed for the marginal g0) and the
realized-`u` truths, `attach_truth()` joins each sampled row back to its
population row on the signature `(cluster, psu_within, A, Y, round(.Q1_int),
round(.Q0_int))`. This is unique per population row in practice; the join is
**asserted** and, if any collision/missing match is detected, the code falls back
to the **u-integrated-only** path (primary Q0/g0 still measured; realized-`u`
columns set NA). The smoke run reports `truth_join=OK` and full-run manifests
record it per task. We also report the **secondary** realized-`u` errors
(`*_real`, against `Q*_real` incl. `uY` and `pscore` incl. `uA`) for completeness,
since the spec names "incl. random effects" — but the rate claim is made on the
primary (integrated) columns.

---

## Files

| file | role |
|---|---|
| `run.R` | driver; sources canonical engine read-only, runs one (scenario×rung×chunk) cell |
| `nuisance_rate_helpers.R` | all NEW logic: truth extraction (`attach_truth`, `.g0_marginal`), CF OOF refit (`cf_nuisance_oof`, reusing `.sl`/`make_cf_folds`), weighted L2 (`.wL2`), per-rep (`one_rep_rate`) |
| `aggregate.R` | combine per-task RDS -> per-rung×scenario rate table |
| `submit.slurm` | SLURM array (mirrors `codes/sim.slurm` resources) |

`_test_join.R` / `_test_l4.R` are throwaway local checks (safe to delete).

---

## Smoke test (run locally before submitting)

Tiny 1-cell (standard × L1_param), 20 reps, serial — finishes in ~1–2 min. It
prints the per-rung rate row and the truth-join status.

```bash
# from the repo root (<repo-root>)
SMOKE=1 \
  REPO_ROOT="<repo-root>" \
  DATA_ROOT="<repo-root>/sim_output" \
  SIM_CODE="<repo-root>/codes" \
  R04_DIR="<repo-root>/codes/arc_runs/R04_nuisance_rate" \
  Rscript codes/arc_runs/R04_nuisance_rate/run.R
```

PowerShell equivalent:
```powershell
$env:SMOKE=1; $env:REPO_ROOT="<repo-root>"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:SIM_CODE="<repo-root>/codes"
$env:R04_DIR="<repo-root>/codes/arc_runs/R04_nuisance_rate"
Rscript codes/arc_runs/R04_nuisance_rate/run.R
```

Writes `sim_output/arc_runs/R04_nuisance_rate/SMOKE_r04_standard_L1_param_chunk001.rds`.
Validated locally (R 4.5.1): `truth_join=OK`, all columns finite,
`mean_prod_sqrtm ≈ 0.064` for L1 (the smooth/parametric baseline).

---

## Submit (full run on ARC)

```bash
# from the SLURM_SUBMIT_DIR (the repo root on scratch)
sbatch codes/arc_runs/R04_nuisance_rate/submit.slurm
```

- Grid = 2 scenarios × 4 rungs × ceil(200/50)=4 chunks = **32 array tasks**
  (`--array=1-32`). 32 cpus, 64 G, 4 h per task (L4 deep-RF is the heaviest rung;
  200 reps is ample — this is a mean-error diagnostic, not a coverage table).
- Resources mirror `codes/sim.slurm` exactly (account, cpus, mem); only
  array/time and the driver path differ.
- Keep `--array`, `SIM_N_REPS=200`, `SIM_CHUNK=50` consistent; `run.R` stops on a
  bad index.

### Aggregate after the run
```bash
SIM_CODE="$PWD/codes" \
  REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R04_OUT="$PWD/sim_output/arc_runs/R04_nuisance_rate" \
  Rscript codes/arc_runs/R04_nuisance_rate/aggregate.R
```

---

## Expected output

- **Per task:** `sim_output/arc_runs/R04_nuisance_rate/r04_<scenario>_<rung>_chunk<NNN>.rds`
  (one `per_rep` data.frame: `eQ_int, eg_int, prod_int, prod_int_sqrtm,
  prod_int_sqrtn`, plus realized-`u` `*_real` columns, `m_psu`, `n`, `V_eff`,
  `truth_join`) + a `manifest_*.rds` (params, truth, seeds, sessionInfo,
  package versions, git rev — like `run_sim.R`).
- **Aggregated:** `results/arc/R04_nuisance_rate_summary.csv` — the per-rung ×
  scenario table. **Key columns to inspect:** `mean_eQ`, `mean_eg`,
  `mean_prod`, `mean_prod_sqrtm` (the demonstration), and
  `prod_sqrtm_vs_L1` (ratio vs the L1 baseline within scenario).

The headline plot/claim: `mean_prod_sqrtm` (y) vs rung L1→L4 (x), one line per
scenario — monotone decrease L1→L3, sharp jump at L4.

---

## Decision rule (smoke gate + interpretation)

**Smoke gate (STOP-and-report vs proceed):**
- **STOP** and report if the smoke run shows `truth_join=FALLBACK_uint_only`
  (the population join failed — the marginal-g and realized-`u` columns would be
  NA; investigate the row signature before burning the full run), OR if any of
  `mean_eQ`, `mean_eg`, `mean_prod_sqrtm` is non-finite (NA/Inf), OR if the L1
  `mean_prod_sqrtm` is absurd (e.g. > 1, since L1 is a correctly-specified GLM
  and should give a small product). The validated smoke gives `truth_join=OK`
  and `mean_prod_sqrtm ≈ 0.064`, so the gate is **PASS** as shipped.
- **PROCEED** to the full 32-task run otherwise.

**Interpreting the full-run table (the scientific claim):**
- **Clean success:** `mean_prod_sqrtm` decreases L1 → L2 → L3 and L4 is clearly
  the largest (ideally ≥ ~2× L3, mirroring the L4 coverage collapse to ~0.17 in
  `results/sim_full_summary.csv`). Report directly: the product rate is met on
  L1–L3 and violated at L4.
- **L3 borderline:** if L3's `mean_prod_sqrtm` is not strictly below L2 (flat or
  slightly above) but L4 still blows up, frame it as the **STRONGER** claim:
  *cross-fitting buffers performance beyond the strict L2 rate* — i.e. CF
  coverage (already ~0.95 at L3 in the locked table) holds even when the raw
  product rate is only borderline, because the cross-fit empirical-process term
  is controlled by sample-splitting rather than by a Donsker/rate condition.
  This is a more favorable story for the method, not a failure.
- **Unexpected:** if L4's `mean_prod_sqrtm` is NOT the largest, STOP and report —
  this would contradict the locked coverage result and indicates a bug in the
  truth target or the OOF refit, not a finding.

---

## Dependencies

- Canonical engine (read-only): `codes/{config,dgp,estimators,diagnostics,learners}.R`.
- R packages (same as the main sim): `SuperLearner, tmle, survey, surveyCV, earth,
  gam, glmnet, ranger`. On ARC use `module load gcc/9.4.0 r/4.4.0` +
  `R_LIBS=$HOME/R/...`. Locally validated on R 4.5.1.
- Reuses canonical seeds (`POP_SEED`, `SAMPLE_SEED_BASE`) and RNG
  (`L'Ecuyer-CMRG`) — population + samples are identical to the locked sim, so
  the measured errors are on the same draws.

## Caveats

- **m for the rate.** We use `m = #sampled PSUs` (the independent design unit) for
  `prod_sqrtm`; `prod_sqrtn` (n = sample size) is reported alongside for
  reference. The clustered design makes #PSUs the right effective sample size.
- **g target.** The primary g0 is the **uA-integrated** marginal propensity
  `E[A|C]`, recomputed by Gauss–Hermite — NOT the realized `pscore`. If a reviewer
  insists on `pscore`, the `*_real` columns already carry that comparison.
- **Outputs are non-destructive:** everything lands under
  `sim_output/arc_runs/R04_nuisance_rate/` and `results/arc/`; the locked
  `results/sim_full_summary.csv` and `sim_output/intermediate/` are never touched.
