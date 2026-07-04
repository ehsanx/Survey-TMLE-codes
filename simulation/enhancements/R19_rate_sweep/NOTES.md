# R19_rate_sweep — √m-scaled product-error diagnostic across the F.5 m-sweep

**Spec item.** `Writing/comments/phase4-arc-sim-specs.md`, section *"√m-scaled
product-error diagnostic across the F.5 m-sweep"* — the re-review's sharpest
evidence ask: show the realized cross-fit nuisance product error scaled by √m
**across** the large-m sweep, to probe whether the assumed product rate (C1)

```
||Qhat - Q0|| * ||ghat - g0||  =  o_P( m^{-1/2} )
```

is plausibly met at each ladder rung. A single-m diagnostic cannot test a
*rate*. The spec says "reuses existing sweeps", but in fact **R04 recorded
nuisance errors only at the default `base_m=6`** and **R02 swept `base_m`
without recording nuisance errors** — so this is a new (small) run that
**crosses** them: R04's truth-extraction/L2 machinery × R02's `base_m` knob.

**Expected result.** Per rung, the trend
`mean_prod_sqrtm(base_m) / mean_prod_sqrtm(base_m=6)` across
`base_m ∈ {6,12,20,30}` (`m_total ∈ {60,120,200,300}`) is **flat or declining
at L2–L3** (rate plausibly met), **growing at L4** (the interpolating deep
RF — rate fails), matching the L4 coverage collapse already in the locked
table, and **plausibly growing at L1**: `SL.glm` is *misspecified* for both
nuisances under the Kang–Schafer complex covariates, so eQ/eg have positive
approximation floors and `sqrt(m)*prod` must eventually grow. That matches
the locked R02 standard/L1 rows — flat Fully-Aware bias floor
0.0183/0.0169/0.0172/0.0188 with Fully-Aware-CF coverage decaying
0.949→0.912→0.888→0.848 across `m_total` 60→300 — and the level gap between
R04's locked L1 `mean_prod_sqrtm = 0.0625` at `m_total=60` and the R19 smoke
L1 value 0.0903 at `m_total=300` (predicted L1 `trend_vs_m6(m30)` ≈ 1.4). A
growing L1 trend is a legitimate, useful additional demonstration for Web-D,
NOT a bug. **Targets:** Web-D rate diagnostic + the A1/`rem:oracle` evidence
paragraph.

---

## Design

- **Cells:** scenario `standard` (FIXED) × `base_m ∈ {6,12,20,30}` ×
  rung L1–L4 = **16 cells**. FULL = 200 reps/cell, `SIM_CHUNK=50` → 4 chunks
  → **`--array=1-64`**, `--time=04:00:00`.
- **m convention (matches R04):** `m` = **total** sampled PSUs,
  `m_total = H*base_m = {60,120,200,300}` (standard `H=10`; `draw_sample`
  keeps `m_h` balanced across strata so `m_psu == m_total` exactly). Both
  `base_m` and `m_total` are recorded per row so the numbers **join the
  existing R02/R04 tables** on `(scenario, rung, base_m)`.
- **Reuse, not copy:** `run.R` `source()`s
  `codes/arc_runs/R04_nuisance_rate/nuisance_rate_helpers.R` from **R04's own
  folder** (env `R04_DIR`, default under `SIM_CODE`), getting `attach_truth`,
  `.g0_marginal`, `cf_nuisance_oof`, `.wL2`. R04's `one_rep_rate()` does NOT
  accept `base_m` (its draw step is hardcoded), so `rate_sweep_helpers.R` adds
  the thin wrapper `one_rep_rate_m()` which replicates **only the draw step**
  (`draw_sample(..., base_m=base_m)`, exactly like R02's `one_rep_m`) and then
  re-states R04's post-draw measurement arithmetic verbatim
  (provenance-marked; R04 does not expose it as a standalone function).
- **Comparability property (preserved):** reps are seeded by the GLOBAL rep
  index, `sample_seed = SAMPLE_SEED_BASE + i`. For a given `base_m` the drawn
  samples are **byte-identical to R02's m-sweep draws**, and at `base_m=6`
  byte-identical to the **locked headline run** and to **R04's draws** — so
  the measured errors are on the same samples as the published coverage rows.
- **RNG streams:** `clusterSetRNGStream(iseed = SAMPLE_SEED_BASE +
  1000L*cell_index + task)` with `cell_index ∈ 1..16`, `task ∈ 1..64 < 1000`
  → no two (cell, chunk) tasks share a stream. (Worker RNG drives only fold
  splits / SL internal CV; the sample draws use scoped per-rep seeds.)
- **Outputs:** per-task `r19_m<NN>_<rung>_chunk###.rds` + manifest RDS in
  `manifest/`. Per-chunk checkpoint `arc_skip_if_done` (full-array resubmits
  are idempotent). SMOKE outputs get a `SMOKE_` filename prefix.
- **Non-destructive:** no existing file is edited; everything lands under
  `sim_output/arc_runs/R19_rate_sweep/` and `results/arc/`.

## Files

| file | role |
|---|---|
| `run.R` | driver; one (base_m × rung × chunk) task; sources engine + R04 helpers read-only |
| `rate_sweep_helpers.R` | the ONLY new logic: `one_rep_rate_m()` (base_m-forwarding wrapper) |
| `aggregate.R` | per-(base_m, rung) rate table + the headline per-rung trend view |
| `submit.slurm` | SLURM array (mirrors R04 resources; `--array=1-64`, 4 h) |

---

## Cost estimate

Per rep = one CF nuisance refit (5 PSU-folds × 2 nuisances; **no** TMLE arms),
i.e. cheaper than a full sim rep. L1/L2 are trivial–cheap; L4 (deep RF,
`min.node.size=1`) at `base_m=30` (n = 7590) is the heavy corner:
~3–8 min/rep single-threaded → 50 reps on 32 cores ≈ 2 waves ≈ 10–30 min per
chunk. All 64 tasks fit easily in the 4 h walltime; total ≈ 64 × (a few
minutes to ~0.5 h) ≈ **well under one R02-sized job** (~30–60 node-hours at
32 cpus, most of it in tasks 49–64).

## Smoke test (run locally before submitting)

ONE cheap cell: `standard × L1_param × base_m=30`, 10 reps, **serial**,
`truth_M=2e5` (~2–4 min). Validates the base_m forwarding at the largest m
(n = 7590, `m_psu = 300`), the R04-helper sourcing, the truth join, and the
output/manifest paths.

PowerShell:
```powershell
$env:SMOKE="1"; $env:REPO_ROOT="<repo-root>"
$env:DATA_ROOT="<repo-root>/sim_output"
$env:SIM_CODE="<repo-root>/R"
$env:R19_DIR="<repo-root>/simulation/enhancements/R19_rate_sweep"
$env:R19_OUT="<repo-root>/sim_output/arc_runs/R19_rate_sweep"
$env:R04_DIR="<repo-root>/simulation/enhancements/R04_nuisance_rate"
$env:SLURM_CPUS_PER_TASK="2"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" "<repo-root>/simulation/enhancements/R19_rate_sweep/run.R"
```

bash:
```bash
# from the repo root
SMOKE=1 \
  REPO_ROOT="<repo-root>" \
  DATA_ROOT="<repo-root>/sim_output" \
  SIM_CODE="<repo-root>/R" \
  R19_DIR="<repo-root>/simulation/enhancements/R19_rate_sweep" \
  R19_OUT="<repo-root>/sim_output/arc_runs/R19_rate_sweep" \
  R04_DIR="<repo-root>/simulation/enhancements/R04_nuisance_rate" \
  Rscript codes/arc_runs/R19_rate_sweep/run.R
```

Writes `sim_output/arc_runs/R19_rate_sweep/SMOKE_r19_m30_L1_param_chunk001.rds`
(+ `manifest/SMOKE_manifest_m30_L1_param_chunk001.rds`) and prints an inline
summary plus an explicit `[SMOKE-GATE] PASS/STOP` line.

### SMOKE-GATE DECISION RULE

- **STOP** (do NOT submit) if any of:
  - `truth_join` ≠ `OK` on any rep (`FALLBACK_uint_only` means the population
    row join failed — the marginal-g/realized-u columns would be NA;
    investigate the row signature first);
  - any error norm non-finite — the primaries `eQ_int`, `eg_int`, `prod_int`,
    `prod_int_sqrtm` or the realized-u secondaries `eQ_real`, `eg_real`,
    `prod_real`, `prod_real_sqrtm` (run.R's `fin_ok` checks all eight);
  - L1 `mean_prod_sqrtm ≥ 1` (L1 is a plain GLM — misspecified under the
    complex covariates, but with a *bounded* approximation floor, so its
    scaled product stays small over this sweep: R04's locked value is 0.0625
    at `base_m=6` and the smoke measured 0.0903 at `base_m=30`; ≥ 1 would
    mean a broken truth join / OOF refit, not the floor).
- **PROCEED** to the full 64-task array otherwise. `run.R` prints the verdict
  itself as `[SMOKE-GATE] PASS: ...` / `[SMOKE-GATE] STOP: ...`.

---

## Submit + aggregate (ARC)

```bash
# from the SLURM_SUBMIT_DIR (repo root on scratch)
sbatch codes/arc_runs/R19_rate_sweep/submit.slurm

# after all 64 tasks finish:
SIM_CODE="$PWD/codes" REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" \
  R19_OUT="$PWD/sim_output/arc_runs/R19_rate_sweep" \
  Rscript codes/arc_runs/R19_rate_sweep/aggregate.R
```

## Expected outputs + KEY NUMBERS

- **Per task:** `sim_output/arc_runs/R19_rate_sweep/r19_m<NN>_<rung>_chunk###.rds`
  (per-rep `eQ_int, eg_int, prod_int, prod_int_sqrtm, prod_int_sqrtn`, the
  `*_real` secondaries, `base_m, m_total, m_psu, n, V_eff, truth_join`) +
  `manifest/manifest_m<NN>_<rung>_chunk###.rds`.
- **Aggregated:** `results/arc/R19_rate_sweep_summary.csv` +
  `<R19_OUT>/R19_rate_sweep_combined.rds`. Inspect:
  - **view 1** per (base_m, rung): `mean_eQ`, `mean_eg`, `mean_prod`,
    `mean_prod_sqrtm`, `sd_prod_sqrtm`, `prod_sqrtm_vs_L1` (vs L1 at the same
    base_m), `n_reps`, `truth_join`;
  - **view 2 (HEADLINE):** per rung, `trend_vs_m6` =
    `mean_prod_sqrtm(base_m)/mean_prod_sqrtm(base_m=6)` across the sweep.
- **Key numbers:** L1 row at `base_m=6` should reproduce R04's standard/L1
  `mean_prod_sqrtm = 0.0625` (identical samples + identical measurement code);
  `truth_join=OK` everywhere; `mean_m_psu == m_total` exactly. **Decision
  numbers:** `trend_vs_m6` ≤ ~1 (flat/declining) at L2–L3 → C1 plausibly met
  there; `trend_vs_m6` increasing in m at L4 → rate fails (the interpolating
  deep RF); `trend_vs_m6` increasing at L1 → the misspecified parametric
  rung's approximation floor showing through (expected — the smoke run
  already measured L1 `mean_prod_sqrtm = 0.0903` at `m_total=300` vs R04's
  locked 0.0625 at `m_total=60`, i.e. predicted L1 `trend_vs_m6(m30)` ≈ 1.4).

**Interpreting the full-run table (the manuscript claim):**
- *Clean success:* L2–L3 trends flat/declining AND L4 trend growing → report
  directly: the assumed product rate (C1) is consistent with the data on
  L2–L3 and demonstrably fails at L4 (matching the L4 single-fit coverage
  collapse + the CF rescue in the locked table).
- *L1 growing (expected, NOT a bug):* `SL.glm` is misspecified for both
  nuisances under the Kang–Schafer complex covariates, so eQ/eg have positive
  approximation floors and `sqrt(m)*prod` must eventually grow. This MATCHES
  the locked R02 standard/L1 rows (flat Fully-Aware bias floor
  0.0183/0.0169/0.0172/0.0188; Fully-Aware-CF coverage decaying
  0.949→0.912→0.888→0.848 across `m_total` 60→300) and the R04-vs-smoke level
  gap (0.0625 at `m_total=60` vs 0.0903 at `m_total=300`, a ~13-combined-MCSE
  difference). Report it as an additional Web-D demonstration: the rate
  diagnostic detects not only the interpolating-RF failure (L4) but also the
  misspecified-parametric failure (L1) — exactly the two rungs whose locked
  coverage degrades.
- *L3 borderline (flat-ish but not declining):* still consistent with C1
  "plausibly met" (the test is necessary-condition style); frame alongside
  R04's single-m levels.
- *Genuinely contradictory → STOP and report:* L4 trend flat/declining, or
  L2 trend clearly growing. Either would contradict the locked table (L4
  single-fit coverage collapse; L2 near-nominal coverage) and indicates a bug
  (truth target / OOF refit), not a finding.

## Caveats / assumptions

- **Diagnostic logic.** "Flat/declining √m-scaled product" is *consistent
  with* (not proof of) the o_P(m^{-1/2}) rate — a finite sweep can only show
  the scaled product is non-increasing over the observed range; growth
  *refutes* it (expected at L4 via the interpolating RF, and plausibly at L1
  via the parametric misspecification floor). This phrasing goes in Web-D.
- **m for the rate** = total sampled PSUs (the independent design unit under
  clustering), exactly R04's convention; `prod_int_sqrtn` (n = unit sample
  size) rides along for reference. n grows proportionally to m_total here
  (`base_n0` fixed), so the sqrt(n) trend is the same shape.
- **Truth targets** are R04's primaries: u-integrated `Q0(a,c)` and
  uA-integrated marginal `g0 = E[A|C]` (what C-only learners are consistent
  for); realized-u secondaries (`*_real`) are carried for completeness.
- **Duplication provenance:** the post-draw measurement block in
  `rate_sweep_helpers.R::one_rep_rate_m` is copied verbatim from R04's
  `one_rep_rate` (R04 doesn't export it separately and its files are locked);
  a comment marks the provenance. If R04's helpers ever change, re-sync.
- **No estimator arms are run** (no tmle targeting/coverage here) — this run
  measures nuisance error rates only; coverage-vs-m is R02's table. The two
  join on `(rung, base_m)`.
- **SMOKE uses `truth_M=2e5`** (cheap truth integral). This changes only the
  Monte-Carlo truth Psi (not used by the rate columns) and NOT the population
  or samples (truth uses a separately-scoped seed), so smoke per-rep error
  numbers are directly comparable to full-run ones.
