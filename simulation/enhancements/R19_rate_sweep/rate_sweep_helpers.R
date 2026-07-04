# =====================================================================
# rate_sweep_helpers.R  —  R19_rate_sweep run-local wrapper (NON-DESTRUCTIVE)
#
# Purpose: cross the R04 nuisance-error instrumentation with the R02 base_m
# sweep, producing the sqrt(m)-scaled product-error TREND per rung — the
# re-review's sharpest evidence ask (a single-m diagnostic cannot test a RATE;
# see Writing/comments/phase4-arc-sim-specs.md, "√m-scaled product-error
# diagnostic across the F.5 m-sweep").
#
# This file adds ONLY new logic; it does NOT edit any canonical codes/*.R nor
# any other arc_run's files. It REQUIRES that run.R has already sourced, in
# order:
#   1. the canonical engine (config/dgp/estimators/diagnostics/learners)
#   2. R04's OWN helpers file
#      codes/arc_runs/R04_nuisance_rate/nuisance_rate_helpers.R
#      (sourced from R04's folder, NOT copied) which provides:
#        attach_truth()     join sampled rows back to population truth columns
#        .g0_marginal()     uA-integrated true propensity E[A|C]
#        cf_nuisance_oof()  CF out-of-fold nuisance fits (reuses .sl/make_cf_folds)
#        .wL2()             design-weighted L2 norm
#        one_rep_rate()     the R04 per-rep driver (FIXED at the scenario-default
#                           base_m — it does NOT accept base_m, hence this wrapper)
#
# WHY A WRAPPER: R04's one_rep_rate() hardcodes
#   draw_sample(population, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
# i.e. the scenario-default base_m (6 for 'standard'). R02's one_rep_m() shows
# the m knob: draw_sample()'s EXISTING `base_m` argument (stage-1 sampled PSUs
# per stratum; m_total = H * base_m). one_rep_rate_m() below replicates ONLY the
# draw step (forwarding base_m) and then reuses R04's internals for everything
# else. The post-draw measurement arithmetic is duplicated VERBATIM from R04
# (provenance-marked) because R04 does not expose it as a standalone function
# and we must not edit R04's files.
#
# m CONVENTION (matches R04): m = m_psu = TOTAL number of sampled PSUs
# (= H * base_m exactly, since draw_sample keeps m_h balanced across strata:
# standard H=10 -> m_total in {60, 120, 200, 300} for base_m in {6, 12, 20, 30}).
# prod_int_sqrtm = prod_int * sqrt(m_psu). base_m AND m_total are recorded per
# row so the numbers join the existing R02/R04 tables directly.
# =====================================================================

# ---- one rep at a chosen sampled-PSU count base_m ----------------------------
# i = GLOBAL rep index. sample_seed = SAMPLE_SEED_BASE + i, so for a given
# base_m the drawn samples are byte-identical to R02's m-sweep draws — and at
# base_m = 6 (the 'standard' default) byte-identical to the locked headline run
# AND to R04's draws. This preserves the cross-run comparability property.
one_rep_rate_m <- function(i, population, learners, rung, scenario, base_m,
                           V_cf = 5L, g_oof_bound = 0.05) {
  # --- draw step: the ONLY logic that differs from R04::one_rep_rate ---------
  # (mirrors R02_largem_sweep/helpers.R::one_rep_m's forwarding of base_m)
  obs <- draw_sample(population, sample_seed = SAMPLE_SEED_BASE + i,
                     model_type = "complex", base_m = base_m)

  # --- measurement block: copied VERBATIM from codes/arc_runs/
  # R04_nuisance_rate/nuisance_rate_helpers.R::one_rep_rate (post-draw body).
  # It CALLS R04's exported pieces (attach_truth, cf_nuisance_oof, .wL2);
  # only the connecting arithmetic is re-stated here. Keep in sync with R04
  # if that file ever changes (it is locked; this is a provenance marker).
  obs <- attach_truth(obs, population)
  oof <- cf_nuisance_oof(obs, learners, V_cf = V_cf, g_oof_bound = g_oof_bound)
  w   <- obs$weight
  n   <- nrow(obs)
  m_psu <- length(unique(obs$cluster))   # total sampled PSUs = H * base_m here

  # PRIMARY targets: u-integrated Q0(a,c) and C-marginal (uA-integrated) g0(1|c)
  eQ1_int <- .wL2(oof$Q1, obs$.Q1_int, w)
  eQ0_int <- .wL2(oof$Q0, obs$.Q0_int, w)
  eQ_int  <- sqrt((eQ1_int^2 + eQ0_int^2) / 2)          # combined Q L2 error (primary)
  eg_marg <- .wL2(oof$g1, obs$.g0_marg, w)              # g L2 error (primary)

  # SECONDARY (realized-u) targets: outcome mean incl uY, propensity incl uA
  eQ1_real <- .wL2(oof$Q1, obs$.Q1_real, w)
  eQ0_real <- .wL2(oof$Q0, obs$.Q0_real, w)
  eQ_real  <- sqrt((eQ1_real^2 + eQ0_real^2) / 2)
  eg_real  <- .wL2(oof$g1, obs$.g0_real, w)

  prod_int  <- eQ_int  * eg_marg                         # ||Qhat-Q0|| * ||ghat-g0|| (primary)
  prod_real <- eQ_real * eg_real
  # --- end copied block -------------------------------------------------------

  # R19 row = R04 row + the swept-knob columns LEADING (base_m, m_total), so
  # the output joins the existing R02/R04 tables on (scenario, rung, base_m).
  data.frame(
    scenario = scenario, rung = rung,
    base_m = base_m, m_total = population$params$H * base_m,
    rep = i, n = n, m_psu = m_psu,
    # primary (u-integrated Q, uA-integrated g) -- the rate object
    eQ_int = eQ_int, eg_int = eg_marg, prod_int = prod_int,
    prod_int_sqrtm = prod_int * sqrt(m_psu),
    prod_int_sqrtn = prod_int * sqrt(n),
    eQ1_int = eQ1_int, eQ0_int = eQ0_int,
    # secondary (realized-u Q incl uY, realized-u g incl uA)
    eQ_real = eQ_real, eg_real = eg_real, prod_real = prod_real,
    prod_real_sqrtm = prod_real * sqrt(m_psu),
    V_eff = oof$V_eff, truth_join = attr(obs, "truth_join"),
    stringsAsFactors = FALSE
  )
}
