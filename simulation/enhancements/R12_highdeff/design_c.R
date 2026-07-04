# =====================================================================
# design_c.R  —  R12_highdeff helper: the "Design C" high-DEFF scenario
#
# PURPOSE -------------------------------------------------------------
# The headline sim (`standard`) realises an EIF clustering design effect
# (deff_clust on the Fully-Aware EIF) of only ~1.25-1.40. Reviewers asked
# for a regime where deff_clust reaches ~2.5-4 to stress-test the
# clustering-aware SE that the method is sold for.
#
# Design C raises the two PSU random-effect variances that drive EIF
# clustering and enlarges the within-PSU (stage-2) sample size:
#   * sigma2_Y : OUTCOME-logit PSU effect   1.5 -> 4.5  (critics: 3-4)
#   * sigma2_A : TREATMENT-logit PSU effect 0.8 -> 2.2  (critics: 1.5-2)
#   * base_n0  : units sampled per PSU       25 -> 55   (critics: 40-50)
# sigma2_C, beta_strat, alpha_g, p_treat_target are left at the canonical
# `standard` defaults so the estimand and overlap regime are unchanged.
#
# TUNED (local GLM-TMLE EIF probe, 20 reps): these land mean deff_clust
# ~2.59 (median 2.53, range ~1.6-4.1) -- inside the requested 2.5-4. The
# headline `standard` defaults give ~1.32 (probe) / ~1.25-1.40 (locked SL
# runs). Note the deep-RF (L4) SL EIF typically realises a SLIGHTLY HIGHER
# DEFF than this GLM proxy, so ~2.5 is a conservative floor; the SMOKE run
# must CONFIRM the realised deff_clust before the full submit.
#
# NON-DESTRUCTIVE: this file is a THIN WRAPPER over the canonical
# make_population()/draw_sample() in R/dgp.R, passing the overrides
# as ORDINARY ARGUMENTS (both functions already expose them). It does
# NOT edit R/dgp.R or define a new scenario inside it. The wrappers
# are sourced ONLY by this run's run.R.
#
# TUNABILITY: every knob is read from an env var (with the Design C
# default) so the smoke-gate retune ("if realized DEFF < ~2.5, raise the
# variances") can be done WITHOUT editing code -- just re-export the var.
# =====================================================================

# null-coalescing helper is defined in dgp.R (`%||%`); dgp.R is sourced
# before this file in run.R, so it is available here.

# ---- Design C knobs (env-overridable for the retune loop) -------------
.envnum <- function(v, default) { x <- Sys.getenv(v); if (nzchar(x)) as.numeric(x) else default }
.enviv  <- function(v, default) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else default }

# Defaults TUNED to land deff_clust in ~2.5-4 (verify in the smoke run).
DESIGNC_SIGMA2_Y <- .envnum("DC_SIGMA2_Y", 4.5)   # was 1.5  -> stronger outcome-PSU clustering
DESIGNC_SIGMA2_A <- .envnum("DC_SIGMA2_A", 2.2)   # was 0.8  -> stronger clever-covariate clustering
DESIGNC_SIGMA2_C <- .envnum("DC_SIGMA2_C", 0.3)   # unchanged from `standard`
DESIGNC_BASE_N0  <- .enviv ("DC_BASE_N0",  55L)   # was 25   -> bigger within-PSU n (more correlated units/PSU)
DESIGNC_BASE_M   <- .enviv ("DC_BASE_M",   6L)    # PSUs/stratum sampled (canonical `standard` default)
DESIGNC_ALPHA_STRAT <- .envnum("DC_ALPHA_STRAT", 2.0)  # weight-informativeness (canonical default)

# ---- Design C population: canonical `standard` geometry, raised variances ----
# pop_seed is threaded from POP_SEED so the population is reproducible and
# identical across workers/tasks (matches run_sim.R's contract).
make_population_designC <- function(pop_seed = 20260606L,
                                    sigma2_Y = DESIGNC_SIGMA2_Y,
                                    sigma2_A = DESIGNC_SIGMA2_A,
                                    sigma2_C = DESIGNC_SIGMA2_C,
                                    truth_M  = 2e6L,
                                    model_type = "complex") {
  make_population(
    scenario   = "standard",        # reuse standard's H/J/M geometry (N=120k, 600 PSUs)
    model_type = model_type,
    pop_seed   = pop_seed,
    sigma2_Y   = sigma2_Y,          # <-- Design C override
    sigma2_A   = sigma2_A,          # <-- Design C override
    sigma2_C   = sigma2_C,
    truth_M    = truth_M
    # beta_strat, alpha_g, p_treat_target, gamma0, gamma_C, te_log_odds:
    # all left at make_population() defaults == the `standard` headline values.
  )
}

# ---- Design C sample: bigger within-PSU n via base_n0 override --------------
# draw_sample() already takes base_m / base_n0 / alpha_strat as args, so we
# just forward the Design C values. (For scenario=="standard" the canonical
# defaults are base_m=6L, base_n0=25L; we raise base_n0 to enlarge stage-2 n.)
draw_sample_designC <- function(population, sample_seed = NULL,
                                model_type = "complex",
                                base_m  = DESIGNC_BASE_M,
                                base_n0 = DESIGNC_BASE_N0,
                                alpha_strat = DESIGNC_ALPHA_STRAT) {
  draw_sample(population, sample_seed = sample_seed, model_type = model_type,
              base_m = base_m, base_n0 = base_n0, alpha_strat = alpha_strat,
              oversample = TRUE)
}

# Convenience: a one-line record of the realized Design C knobs (for manifest/logs).
designC_knobs <- function() list(
  scenario_base = "standard",
  sigma2_Y = DESIGNC_SIGMA2_Y, sigma2_A = DESIGNC_SIGMA2_A, sigma2_C = DESIGNC_SIGMA2_C,
  base_n0 = DESIGNC_BASE_N0, base_m = DESIGNC_BASE_M, alpha_strat = DESIGNC_ALPHA_STRAT,
  note = "high-DEFF Design C: raised PSU RE variances + within-PSU n to push deff_clust ~2.5-4"
)
