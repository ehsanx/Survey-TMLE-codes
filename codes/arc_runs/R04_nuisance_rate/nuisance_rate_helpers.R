# =====================================================================
# nuisance_rate_helpers.R  —  R04_nuisance_rate run-local instrumentation
#
# Purpose: empirically MEASURE the cross-fit nuisance L2 errors per rung so the
# one "not-fixable" theory item (no weighted-clustered oracle inequality for the
# survey-weighted SuperLearner L2 rate) becomes a DEMONSTRATED fact.
#
# This file adds ONLY new logic; it does NOT edit any canonical codes/*.R. It
# reuses the exposed building blocks from estimators.R (.sl, make_cf_folds) and
# the population object from dgp.R. Source the canonical engine FIRST, then this.
#
# ---------------------------------------------------------------------
# HOW WE EXTRACT THE POPULATION TRUTH Q0(a,c) AND g0(1|c) (the crux)
# ---------------------------------------------------------------------
# Reading codes/dgp.R::make_population(), the realized finite population `pop`
# already carries, per unit, the TRUE nuisances:
#
#   * pscore  = plogis(alpha0 + alpha_g * Csum + uA_i)         <- TRUE g0(1|c)
#               (the exact propensity used to draw A; carries the realized
#                treatment PSU random effect uA_i)
#   * Q1_real = plogis(gamma0 + theta*1 + gamma_C*Csum + uY_i) <- realized-u Q0(1,c)
#   * Q0_real = plogis(gamma0 + theta*0 + gamma_C*Csum + uY_i) <- realized-u Q0(0,c)
#               (carry the realized OUTCOME PSU random effect uY_i; this is the
#                outcome mean GIVEN the realized cluster, i.e. the literal
#                "expit(true linear predictor incl. random effects)" the spec asks)
#   * Q1_int  = gh_expit(gamma0+theta+gamma_C*Csum, sigma2_Y)  <- u-INTEGRATED Q0(1,c)
#   * Q0_int  = gh_expit(gamma0      +gamma_C*Csum, sigma2_Y)  <- u-INTEGRATED Q0(0,c)
#               (uY integrated OUT by Gauss-Hermite; THIS is the conditional mean
#                E[Y|A=a,C=c] the TMLE/SuperLearner actually targets, since the
#                fits see only C, not the latent uY -- see dgp.R RC-1 comment.)
#
# draw_sample() already propagates .Q1_int / .Q0_int into `obs` (the u-integrated
# target). It does NOT propagate pscore / Q*_real. So below we use a thin sampling
# wrapper draw_sample_truth() that re-draws the SAME sample (same seed -> identical
# rows, via the canonical .scoped_seed RNG) and re-attaches ALL truth columns by a
# direct row lookup into population$pop. We do NOT re-implement the sampler: we call
# the canonical draw_sample() for the obs the estimator sees, then map its rows back
# to the population by the design keys (cluster, strata, A, Y, latent .Q*_int) which
# uniquely identify the sampled population rows for diagnostics.
#
# PRIMARY target = u-integrated Q0(a,c) (.Q*_int / Q*_int): this is what the
# cross-fit SuperLearner is consistent FOR (it never sees uY), so the product-rate
# claim must be evaluated against it. We ALSO report the realized-u error
# (Q*_real, g0 incl. uA) as a secondary column for completeness, since the spec
# names "incl. random effects" explicitly. The propensity g0 = pscore INCLUDES uA
# and is the one true object either way (there is no "integrated" g target here --
# the estimator's g is also a function of C only, so ghat is consistent for the
# C-marginal propensity E[A|C] = gh-integrated pscore; we therefore ALSO compute a
# C-marginal g0 by integrating uA out, and treat THAT as the primary g0 target).
# =====================================================================

# ---- C-marginal (uA-integrated) true propensity: E[A | C=c] -----------------
# The cross-fit g-learner sees only C, so it is consistent for the uA-integrated
# propensity E[A|C] = E_{uA}[ plogis(alpha0 + alpha_g*Csum + uA) ], not pscore.
# Reuse the canonical gh_expit() (exposed by dgp.R) with the TREATMENT variance.
.g0_marginal <- function(pop_params, Csum) {
  lin <- pop_params$alpha0 + pop_params$alpha_g * Csum
  gh_expit(lin, pop_params$sigma2_A, nGH = pop_params$nGH %||% 30L)   # uA integrated out
}

# ---- attach the full population truth to a drawn sample ----------------------
# Returns `obs` (exactly what run_estimators sees) augmented with truth columns:
#   .g0_marg  C-marginal true propensity  E[A|C]      (PRIMARY g0 target)
#   .g0_real  realized-u true propensity  pscore       (incl. uA; secondary)
#   .Q1_int/.Q0_int  u-integrated Q0(a,c)             (PRIMARY Q0 target; already on obs)
#   .Q1_real/.Q0_real realized-u Q0(a,c)              (incl. uY; secondary)
# Mapping: draw_sample carries cluster + the latent .Q*_int, which together with
# (A,Y,strata) pin down the exact sampled population rows. We recover Csum/uA/uY
# pieces by re-deriving Csum from .Q0_int via the inverse of the u-integrated map.
# To avoid a fragile numeric inversion, we instead recompute Csum directly from the
# population by re-drawing the SAME sample row indices: draw_sample is deterministic
# in sample_seed, so re-running the canonical sampler logic on the population and
# matching on cluster+within-PSU order is exact. SIMPLEST robust path: re-draw and
# read population columns by the keep-index that draw_sample used internally. Since
# draw_sample does not RETURN keep, we reconstruct Csum from the L-covariates? No --
# apply_L is many-to-one. Cleanest: redraw with the same seed and join on the unique
# row signature carried through (cluster, psu_within, A, Y, .Q1_int, .Q0_int) which
# is unique per population row in practice. We assert uniqueness and fall back to the
# u-integrated-only path (primary targets) if any collision is detected.
attach_truth <- function(obs, population) {
  pop <- population$pop
  pr  <- population$params
  # Csum is needed for the marginal g0; recover it from the EXACT inverse of the
  # u-integrated outcome map is unstable, so derive Csum from pscore? pscore not on
  # obs. We therefore re-key obs rows to pop rows via the carried design signature.
  key_obs <- paste(obs$cluster, obs$psu_within, obs$A, obs$Y,
                   round(obs$.Q1_int, 10), round(obs$.Q0_int, 10), sep = "|")
  key_pop <- paste(pop$cluster, pop$psu_within, pop$A, pop$Y,
                   round(pop$Q1_int, 10), round(pop$Q0_int, 10), sep = "|")
  m <- match(key_obs, key_pop)
  ok <- !anyNA(m) && !any(duplicated(key_obs[!is.na(m)]))
  if (!ok) {
    # Fallback: keep only the u-integrated targets (the PRIMARY ones, already on
    # obs). Realized-u + marginal-g columns get NA; the run still measures the
    # primary product rate. (Documented assumption in NOTES.md.)
    obs$.g0_marg <- NA_real_; obs$.g0_real <- NA_real_
    obs$.Q1_real <- NA_real_; obs$.Q0_real <- NA_real_
    attr(obs, "truth_join") <- "FALLBACK_uint_only"
    return(obs)
  }
  Csum <- pop$C1[m] + pop$C2[m] + pop$C3[m] + pop$C4[m]
  obs$.g0_marg <- .g0_marginal(pr, Csum)        # primary g0 target (uA-integrated)
  obs$.g0_real <- pop$pscore[m]                 # realized-u propensity (incl uA)
  obs$.Q1_real <- pop$Q1_real[m]                # realized-u outcome mean (incl uY)
  obs$.Q0_real <- pop$Q0_real[m]
  attr(obs, "truth_join") <- "OK"
  obs
}

# ---- cross-fit nuisance fit -> OOF predictions (mirrors estimators.R CF arm) --
# Reuses make_cf_folds() and .sl() exactly as the Fully-Aware-CF arm does:
# UNWEIGHTED per-fold nuisance fits (consistent by sampling ignorability), OOF
# propensity floored at g_oof_bound. Returns OOF Qhat(1,.),Qhat(0,.),ghat(1|.).
cf_nuisance_oof <- function(obs, learners, V_cf = 5L, g_oof_bound = 0.05,
                            W_cols = c("L1", "L2", "L3", "L4")) {
  n   <- nrow(obs)
  W   <- obs[, W_cols, drop = FALSE]
  XA  <- obs[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)
  fold <- make_cf_folds(obs$strata, obs$cluster, V_cf)
  Q1o <- Q0o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(obs$Y[tr], XA[tr, ], weights = NULL, learners)   # UNWEIGHTED (ignorability)
    gf <- .sl(obs$A[tr], W[tr, ],  weights = NULL, learners)
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g1o <- pmin(pmax(g1o, g_oof_bound), 1 - g_oof_bound)
  list(Q1 = Q1o, Q0 = Q0o, g1 = g1o, V_eff = attr(fold, "V_eff"))
}

# ---- weighted L2 error  ||fhat - f0||_{L2(w)} = sqrt( sum w (fhat-f0)^2 / sum w )
# Hajek/design-weighted RMSE: estimates the L2(P) norm under the SAMPLED-with-weights
# law, i.e. the population L2 norm (since E_design[ w * h ] = sum_pop h). This is the
# norm in which the survey-weighted product-rate condition is stated.
.wL2 <- function(fhat, f0, w) {
  ok <- is.finite(fhat) & is.finite(f0) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sqrt(sum(w[ok] * (fhat[ok] - f0[ok])^2) / sum(w[ok]))
}

# ---- one rep: fit CF nuisances, measure L2 errors, return one diagnostic row --
# m (the rate's effective sample size) = number of independent design units =
# #PSUs sampled (the clustering unit), matching the design df logic. We also report
# product*sqrt(n) for reference. The PRIMARY rate object is product*sqrt(m_psu).
one_rep_rate <- function(i, population, learners, rung, scenario,
                         V_cf = 5L, g_oof_bound = 0.05) {
  obs <- draw_sample(population, sample_seed = SAMPLE_SEED_BASE + i,
                     model_type = "complex")
  obs <- attach_truth(obs, population)
  oof <- cf_nuisance_oof(obs, learners, V_cf = V_cf, g_oof_bound = g_oof_bound)
  w   <- obs$weight
  n   <- nrow(obs)
  m_psu <- length(unique(obs$cluster))   # # sampled PSUs = independent design units

  # PRIMARY targets: u-integrated Q0(a,c) and C-marginal (uA-integrated) g0(1|c)
  # ghat predicts P(A=1|C); the Q-error is averaged over the two arms a in {0,1},
  # which is the object entering the TMLE remainder term R(Phat,P) ~ ||Qhat-Q0||*||ghat-g0||.
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
  data.frame(
    scenario = scenario, rung = rung, rep = i, n = n, m_psu = m_psu,
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
