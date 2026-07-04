# =====================================================================
# fpc_helpers.R  —  R10_fpc run-local helpers (NON-CANONICAL)
#
# Purpose: support the FPC / first-stage-fraction sensitivity run. Nothing here
# edits R/*.R; everything reuses the building blocks that sourcing
# estimators.R / dgp.R / diagnostics.R already exposes:
#   .eif_from_tmle(), make_cf_folds(), .sl(), .se_des(), deff_clust(),
#   draw_sample(), make_population().
#
# Two pieces of NEW logic live here:
#   (1) .se_des_fpc()    — a finite-population-correction (fpc=) variant of the
#                          canonical .se_des() CLUSTERED design SE. Same Hajek
#                          svymean of the same EIF, but the stage-1 svydesign is
#                          built WITH fpc = (population PSUs-per-stratum), so the
#                          between-PSU variance is shrunk by the realized
#                          first-stage sampling fraction f = m/M_pop.
#   (2) cf_eif_for_rep() — recomputes the Fully-Aware-CF arm's per-unit EIF and
#                          point estimate for ONE sample, by replaying exactly the
#                          CF recipe in run_estimators() (UNWEIGHTED out-of-fold
#                          nuisance fits, g floored at g_oof_bound, weighted pooled
#                          tmle targeting). This is needed because run_estimators()
#                          returns eif_fa but NOT the CF EIF, and we must apply the
#                          two SE formulas (with / without fpc) to the SAME EIF.
#                          The Fully-Aware (single weighted fit) EIF is taken
#                          straight from run_estimators()'s diagnostics (eif_fa).
# =====================================================================

# ---- (1) FPC clustered design SE ---------------------------------------------
# Mirrors R/estimators.R::.se_des(..., clustered = TRUE) EXACTLY, except the
# stage-1 design carries fpc = number-of-PSUs-in-the-stratum-population (one value
# per row, constant within a stratum). With a single-stage ids=~cluster design,
# survey applies the stage-1 FPC (1 - m_h/M_h) to the between-PSU variance.
#
# `fpc1` must be a numeric vector the same length as `eif`, equal to the
# POPULATION PSU count of each unit's stratum (J_per_stratum for the sim).
.se_des_fpc <- function(eif, strata, cluster, weight, fpc1) {
  stopifnot(length(fpc1) == length(eif))
  dat <- data.frame(eif = eif, strata = strata, cluster = cluster,
                    weight = weight, fpc1 = fpc1)
  des <- survey::svydesign(ids = ~cluster, strata = ~strata,
                           weights = ~weight, fpc = ~fpc1, data = dat)
  list(se = as.numeric(survey::SE(survey::svymean(~eif, des))),
       df = survey::degf(des))
}

# ---- (2) replay the Fully-Aware-CF arm to obtain its per-unit EIF -------------
# Byte-faithful to run_estimators()'s CF block: same fold construction, same
# UNWEIGHTED .sl() out-of-fold fits, same g_oof_bound flooring, same Q clamp,
# same weighted POOLED tmle() targeting, same .eif_from_tmle() reconstruction.
# Returns list(psi, eif, fold) so the caller can apply .se_des / .se_des_fpc.
cf_eif_for_rep <- function(obs, learners, V_cf = 5L, g_oof_bound = 0.05,
                           W_cols = c("L1", "L2", "L3", "L4")) {
  d   <- obs
  n   <- nrow(d); w <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)

  fold <- make_cf_folds(d$strata, d$cluster, V_cf)
  Q0o <- Q1o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = NULL, learners)   # UNWEIGHTED (ignorability)
    gf <- .sl(d$A[tr], W[tr, ],  weights = NULL, learners)
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g1o <- pmin(pmax(g1o, g_oof_bound), 1 - g_oof_bound)
  Qoo <- pmin(pmax(cbind(Q0o, Q1o), 1e-3), 1 - 1e-3)
  cf <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
             Q = Qoo, g1W = g1o)                                # weighted POOLED targeting
  e_cf <- .eif_from_tmle(cf, d$Y, d$A, w)
  list(psi = e_cf$psi, eif = e_cf$eif, fold = fold)
}

# ---- small summary helper: bias / coverage / se_ratio / mcse (run_sim convention)
# crit uses a t reference with the supplied (mean) design df, matching
# aggregate_sim.R::z_or_t().  b, se are per-rep vectors; df may be a vector.
summarize_arm <- function(b, se, df, Psi, scenario, rung, method, fpc) {
  nrep <- length(b)
  crit <- qt(0.975, pmax(1, df))
  cov  <- mean(abs(b - Psi) <= crit * se)
  data.frame(
    scenario = scenario, rung = rung, method = method, fpc = fpc,
    n_reps = nrep, Psi = Psi,
    bias = mean(b) - Psi, emp_sd = sd(b), mean_se = mean(se),
    se_ratio = mean(se) / sd(b),
    coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrep),
    stringsAsFactors = FALSE
  )
}
