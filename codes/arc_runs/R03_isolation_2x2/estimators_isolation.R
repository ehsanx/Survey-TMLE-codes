# =====================================================================
# estimators_isolation.R  —  R03 helper: the 2x2 {cross-fit} x {weighting}
# isolation, all at a HARMONIZED propensity floor g_floor.
#
# PURPOSE (see NOTES.md): the locked engine's Fully-Aware-vs-CF contrast
# CONFOUNDS two factors -- the FA arm is single-fit + weighted-nuisance +
# EIF floor 1e-3, while the CF arm is cross-fit + UNWEIGHTED-OOF-nuisance +
# OOF/EIF floor 0.05. We cannot tell whether the coverage rescue comes from
# CROSS-FITTING or from DE-WEIGHTING the nuisances. This helper crosses the
# two factors and harmonizes the floor so only one thing varies per cell:
#
#        |  weighted-OOF nuisance   |  unweighted-OOF nuisance
#  ------+--------------------------+--------------------------
#  single|  SF-W   (FA-like)        |  SF-U   (the DANGEROUS tail)
#  cross |  CF-W                    |  CF-U   (engine FA-CF-like)
#
# All four arms:
#   * estimate Q (outcome) and g (propensity) themselves, then hand
#     pre-computed (Q, g1W) to tmle() for the SAME weighted POOLED targeting
#     the engine's CV/CF arms use -- so the ONLY differences across cells are
#     (a) in-sample vs out-of-fold nuisance fits and (b) obsWeights to .sl().
#   * floor g at the SAME g_floor in BOTH the OOF/in-sample predictions AND
#     the EIF reconstruction (.eif_from_tmle(..., gbound = g_floor)) -- this
#     removes the 1e-3-vs-0.05 asymmetry that the engine bakes in.
#   * compute the clustered design SE (Eq 8) via the engine's .se_des().
#
# REUSES the building blocks exposed by sourcing estimators.R:
#   .sl(Y, X, weights, learners)            -> predictor closure
#   .eif_from_tmle(fit, Y, A, w, gbound)    -> reconstructed per-unit EIF
#   .se_des(eif, strata, cluster, weight, clustered)
#   make_cf_folds(strata, cluster, V)       -> whole-PSU-within-stratum folds
# It does NOT edit codes/*.R.
#
# Pre-stated expectation (the demonstration this is meant to confirm):
#   single-fit  (EITHER weighting) UNDER-covers at L4_aggressive;
#   cross-fit   (EITHER weighting) RECOVERS ~nominal coverage,
#   => CROSS-FITTING is the active ingredient, not de-weighting.
# The decision rule lives in run.R / NOTES.md: if SF-U already covers well at
# L4, de-weighting (not cross-fitting) is doing the work -> STOP and report.
# =====================================================================

# Q clamp used by the engine before handing Q to tmle() (estimators.R lines
# 162, 178). Keep identical so targeting behaves the same as the locked arms.
.Q_CLAMP <- 1e-3

# ---- fit nuisances either single-fit (in-sample) or cross-fit (out-of-fold) --
# Returns a list with the n x 2 outcome matrix Q = [Q(A=0,W), Q(A=1,W)] and the
# propensity vector g1 (already floored at g_floor), plus the realized CF V.
# weighted  : TRUE -> pass obsWeights = w to .sl (weighted-OOF / weighted-SF)
#             FALSE -> .sl(weights = NULL) (UNWEIGHTED, sampling-ignorability)
# crossfit  : TRUE -> PSU-level out-of-fold fits; FALSE -> full-sample in-sample
.fit_nuisances <- function(d, learners, weighted, crossfit, g_floor, V_cf = 5L) {
  n   <- nrow(d)
  w   <- d$weight
  W   <- d[, c("L1", "L2", "L3", "L4"), drop = FALSE]
  XA  <- d[, c("A", "L1", "L2", "L3", "L4"), drop = FALSE]
  Xa1 <- data.frame(A = 1, W)
  Xa0 <- data.frame(A = 0, W)
  wt  <- if (weighted) w else NULL          # the ONLY weighting lever

  if (!crossfit) {
    ## ---- SINGLE-FIT: one Q-fit + one g-fit on the FULL sample, predict in-sample
    qf <- .sl(d$Y, XA, weights = wt, learners)
    gf <- .sl(d$A, W,  weights = wt, learners)
    Q1 <- qf(Xa1); Q0 <- qf(Xa0); g1 <- gf(W)
    V_eff <- NA_integer_
  } else {
    ## ---- CROSS-FIT: whole-PSU-within-stratum folds; fit on train, predict OOF
    fold <- make_cf_folds(d$strata, d$cluster, V_cf)
    Q1 <- Q0 <- g1 <- numeric(n)
    for (v in sort(unique(fold))) {
      tr <- which(fold != v); ho <- which(fold == v)
      wtr <- if (weighted) w[tr] else NULL
      qf  <- .sl(d$Y[tr], XA[tr, ], weights = wtr, learners)
      gf  <- .sl(d$A[tr], W[tr, ],  weights = wtr, learners)
      Q1[ho] <- qf(Xa1[ho, ]); Q0[ho] <- qf(Xa0[ho, ]); g1[ho] <- gf(W[ho, ])
    }
    V_eff <- attr(fold, "V_eff")
  }

  # HARMONIZED floor on the propensity (same g_floor for every arm)
  g1 <- pmin(pmax(g1, g_floor), 1 - g_floor)
  Q  <- pmin(pmax(cbind(Q0, Q1), .Q_CLAMP), 1 - .Q_CLAMP)
  list(Q = Q, g1 = g1, V_eff = V_eff)
}

# ---- one arm: pre-computed (Q, g1) -> weighted pooled tmle targeting -> SE ----
.one_arm <- function(d, nu, method_label, g_floor) {
  w <- d$weight
  W <- d[, c("L1", "L2", "L3", "L4"), drop = FALSE]
  fit <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
              Q = nu$Q, g1W = nu$g1)               # weighted POOLED targeting
  # HARMONIZED floor in the EIF reconstruction too (gbound = g_floor, not 1e-3)
  e <- .eif_from_tmle(fit, d$Y, d$A, w, gbound = g_floor)
  s <- .se_des(e$eif, d$strata, d$cluster, w, clustered = TRUE)
  # DIVERGENCE GUARD: weighted per-fold nuisance fits can separate and the pooled
  # targeting then diverges (eps -> ~1e13; observed in the R05 NHANES smoke). A
  # diverged fluctuation yields a meaningless point estimate, so mark such reps
  # NON-CONVERGENT (NA) -- a single 1e13 rep would otherwise wreck the 1000-rep
  # CF-W/SF-W mean. The aggregator counts these (n_diverged) and uses na.rm.
  eps_max  <- max(abs(fit$epsilon))
  diverged <- !is.finite(eps_max) || eps_max > 20
  list(
    row  = data.frame(method = method_label,
                      b  = if (diverged) NA_real_ else e$psi,
                      se = if (diverged) NA_real_ else s$se,
                      df = s$df, diverged = diverged,
                      stringsAsFactors = FALSE),
    eif  = e$eif,
    drow = data.frame(method = method_label,
                      eps = eps_max,
                      g_min = min(nu$g1), g_max = max(nu$g1),
                      V_eff = nu$V_eff, stringsAsFactors = FALSE)
  )
}

# =====================================================================
# run_isolation(obs, learners, g_floor = 0.05, V_cf = 5L)
#   Returns the SAME shape as run_estimators(): list(results, diagnostics).
#   results: data.frame(method, b, se, df) with the four arm labels
#     SF-W  single-fit / weighted-nuisance
#     SF-U  single-fit / unweighted-nuisance   (the dangerous tail)
#     CF-W  cross-fit  / weighted-OOF-nuisance
#     CF-U  cross-fit  / unweighted-OOF-nuisance
#   diagnostics: list(eif_sfw, drow) -- eif of the SF-W arm is returned for the
#     DEFF/ICC audit (deff_clust), matching how run_sim uses eif_fa.
# =====================================================================
run_isolation <- function(obs, learners = "SL.glm", g_floor = 0.05, V_cf = 5L) {
  d <- obs                                  # whole-sample (no domain in the sim)

  cells <- list(
    SF_W = list(label = "SF-W", weighted = TRUE,  crossfit = FALSE),
    SF_U = list(label = "SF-U", weighted = FALSE, crossfit = FALSE),
    CF_W = list(label = "CF-W", weighted = TRUE,  crossfit = TRUE),
    CF_U = list(label = "CF-U", weighted = FALSE, crossfit = TRUE)
  )

  arms <- lapply(cells, function(c) {
    nu <- .fit_nuisances(d, learners, weighted = c$weighted,
                         crossfit = c$crossfit, g_floor = g_floor, V_cf = V_cf)
    .one_arm(d, nu, c$label, g_floor)
  })

  results <- do.call(rbind, lapply(arms, `[[`, "row"))
  drow    <- do.call(rbind, lapply(arms, `[[`, "drow"))
  rownames(results) <- NULL; rownames(drow) <- NULL

  list(results = results,
       diagnostics = list(eif_sfw = arms$SF_W$eif, drow = drow))
}
