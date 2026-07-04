# =====================================================================
# cf_arms.R  —  R07_informative helper: the two cross-fitted arms compared
#
# Lives in THIS run folder. Reuses the building blocks EXPOSED by sourcing
# R/estimators.R (do NOT edit that file):
#   make_cf_folds(strata, cluster, V)         PSU-level CF fold split
#   .sl(Y, X, weights, learners)              SuperLearner -> predictor closure
#   .eif_from_tmle(fit, Y, A, w, gbound)      reconstruct per-unit EIF
#   .se_des(eif, strata, cluster, weight, clustered)  Eq-8 design SE + df
# and the `tmle` package for the (weighted, pooled) targeting step.
#
# Both arms are IDENTICAL except for ONE line: whether the per-fold out-of-fold
# (OOF) nuisance fits use survey weights.
#   * "Fully-Aware-CF-unwt"  : .sl(... weights = NULL ...)   <- the paper DEFAULT
#       (consistent ONLY under sampling-ignorability  S _||_ (A,Y) | C).
#   * "Fully-Aware-CF-wt"    : .sl(... weights = w[tr] ...)  <- the defended
#       alternative that remains consistent under INFORMATIVE-beyond-C sampling.
# Everything downstream (g flooring, pooled weighted tmle() targeting, design SE)
# is shared, so any bias DIFFERENCE between the two arms is attributable solely
# to de-weighting the OOF fits -- exactly the quantity R07 measures vs rho.
# =====================================================================

# ---- one cross-fitted arm; `oof_weighted` toggles the single contested choice -
# Returns data.frame(method, b, se, df) plus a small diagnostics list, matching
# the canonical run_estimators() row schema (method,b,se,df).
.cf_arm <- function(d, learners, oof_weighted, method_name, fold,
                    g_oof_bound = 0.05,
                    W_cols = c("L1", "L2", "L3", "L4")) {
  n   <- nrow(d); w <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)

  # `fold` is SHARED between the two arms (built once in run_cf_pair) so the only
  # difference between unwt/wt is the OOF weighting -- a clean PAIRED comparison.
  Q0o <- Q1o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    # >>> THE ONE CONTESTED LINE <<<  OOF nuisance fits weighted vs unweighted
    wtr <- if (oof_weighted) w[tr] else NULL
    qf <- .sl(d$Y[tr], XA[tr, ], weights = wtr, learners)
    gf <- .sl(d$A[tr], W[tr, ],  weights = wtr, learners)
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g1o <- pmin(pmax(g1o, g_oof_bound), 1 - g_oof_bound)
  Qoo <- pmin(pmax(cbind(Q0o, Q1o), 1e-3), 1 - 1e-3)
  # pooled WEIGHTED targeting (shared by both arms; this is the canonical CF step)
  cf <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
             Q = Qoo, g1W = g1o)
  e_cf <- .eif_from_tmle(cf, d$Y, d$A, w)
  s_cf <- .se_des(e_cf$eif, d$strata, d$cluster, w, clustered = TRUE)
  list(
    row = data.frame(method = method_name, b = e_cf$psi, se = s_cf$se, df = s_cf$df,
                     stringsAsFactors = FALSE),
    diag = data.frame(method = method_name,
                      eps_cf = max(abs(cf$epsilon)),
                      g_cf_min = min(g1o), g_cf_max = max(g1o),
                      cf_V_eff = attr(fold, "V_eff"),
                      stringsAsFactors = FALSE),
    eif = e_cf$eif
  )
}

# ---- run BOTH CF arms on one sample and return stacked rows + diagnostics ----
# `obs` is a draw_sample_infsamp() output. We do NOT call run_estimators() here
# because we only need the two CF arms for the rho sweep (cheaper, focused).
run_cf_pair <- function(obs, learners, V_cf = 5L, g_oof_bound = 0.05,
                        W_cols = c("L1", "L2", "L3", "L4")) {
  d <- obs
  # build the PSU-level CF folds ONCE and share across both arms (paired design)
  fold <- make_cf_folds(d$strata, d$cluster, V_cf)
  unwt <- .cf_arm(d, learners, oof_weighted = FALSE,
                  method_name = "Fully-Aware-CF-unwt", fold = fold,
                  g_oof_bound = g_oof_bound, W_cols = W_cols)
  wt   <- .cf_arm(d, learners, oof_weighted = TRUE,
                  method_name = "Fully-Aware-CF-wt", fold = fold,
                  g_oof_bound = g_oof_bound, W_cols = W_cols)
  list(
    results = rbind(unwt$row, wt$row),
    diag    = rbind(unwt$diag, wt$diag),
    eif_unwt = unwt$eif    # used for the clustering DEFF diagnostic
  )
}
