# =====================================================================
# estimators_nhanes_iso.R  —  R05 helper: harmonized-floor + weighted-OOF-CF
#
# PURPOSE (de-confound the real-data ladder)
# The locked NHANES ladder (nhanes/nhanes_arc.R -> codes/run_estimators) imports
# the SAME Fully-Aware-vs-Cross-Fit confound as the simulation:
#   * Fully-Aware (FA)  : SINGLE weighted tmle() fit; EIF re-bounded at gbound=1e-3
#   * Fully-Aware-CF    : CROSS-FIT, but per-fold nuisances are UNWEIGHTED
#                         (.sl(weights=NULL)) AND the OOF propensity is floored at
#                         g_oof_bound=0.05.
# So FA-vs-CF differs on THREE axes at once: (single-fit vs cross-fit),
# (weighted vs unweighted nuisance), (1e-3 vs 0.05 floor). This helper isolates
# them on the E3/E4 ladders by adding two NEW arms that change ONE axis at a time,
# WITHOUT editing any canonical R/*.R file:
#
#   A. "Fully-Aware-h05"      = the canonical FA targeting, but the reconstructed
#                               EIF is re-bounded at 0.05 (HARMONIZED floor).
#                               -> same axis as CF's floor; isolates floor effect.
#   B. "Fully-Aware-CF-wOOF"  = cross-fit like FA-CF, but per-fold nuisances are
#                               fit WEIGHTED (.sl(weights = w_train)); OOF
#                               propensity floored at 0.05 (same as FA-CF).
#                               -> isolates de-weighting from cross-fitting.
#
# Together with the canonical FA (1e-3) and FA-CF (unweighted-OOF, 0.05) that we
# also re-emit unchanged, this lets the manuscript attribute the divergence to
# CROSS-FITTING, not to the floor or the de-weighting.
#
# REUSE: this is a thin reimplementation of run_estimators()'s FA + CF arms that
# REUSES the exposed building blocks from R/estimators.R:
#   .eif_from_tmle(fit,Y,A,w,gbound)  — EIF reconstruction with a TUNABLE floor
#   make_cf_folds(strata,cluster,V)   — PSU-within-stratum cross-fit folds
#   .sl(Y,X,weights,learners)         — per-fold SuperLearner predictor closure
#   .se_des(eif,strata,cluster,weight,clustered,nest,inpop) — design SE (Eq 8),
#                                       incl. the FULL-design-subset domain path.
# The inpop / nest / W_cols sub-population + SDMVPSU handling is IDENTICAL to
# run_estimators() (same `sub <- which(inpop)`, same se_des closure), so design df
# and the domain variance match the locked pipeline exactly.
#
# NOTE: the FA POINT estimate psi = weighted.mean(Q1*-Q0*, w) is computed from the
# tmle targeted Qstar and is INDEPENDENT of the EIF floor. Hence Fully-Aware and
# Fully-Aware-h05 share an identical b and differ ONLY in se/df-driven CI width.
# That is the spec's "floor-robust point-estimate divergence": re-flooring cannot
# move the FA point estimate, so the FA-vs-CF point divergence is floor-invariant.
# =====================================================================

# Sourced AFTER R/estimators.R (which defines the building blocks + loads
# SuperLearner/survey/tmle). We do not re-library here to avoid masking.

# ---- iso runner: returns the SAME shape as run_estimators() -----------------
# Args mirror run_estimators(); ADDED: fa_gbound (the harmonized FA floor).
# Default fa_gbound = 0.05 to match the CF OOF floor (g_oof_bound).
run_estimators_iso <- function(obs, learners = "SL.glm", V_cf = 5L,
                               g_oof_bound = 0.05, fa_gbound = 0.05,
                               W_cols = c("L1", "L2", "L3", "L4"),
                               nest = FALSE, inpop = NULL) {
  # ----- domain slice: identical contract to run_estimators() -----
  sub <- if (is.null(inpop)) seq_len(nrow(obs)) else which(inpop)
  d   <- obs[sub, , drop = FALSE]
  n   <- nrow(d); w <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)
  # design-SE closure: simulation path uses domain rows; domain path builds the
  # FULL design and subset()s it to inpop (PSUs absent from the domain still count
  # toward the df). Byte-identical to run_estimators()'s se_des.
  se_des <- function(eif, clustered)
    if (is.null(inpop)) .se_des(eif, d$strata, d$cluster, w, clustered)
    else .se_des(eif, obs$strata, obs$cluster, obs$weight, clustered, nest = nest, inpop = inpop)
  res <- list()

  ## ================= Fully-Aware (canonical, 1e-3 floor) =================
  ## Re-emit the locked FA arm UNCHANGED so the harmonized arm has a same-fit
  ## comparator (identical tmle targeting; only the EIF floor differs).
  fa <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
             Q.SL.library = learners, g.SL.library = learners)
  e_fa_1e3 <- .eif_from_tmle(fa, d$Y, d$A, w, gbound = 1e-3)   # canonical floor
  s_fa_1e3 <- se_des(e_fa_1e3$eif, TRUE)
  res$FullyAware <- data.frame(method = "Fully-Aware",
                               b = e_fa_1e3$psi, se = s_fa_1e3$se, df = s_fa_1e3$df)

  ## ================= Fully-Aware-h05 (HARMONIZED floor) =================
  ## SAME tmle fit `fa`; ONLY re-bound the reconstructed EIF at fa_gbound (0.05).
  ## psi is unchanged (floor-invariant); se/df may change via the H clever cov.
  e_fa_h <- .eif_from_tmle(fa, d$Y, d$A, w, gbound = fa_gbound)
  s_fa_h <- se_des(e_fa_h$eif, TRUE)
  res$FullyAware_h05 <- data.frame(method = "Fully-Aware-h05",
                                   b = e_fa_h$psi, se = s_fa_h$se, df = s_fa_h$df)

  ## ================= Fully-Aware-CF (canonical: UNWEIGHTED OOF) =================
  ## Re-emit the locked CF arm UNCHANGED: per-fold nuisances UNWEIGHTED, OOF
  ## propensity floored at g_oof_bound, weighted POOLED tmle targeting.
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
             Q = Qoo, g1W = g1o)
  e_cf <- .eif_from_tmle(cf, d$Y, d$A, w)                     # canonical 1e-3 EIF floor
  s_cf <- se_des(e_cf$eif, TRUE)
  res$FullyAwareCF <- data.frame(method = "Fully-Aware-CF",
                                 b = e_cf$psi, se = s_cf$se, df = s_cf$df)

  ## ================= Fully-Aware-CF-wOOF (WEIGHTED OOF) =================
  ## NEW ARM: same cross-fit STRUCTURE as FA-CF (same fold object), but per-fold
  ## nuisances are fit WEIGHTED via .sl(weights = w[tr]). OOF propensity floored at
  ## g_oof_bound (same as FA-CF). Isolates de-weighting from cross-fitting: if this
  ## arm tracks FA-CF (not single-fit FA), cross-fitting -- not de-weighting -- is
  ## what stabilizes the ladder.
  Q0w <- Q1w <- g1w <- numeric(n)
  for (v in sort(unique(fold))) {                            # REUSE the same folds
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = w[tr], learners)  # WEIGHTED OOF
    gf <- .sl(d$A[tr], W[tr, ],  weights = w[tr], learners)  # WEIGHTED OOF
    Q1w[ho] <- qf(Xa1[ho, ]); Q0w[ho] <- qf(Xa0[ho, ]); g1w[ho] <- gf(W[ho, ])
  }
  g1w <- pmin(pmax(g1w, g_oof_bound), 1 - g_oof_bound)
  Qww <- pmin(pmax(cbind(Q0w, Q1w), 1e-3), 1 - 1e-3)
  cfw <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
              Q = Qww, g1W = g1w)
  ## DIVERGENCE GUARD. Weighted per-fold nuisance fits re-introduce propensity
  ## separation on rare-exposure survey domains (e.g. E4, ~6% exposed), and the
  ## pooled weighted targeting then DIVERGES (eps -> ~1e13; verified in the R05
  ## smoke). A diverged fluctuation yields a meaningless point estimate (e.g.
  ## -0.50). We therefore report this arm as NON-CONVERGENT rather than emitting a
  ## garbage number; the divergence itself is the empirical justification for the
  ## primary estimator's UNWEIGHTED out-of-fold fits (valid under S _||_ (A,Y)|C).
  eps_cfw_max <- max(abs(cfw$epsilon))
  EPS_DIVERGE <- 20
  if (is.finite(eps_cfw_max) && eps_cfw_max <= EPS_DIVERGE) {
    e_cfw <- .eif_from_tmle(cfw, d$Y, d$A, w)
    s_cfw <- se_des(e_cfw$eif, TRUE)
    res$FullyAwareCFw <- data.frame(method = "Fully-Aware-CF-wOOF",
                                    b = e_cfw$psi, se = s_cfw$se, df = s_cfw$df)
  } else {
    res$FullyAwareCFw <- data.frame(method = "Fully-Aware-CF-wOOF",
                                    b = NA_real_, se = NA_real_, df = s_fa_1e3$df)
  }

  # ---- per-run diagnostics (extends run_estimators()'s drow with iso fields) ----
  g_fa <- fa$g$g1W
  drow <- data.frame(
    eps_fa       = max(abs(fa$epsilon)),
    g_fa_min     = min(g_fa), g_fa_max = max(g_fa),
    eps_cf       = max(abs(cf$epsilon)),
    g_cf_min     = min(g1o),  g_cf_max = max(g1o),         # UNWEIGHTED-OOF propensity
    eps_cfw      = max(abs(cfw$epsilon)),
    g_cfw_min    = min(g1w),  g_cfw_max = max(g1w),        # WEIGHTED-OOF propensity
    cf_V_eff     = attr(fold, "V_eff"),
    fa_gbound    = fa_gbound, g_oof_bound = g_oof_bound,
    # floor-SENSITIVE near-bound mass (reported but NOT the decision statistic):
    g_fa_near_1e3 = mean(g_fa < 1e-3 | g_fa > 1 - 1e-3),
    g_fa_near_h05 = mean(g_fa < fa_gbound | g_fa > 1 - fa_gbound)
  )
  diag <- list(eif_fa = e_fa_1e3$eif,      # canonical FA EIF for deff_clust (matches locked)
               eif_fa_h05 = e_fa_h$eif,    # harmonized-floor FA EIF (for sensitivity)
               g_fa = g_fa, g_cf = g1o, g_cfw = g1w, drow = drow)
  list(results = do.call(rbind, res), diagnostics = diag)
}
