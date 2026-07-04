# =====================================================================
# estimators_harmonized.R  —  R16 helper: the FIVE paper arms under a single
# HARMONIZED truncation rule (spec item A9a, harmonized-truncation headline
# simulation; Writing/comments/phase4-arc-sim-specs.md).
#
# PURPOSE. The locked headline grid (codes/estimators.R::run_estimators) lets
# the truncation rule vary WITH the arm: the single-fit arms reconstruct the
# EIF at gbound=1e-3 from tmle()'s internally-bounded g, while the CV/CF arms
# floor their OOF/internal-CV propensity at g_oof_bound=0.05 before targeting.
# So the Figure 1 / Table 1 single-fit-vs-CF contrast confounds CROSS-FITTING
# with the FLOOR. This helper re-runs ALL FIVE paper arms with the SAME rule:
#
#   * g predictions floored to [g_floor, 1-g_floor] BEFORE targeting,
#   * Q predictions truncated to [q_lo, 1-q_lo] BEFORE targeting,
#   * targeting always via tmle(Y, A, W, family="binomial",
#       obsWeights=<arm>, Q=<pre-fit>, g1W=<pre-fit>)   (pre-fit nuisances; the
#       pattern validated by R03_isolation_2x2 / estimators_isolation.R),
#   * EIF via .eif_from_tmle(fit, Y, A, w_arm, gbound = g_floor),
#   * SEs via the engine's .se_des().
#
# Arms (results$method labels carry the "-h" suffix so the aggregate can join
# them to the locked arms without name collisions):
#   1. Fully-Aware-h     weighted single-fit nuisances (.sl(weights=w)) for Q on
#                        (A,W) and g on W, in-sample predictions, WEIGHTED
#                        targeting, clustered design SE (Eq 8).
#   2. Partially-Aware-h SAME fit/EIF as arm 1; .se_des(clustered=FALSE).
#   3. Non-Aware-h       UNWEIGHTED nuisances + UNWEIGHTED targeting
#                        (obsWeights=rep(1,n)); EIF with w=1;
#                        SE = sqrt(var(eif)/n), df = n-1.
#   4. Fully-Aware-CV-h  (only if length(learners) > 1) the engine's CV block
#                        (estimators.R ~160-182: surveyCV::folds.svy validRows,
#                        WEIGHTED internal-CV nuisances = the sim-path foil) at
#                        the harmonized g_floor/q_lo; weighted targeting;
#                        clustered SE.
#   5. Fully-Aware-CF-h  the engine's CF logic (make_cf_folds, UNWEIGHTED
#                        per-fold .sl fits, OOF predictions) with OOF g floored
#                        at g_floor; weighted POOLED targeting; clustered SE.
#                        Identical to the locked CF arm except the EIF gbound
#                        1e-3 -> g_floor (a no-op at the default g_floor=0.05,
#                        since the OOF g handed to tmle is already in
#                        [0.05, 0.95] -- stated for harmonization completeness).
#
# REUSES the building blocks exposed by sourcing codes/estimators.R:
#   .sl(Y, X, weights, learners)             -> SuperLearner predictor closure
#   .eif_from_tmle(fit, Y, A, w, gbound)     -> reconstructed per-unit EIF
#   .se_des(eif, strata, cluster, weight, clustered)
#   make_cf_folds(strata, cluster, V)        -> whole-PSU-within-stratum folds
# It does NOT edit any codes/*.R.
#
# Pre-stated expectation (per R03_isolation_2x2): the L4 single-fit
# undercoverage SURVIVES harmonization (Fully-Aware-h still collapses at
# L4_aggressive while Fully-Aware-CF-h holds ~nominal) -> the headline contrast
# is NOT a truncation artifact. The decision rule lives in aggregate.R/NOTES.md.
# =====================================================================

# Sourced AFTER codes/estimators.R (building blocks + SuperLearner/survey/
# surveyCV/tmle already loaded). We do not re-library here to avoid masking.

# ---- harmonized clamps (the ONE rule every arm uses) -------------------------
.q_clamp <- function(Q, q_lo)     pmin(pmax(Q, q_lo),    1 - q_lo)
.g_clamp <- function(g, g_floor)  pmin(pmax(g, g_floor), 1 - g_floor)

# ---- one targeting step: pre-fit (Q, g1) -> tmle(Q=, g1W=) -> harmonized EIF --
# obsw is the ARM's targeting weight (survey w for arms 1/2/4/5; rep(1,n) for 3)
# and is also the EIF weight handed to .eif_from_tmle (w_arm).
# DIVERGENCE GUARD (house pattern, R03/R05): weighted nuisance fits can separate
# and the fluctuation then diverges (eps -> ~1e13), yielding a meaningless point
# estimate; such reps are marked NON-CONVERGENT (b/se = NA) rather than emitting
# a garbage number. aggregate.R counts them (n_diverged) and drops them.
.target_h <- function(d, W, Q, g1, obsw, g_floor) {
  fit <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial",
              obsWeights = obsw, Q = Q, g1W = g1)
  e   <- .eif_from_tmle(fit, d$Y, d$A, obsw, gbound = g_floor)
  eps_max  <- max(abs(fit$epsilon))
  diverged <- !is.finite(eps_max) || eps_max > 20
  list(e = e, eps_max = eps_max, diverged = diverged)
}

# =====================================================================
# run_harmonized(obs, learners, g_floor=0.05, q_lo=1e-3, V_cf=5L,
#                inner_cv_folds=5L)
#   Returns the SAME shape as run_estimators(): list(results, diagnostics).
#   results: data.frame(method, b, se, df, diverged) with the five arm labels
#     (Fully-Aware-CV-h omitted on single-learner rungs, like the engine).
#   diagnostics: list(eif_fah = the Fully-Aware-h EIF (for the deff_clust /
#     ICC audit, the analogue of run_sim.R's eif_fa), drow = per-arm
#     eps / g-range / V_eff table).
# =====================================================================
run_harmonized <- function(obs, learners = "SL.glm", g_floor = 0.05, q_lo = 1e-3,
                           V_cf = 5L, inner_cv_folds = 5L) {
  d    <- obs                               # whole-sample (no domain in the sim)
  n    <- nrow(d); w <- d$weight
  W    <- d[, c("L1", "L2", "L3", "L4"), drop = FALSE]
  XA   <- d[, c("A", "L1", "L2", "L3", "L4"), drop = FALSE]
  Xa1  <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)
  ones <- rep(1, n)

  row_of <- function(label, t, se, df) data.frame(
    method = label,
    b  = if (t$diverged) NA_real_ else t$e$psi,
    se = if (t$diverged) NA_real_ else se,
    df = df, diverged = t$diverged, stringsAsFactors = FALSE)
  drow_of <- function(label, t, g1, V_eff = NA_integer_) data.frame(
    method = label, eps = t$eps_max, g_min = min(g1), g_max = max(g1),
    V_eff = V_eff, stringsAsFactors = FALSE)
  res <- list(); dr <- list()

  ## ---- arms 1 + 2: WEIGHTED single-fit nuisances, shared fit/EIF -------------
  ## (the harmonized analogue of the engine's Fully-Aware/Partially-Aware pair;
  ##  identical to R03's SF-W cell but emitting BOTH design SEs)
  qf_w <- .sl(d$Y, XA, weights = w, learners)      # Q on (A, W), weighted
  gf_w <- .sl(d$A, W,  weights = w, learners)      # g on W, weighted
  Q_sf <- .q_clamp(cbind(qf_w(Xa0), qf_w(Xa1)), q_lo)
  g_sf <- .g_clamp(gf_w(W), g_floor)               # floored BEFORE targeting
  t_fa <- .target_h(d, W, Q_sf, g_sf, w, g_floor)  # weighted targeting
  s_fa <- .se_des(t_fa$e$eif, d$strata, d$cluster, w, clustered = TRUE)
  s_pa <- .se_des(t_fa$e$eif, d$strata, d$cluster, w, clustered = FALSE)
  res$fa <- row_of("Fully-Aware-h",     t_fa, s_fa$se, s_fa$df)
  res$pa <- row_of("Partially-Aware-h", t_fa, s_pa$se, s_pa$df)   # same psi/EIF
  dr$fa  <- drow_of("Fully-Aware-h",     t_fa, g_sf)
  dr$pa  <- drow_of("Partially-Aware-h", t_fa, g_sf)

  ## ---- arm 3: Non-Aware-h (UNWEIGHTED nuisances + UNWEIGHTED targeting) ------
  qf_u <- .sl(d$Y, XA, weights = NULL, learners)
  gf_u <- .sl(d$A, W,  weights = NULL, learners)
  Q_na <- .q_clamp(cbind(qf_u(Xa0), qf_u(Xa1)), q_lo)
  g_na <- .g_clamp(gf_u(W), g_floor)
  t_na <- .target_h(d, W, Q_na, g_na, ones, g_floor)   # obsWeights = 1, EIF w = 1
  res$na <- row_of("Non-Aware-h", t_na, sqrt(var(t_na$e$eif) / n), n - 1L)
  dr$na  <- drow_of("Non-Aware-h", t_na, g_na)

  ## ---- arm 4: Fully-Aware-CV-h (multi-learner rungs only, like the engine) ---
  ## Mirrors estimators.R lines ~160-182 exactly (surveyCV::folds.svy validRows,
  ## WEIGHTED internal-CV nuisances = the sim-path foil), with the engine's
  ## hardcoded 0.05/1e-3 replaced by the harmonized g_floor/q_lo.
  if (length(learners) > 1) {
    vr <- tryCatch(split(seq_len(n),
                         surveyCV::folds.svy(d, nfolds = inner_cv_folds, clusterID = "cluster")),
                   error = function(e) NULL)
    cvc <- if (is.null(vr)) list() else list(V = length(vr), validRows = vr)
    qg <- SuperLearner(d$Y, XA, family = binomial(), SL.library = learners,
                       obsWeights = w, cvControl = cvc)
    gg <- SuperLearner(d$A, W,  family = binomial(), SL.library = learners,
                       obsWeights = w, cvControl = cvc)
    Q_cv <- .q_clamp(cbind(as.numeric(predict(qg, Xa0)$pred),
                           as.numeric(predict(qg, Xa1)$pred)), q_lo)
    g_cv <- .g_clamp(as.numeric(predict(gg, W)$pred), g_floor)
    t_cv <- .target_h(d, W, Q_cv, g_cv, w, g_floor)
    s_cv <- .se_des(t_cv$e$eif, d$strata, d$cluster, w, clustered = TRUE)
    res$cv <- row_of("Fully-Aware-CV-h", t_cv, s_cv$se, s_cv$df)
    dr$cv  <- drow_of("Fully-Aware-CV-h", t_cv, g_cv)
  }

  ## ---- arm 5: Fully-Aware-CF-h (engine CF logic; UNWEIGHTED per-fold fits) ---
  fold <- make_cf_folds(d$strata, d$cluster, V_cf)
  Q0o <- Q1o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = NULL, learners)   # UNWEIGHTED (ignorability)
    gf <- .sl(d$A[tr], W[tr, ],  weights = NULL, learners)
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g1o  <- .g_clamp(g1o, g_floor)                 # OOF g floored at g_floor
  Qoo  <- .q_clamp(cbind(Q0o, Q1o), q_lo)
  t_cf <- .target_h(d, W, Qoo, g1o, w, g_floor)  # weighted POOLED targeting
  s_cf <- .se_des(t_cf$e$eif, d$strata, d$cluster, w, clustered = TRUE)
  res$cf <- row_of("Fully-Aware-CF-h", t_cf, s_cf$se, s_cf$df)
  dr$cf  <- drow_of("Fully-Aware-CF-h", t_cf, g1o, V_eff = attr(fold, "V_eff"))

  results <- do.call(rbind, res); rownames(results) <- NULL
  drow    <- do.call(rbind, dr);  rownames(drow)    <- NULL
  list(results = results,
       diagnostics = list(eif_fah = t_fa$e$eif, drow = drow))
}
