# =====================================================================
# r21_helpers.R  —  R21_deployable_cvcf run-local helpers (NON-DESTRUCTIVE)
#
# Two new pieces, both defined HERE (never edits engine code):
#
#   (1) SL.xgboost.tuned  — a tuned-xgboost SuperLearner wrapper that carries
#       the "aggressive-but-DEPLOYABLE" load of the L6 library. It sits next to
#       the commented SL.xgb.deep slot at learners.R:48-51 in spirit, but is
#       deliberately NOT interpolating-by-construction (max_depth modest,
#       shrinkage 0.05, minobspernode>=5) — the point of R21 is that L6 is the
#       library one would ACTUALLY deploy, not a strawman lone interpolator.
#
#   (2) one_arm_cv() / aipw_cv_arm() — the new AIPW cluster-CV arm: the AIPW
#       counterpart of the engine's cluster-aware Fully-Aware-CV arm
#       (estimators.R:177-199). It fits Q,g with ONE SuperLearner call whose
#       cvControl uses cluster-aware (whole-PSU) folds via
#       surveyCV::folds.svy(d, clusterID='cluster'), predicts IN-SAMPLE, clamps
#       g/Q with the R15 clamps VERBATIM, then calls R15's one_arm() machinery
#       (re-implemented here against aipw_helpers.R's .aipw_jkn + the engine's
#       .se_des) for the Hajek point + JKn replicate SE + Eq-8 linearized SE.
#
# ⚠ WEIGHTING-CONVENTION GUARD (critical for an unconfounded CV-vs-CF read):
#   The engine's CV arm in the SIMULATION path is WEIGHTED ("as a deliberate
#   foil", estimators.R:186-187: w_cv <- if (is.null(inpop)) w else rep(1,n)).
#   The AIPW-CV arm MUST match that convention or the AIPW SF/CV/CF comparison
#   would differ in the weighting protocol rather than on the CV-vs-CF axis.
#   So aipw_cv_arm() fits its CV nuisances WEIGHTED, using R15's mean-1 weight
#   normalization (w/mean(w)) exactly as AIPW-SF does (aipw_helpers.R:132) —
#   raw survey weights spuriously diverge glm's IRLS; tmle()/SuperLearner apply
#   the same internal normalization. The Hajek point estimate and BOTH survey
#   variances keep the RAW weights, identical to AIPW-SF. Net: AIPW-CV is the
#   weighted-single-SuperLearner-fit-with-cluster-CV-folds analogue of the
#   engine's weighted Fully-Aware-CV foil — same weighting, only CV-folds vs
#   single-fit vs cross-fit differs across the three AIPW arms.
#
# REUSES (read-only, sourced by run.R BEFORE this file):
#   codes/estimators.R   -> .sl(), .se_des(), make_cf_folds()
#   R15/aipw_helpers.R   -> .aipw_jkn(), .AIPW_Q_CLAMP, aipw_arms()
#   codes/dgp.R, config.R, diagnostics.R, learners.R
# =====================================================================

suppressMessages({ library(SuperLearner) })

# ---------------------------------------------------------------------
# (1) Tuned-but-deployable xgboost wrapper.
#   ntrees=500, shrinkage=0.05, minobspernode=5; max_depth fixed at 6 by
#   default (cheaper + reproducible), or genuinely tuned over {3,6,8} via a
#   tiny internal grid when R21_XGB_TUNE=1. Both modes are "aggressive but not
#   interpolating": minobspernode>=5 and shrinkage 0.05 keep it off the
#   interpolation boundary that the lone deep RF sits on. num.threads is forced
#   to 1 because parLapply already parallelizes over reps (same discipline as
#   SL.ranger.deep at learners.R:26).
# ---------------------------------------------------------------------
.R21_XGB_MAXDEPTH <- local({
  d <- Sys.getenv("R21_XGB_MAXDEPTH")
  if (nzchar(d)) as.integer(d) else 6L            # fixed config default
})
.R21_XGB_TUNE <- identical(Sys.getenv("R21_XGB_TUNE"), "1")

SL.xgboost.tuned <- function(Y, X, newX, family, obsWeights, ...) {
  if (!.R21_XGB_TUNE) {
    return(SuperLearner::SL.xgboost(
      Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
      ntrees = 500, max_depth = .R21_XGB_MAXDEPTH, shrinkage = 0.05,
      minobspernode = 5, nthread = 1))
  }
  # genuine (cheap) depth tuning over {3,6,8} via xgboost's internal CV-free
  # holdout: pick the depth minimizing a 5-fold xgb.cv logloss/RMSE. Kept small
  # so it stays "deployable" without exploding wall-time; OFF by default.
  depths <- c(3L, 6L, 8L)
  best <- depths[1]; best_err <- Inf
  for (md in depths) {
    fit <- tryCatch(SuperLearner::SL.xgboost(
      Y = Y, X = X, newX = X, family = family, obsWeights = obsWeights,
      ntrees = 500, max_depth = md, shrinkage = 0.05,
      minobspernode = 5, nthread = 1), error = function(e) NULL)
    if (is.null(fit)) next
    err <- mean((as.numeric(fit$pred) - Y)^2)        # in-sample proxy (cheap)
    if (is.finite(err) && err < best_err) { best_err <- err; best <- md }
  }
  SuperLearner::SL.xgboost(
    Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
    ntrees = 500, max_depth = best, shrinkage = 0.05,
    minobspernode = 5, nthread = 1)
}

# ---------------------------------------------------------------------
# (1b) SL.hal9001 wrapper. SuperLearner does NOT ship one and hal9001 does NOT
#   export one (passing the string "SL.hal9001" errors with "object not found"),
#   so define it HERE -- sourced on every parallel worker via run.R's
#   clusterEvalQ. Namespaced to hal9001::fit_hal, so hal9001 only needs to be
#   INSTALLED, not attached. The basis is CAPPED by default (max_degree=1
#   additive, num_knots=5 per covariate) so HAL stays tractable inside the
#   SuperLearner CV loop; raise R21_HAL_DEGREE / R21_HAL_KNOTS to enrich it.
# ---------------------------------------------------------------------
.R21_HAL_DEGREE <- local({ d <- Sys.getenv("R21_HAL_DEGREE"); if (nzchar(d)) as.integer(d) else 1L })
.R21_HAL_KNOTS  <- local({ k <- Sys.getenv("R21_HAL_KNOTS");  if (nzchar(k)) as.integer(k) else 5L })

SL.hal9001 <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)),
                       max_degree = .R21_HAL_DEGREE, smoothness_orders = 1L,
                       num_knots = NULL, ...) {
  if (is.null(num_knots)) num_knots <- rep(.R21_HAL_KNOTS, max_degree)
  fam <- if (family$family == "binomial") "binomial" else "gaussian"
  fit <- hal9001::fit_hal(X = as.matrix(X), Y = as.numeric(Y), family = fam,
                          max_degree = max_degree, smoothness_orders = smoothness_orders,
                          num_knots = num_knots, weights = obsWeights)
  pred <- as.numeric(stats::predict(fit, new_data = as.matrix(newX)))
  out  <- list(pred = pred, fit = list(object = fit))
  class(out$fit) <- "SL.hal9001"
  out
}
predict.SL.hal9001 <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, new_data = as.matrix(newdata)))
}

# packages the L6 library + AIPW-CV arm need on each parallel worker
R21_PKGS <- c("SuperLearner", "earth", "glmnet", "xgboost", "ranger",
              "survey", "surveyCV", "tmle", "hal9001")

# ---------------------------------------------------------------------
# (2) AIPW cluster-CV arm.
#   aipw_cv_arm(obs, learners, g_floor, inner_cv_folds, nest, inpop)
#   Returns a one-row results data.frame in the SAME schema as
#   aipw_helpers.R::aipw_arms (method, b, se_jkn, se_lin, df) plus the centered
#   pseudo-outcome + diagnostics, so run.R can rbind it onto the AIPW-SF/-CF
#   rows from aipw_arms().
#
#   Mirrors estimators.R:178-198 on the AIPW side:
#     - validRows = split(seq_len(n), surveyCV::folds.svy(d, nfolds, 'cluster'))
#       -> ONE weighted SuperLearner fit per nuisance whose internal CV honors
#          the cluster structure (whole-PSU folds), exactly like the engine arm.
#     - in-sample predict for Q1/Q0/g  (engine does predict(qg, Xa1) etc.)
#     - clamp g -> [g_floor, 1-g_floor], Q -> [.AIPW_Q_CLAMP, 1-.AIPW_Q_CLAMP]
#       (R15 clamps verbatim: aipw_helpers.R:148-154)
#     - Hajek AIPW point + .aipw_jkn() (JKn) + .se_des() (Eq-8) via the same
#       one_arm() logic as R15.
# ---------------------------------------------------------------------
aipw_cv_arm <- function(obs, learners, g_floor = 0.05, inner_cv_folds = 5L,
                        W_cols = c("L1", "L2", "L3", "L4"),
                        nest = FALSE, inpop = NULL) {
  sub <- if (is.null(inpop)) seq_len(nrow(obs)) else which(inpop)
  d   <- obs[sub, , drop = FALSE]
  n   <- nrow(d); w <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)

  ## ---- cluster-aware (whole-PSU) internal-CV folds, engine recipe ----------
  ## surveyCV::folds.svy needs a data.frame carrying the cluster column; d has
  ## `cluster`. tryCatch -> default folds if folds.svy errors (mirrors the
  ## engine's `vr <- tryCatch(...)` guard at estimators.R:179-181).
  vr  <- tryCatch(split(seq_len(n),
                        surveyCV::folds.svy(d, nfolds = inner_cv_folds,
                                            clusterID = "cluster")),
                  error = function(e) NULL)
  cvc <- if (is.null(vr)) list(V = inner_cv_folds) else list(V = length(vr), validRows = vr)

  ## ---- WEIGHTING GUARD: weighted CV fits (mean-1 normalized), matching the
  ## engine's weighted CV foil (estimators.R:187 sim path) AND R15's AIPW-SF
  ## mean-1 normalization (aipw_helpers.R:132). Point estimate + variances keep
  ## RAW weights below. -----------------------------------------------------
  w_fit <- w / mean(w)
  qg <- SuperLearner(Y = d$Y, X = XA, family = binomial(), SL.library = learners,
                     obsWeights = w_fit, cvControl = cvc)
  gg <- SuperLearner(Y = d$A, X = W,  family = binomial(), SL.library = learners,
                     obsWeights = w_fit, cvControl = cvc)
  Q1_raw <- as.numeric(predict(qg, Xa1)$pred)
  Q0_raw <- as.numeric(predict(qg, Xa0)$pred)
  g_raw  <- as.numeric(predict(gg, W)$pred)

  ## ---- R15 clamps VERBATIM ------------------------------------------------
  clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  g1 <- clamp(g_raw,  g_floor, 1 - g_floor)
  Q1 <- clamp(Q1_raw, .AIPW_Q_CLAMP, 1 - .AIPW_Q_CLAMP)
  Q0 <- clamp(Q0_raw, .AIPW_Q_CLAMP, 1 - .AIPW_Q_CLAMP)

  ## ---- Hajek AIPW point + JKn replicate SE + Eq-8 linearized SE -----------
  ## (the one_arm() body from aipw_helpers.R:157-172, raw weights throughout)
  D   <- d$A * (d$Y - Q1) / g1 - (1 - d$A) * (d$Y - Q0) / (1 - g1) + Q1 - Q0
  psi <- sum(w * D) / sum(w)
  jk <- if (is.null(inpop))
    .aipw_jkn(D, d$strata, d$cluster, w, nest = nest)
  else
    .aipw_jkn(D, obs$strata, obs$cluster, obs$weight, nest = nest, inpop = inpop)
  sl <- if (is.null(inpop))
    .se_des(D - psi, d$strata, d$cluster, w, clustered = TRUE)
  else
    .se_des(D - psi, obs$strata, obs$cluster, obs$weight, clustered = TRUE,
            nest = nest, inpop = inpop)

  row <- data.frame(method = "AIPW-CV", b = psi, se_jkn = jk$se,
                    se_lin = sl$se, df = sl$df, stringsAsFactors = FALSE)
  drow <- data.frame(
    g_cv_min = min(g1), g_cv_max = max(g1),
    g_cv_outside_floor = mean(g_raw < g_floor | g_raw > 1 - g_floor),
    cv_V_eff = if (is.null(vr)) NA_integer_ else length(vr),
    df_jkn_cv = jk$df, jkn_mode_cv = jk$mode,
    g_floor = g_floor, stringsAsFactors = FALSE)

  list(results = row, D_centered = D - psi, drow = drow)
}

# ---------------------------------------------------------------------
# aipw_arms_all(obs, learners, g_floor, V_cf, inner_cv_folds, ...)
#   Convenience wrapper that returns the THREE AIPW arms (SF, CF, CV) in one
#   results data.frame, schema-compatible with aipw_arms(). Reuses R15's
#   aipw_arms() for SF+CF (unchanged) and appends the new CV arm.
# ---------------------------------------------------------------------
aipw_arms_all <- function(obs, learners, g_floor = 0.05, V_cf = 5L,
                          inner_cv_folds = 5L,
                          W_cols = c("L1", "L2", "L3", "L4"),
                          nest = FALSE, inpop = NULL) {
  base <- aipw_arms(obs, learners = learners, g_floor = g_floor, V_cf = V_cf,
                    W_cols = W_cols, nest = nest, inpop = inpop)   # AIPW-SF, AIPW-CF
  cv   <- aipw_cv_arm(obs, learners = learners, g_floor = g_floor,
                      inner_cv_folds = inner_cv_folds, W_cols = W_cols,
                      nest = nest, inpop = inpop)                  # AIPW-CV (new)
  results <- rbind(base$results, cv$results)
  list(results = results,
       diagnostics = list(drow = base$diagnostics$drow, drow_cv = cv$drow,
                          D_sf = base$diagnostics$D_sf,
                          D_cf = base$diagnostics$D_cf,
                          D_cv = cv$D_centered))
}
