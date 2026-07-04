# =====================================================================
# cvu_arm.R  —  R18_cvu helper: UNWEIGHTED-nuisance internal-CV arm (CV-u)
#
# SPEC A4 (Writing/comments/phase4-arc-sim-specs.md): the APPLICATION's CV arm
# fits its nuisances UNWEIGHTED (the real-data path in run_estimators(), taken
# when `inpop` is supplied; see codes/estimators.R header note on arm 4), but
# the locked SIMULATION CV foil fits them WEIGHTED. This helper adds the
# unweighted variant to the simulation so "internal CV != cross-fitting" is
# demonstrated for the object the application actually reports.
#
# run_cvu() mirrors the engine's Fully-Aware-CV block
# (codes/estimators.R, lines ~160-182) EXACTLY, with ONE change:
#   * obsWeights for BOTH nuisance SuperLearner fits = rep(1, n)  <- real-data path
# Everything downstream is IDENTICAL to the engine CV arm:
#   * the same surveyCV::folds.svy cluster-aware internal-CV folds,
#   * the same Q truncation to [1e-3, 1-1e-3] and g floor at g_oof_bound (0.05),
#   * the same WEIGHTED tmle() targeting (obsWeights = w, Q + g1W supplied),
#   * the same EIF reconstruction .eif_from_tmle(..., gbound = 1e-3) — the
#     engine DEFAULT, kept so rows are comparable to the locked CV-w/CF rows,
#   * the same clustered design SE via .se_des (Eq 8) with its design df.
# Method label: "Fully-Aware-CVu". Valid here because the simulation's design
# is ignorable given C (S _||_ (A,Y) | C: selection depends only on the
# stratum/PSU, and A,Y depend on the stratum only through C), so the
# unweighted nuisance fits are consistent — any CVu under-coverage at L3 is
# attributable to the LACK OF SAMPLE-SPLITTING, not to weighting bias.
#
# REUSES (read-only, from codes/estimators.R): .eif_from_tmle, .se_des, and the
# packages it loads (SuperLearner, survey, surveyCV, tmle). NEVER edits engine
# files; sourced only by codes/arc_runs/R18_cvu/run.R.
# =====================================================================

run_cvu <- function(obs, learners, inner_cv_folds = 5L, g_oof_bound = 0.05) {
  stopifnot(length(learners) > 1L)   # CV arm exists only on multi-learner rungs (L2/L3)
  d   <- obs
  n   <- nrow(d); w <- d$weight
  W_cols <- c("L1", "L2", "L3", "L4")               # simulation covariates
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)

  ## ---- cluster-aware internal-CV folds (identical to the engine CV block) ----
  vr <- tryCatch(split(seq_len(n),
                       surveyCV::folds.svy(d, nfolds = inner_cv_folds, clusterID = "cluster")),
                 error = function(e) NULL)
  cvc <- if (is.null(vr)) list() else list(V = length(vr), validRows = vr)

  ## ---- nuisance fits: UNWEIGHTED — the ONLY change vs codes/estimators.R ----
  w_cvu <- rep(1, n)                                 # real-data (application) path
  qg <- SuperLearner(d$Y, XA, family = binomial(), SL.library = learners,
                     obsWeights = w_cvu, cvControl = cvc)
  gg <- SuperLearner(d$A, W,  family = binomial(), SL.library = learners,
                     obsWeights = w_cvu, cvControl = cvc)
  Qcv <- cbind(as.numeric(predict(qg, Xa0)$pred), as.numeric(predict(qg, Xa1)$pred))
  gcv <- pmin(pmax(as.numeric(predict(gg, W)$pred), g_oof_bound), 1 - g_oof_bound)

  ## ---- WEIGHTED targeting + EIF + clustered design SE (identical to engine) ----
  cv <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
             Q = pmin(pmax(Qcv, 1e-3), 1 - 1e-3), g1W = gcv)
  e_cv <- .eif_from_tmle(cv, d$Y, d$A, w, gbound = 1e-3)   # engine default gbound
  s_cv <- .se_des(e_cv$eif, d$strata, d$cluster, w, clustered = TRUE)

  list(
    results = data.frame(method = "Fully-Aware-CVu",
                         b = e_cv$psi, se = s_cv$se, df = s_cv$df),
    diagnostics = list(
      eif_cvu = e_cv$eif,
      drow = data.frame(
        eps_cvu   = max(abs(cv$epsilon)),            # targeting fluctuation size
        g_cvu_min = min(gcv), g_cvu_max = max(gcv),  # floored propensity range
        Q_cvu_min = min(Qcv), Q_cvu_max = max(Qcv),  # pre-truncation Q range
        cv_V      = if (is.null(vr)) NA_integer_ else length(vr)))
  )
}
