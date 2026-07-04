# =====================================================================
# estimators.R  —  the five survey-TMLE arms (tmle-package engine)
#
# Engine decision (Phase-2 smoke, 2026-06-06): the hand-rolled single-step
# fluctuation hit perfect separation (eps -> 1e14) under realistic informative
# weights. The mature `tmle` package is robust and is what the project's NHANES
# analysis already uses. We use tmle() for the (robust) TARGETING, then
# RECONSTRUCT the per-unit EIF from its targeted Qstar + g1W and apply Eq 8 via
# svymean -- avoiding the double-weighting that inflates the SE if IC.ATE
# (already weight-incorporated) is fed straight to svymean.
#
# Arms:
#   1. Fully-Aware     weighted tmle() fit + targeting; design SE (Eq 8)
#   2. Partially-Aware same psi/EIF; SE = WEIGHTED design SE WITHOUT clustering
#   3. Non-Aware       unweighted tmle(); iid SE
#   4. Fully-Aware-CV  cluster-aware INTERNAL CV, full-sample nuisance fit -> tmle
#                      targeting; design SE  (NOT cross-fitted; SL tier only).
#                      Nuisances fit WEIGHTED on the simulation path (a deliberate
#                      foil) but UNWEIGHTED on the real-data domain path (inpop
#                      supplied) -- identical to the CF arm, valid under
#                      S _||_ (A,Y) | C -- so on real data CV and CF differ ONLY
#                      in internal-CV vs PSU-level cross-fitting.
#   5. Fully-Aware-CF  PSU-level cross-fitted: UNWEIGHTED per-fold nuisance fits
#                      (consistent by sampling-ignorability; avoids per-fold
#                      separation) -> tmle() weighted POOLED targeting; design SE
#
# Data contract (draw_sample): L1..L4, A, Y, strata, cluster (globally-unique
# PSU id), weight.  svydesign uses ids=~cluster (no nest needed).
#
# Generality (added 2026-06-07 for the NHANES application): run_estimators() also
# takes W_cols (arbitrary covariate names; default the simulation's L1..L4), nest
# (TRUE when PSU ids repeat across strata, as in NHANES SDMVPSU), and inpop (a
# logical over rows of `obs` selecting a SUB-POPULATION/domain). With inpop given,
# nuisances + targeting are fit on the sub-population, but the design variance is
# computed on the FULL design then subset()-ed to the domain -- i.e. subset the
# DESIGN, not the data, preserving the full stratum/PSU structure so PSUs that
# partially intersect the domain contribute correctly to the domain variance
# (EpiMethods surveydata8). NOTE: survey::degf excludes PSUs with no in-domain
# observations (subset() zero-weights them), so a wholly-absent PSU does NOT add
# df; in the four NHANES examples every PSU intersects every domain, so the
# domain df equal the full-design df. With inpop = NULL (the
# default) the behaviour is byte-identical to the locked simulation.
# =====================================================================

suppressMessages({
  library(SuperLearner); library(survey); library(surveyCV); library(tmle)
})

# ---- reconstruct the per-unit EIF from a tmle fit's targeted nuisances -------
# Uses the fit's OWN g-bounds (fit$g$bound) and replicates tmle's internal
# bounding exactly: fit$g$g1W is stored RAW; the fluctuation floored A's
# denominator at .bound(g1W, gbound) and (1-A)'s at .bound(1-g1W, gbound)
# INDEPENDENTLY. Reconstructing with any other truncation (the old hard-coded
# 1e-3) yields an EIF the targeting step did not solve, exactly where raw g
# escapes the bound (aggressive learners / limited overlap).
.eif_from_tmle <- function(fit, Y, A, w, gbound = NULL) {
  Qst <- fit$Qstar                          # n x 2 : [Q*(A=0,W), Q*(A=1,W)]
  if (is.null(gbound)) {                    # DEFAULT: replicate the fit's own bounding
    gb <- fit$g$bound
    if (is.null(gb)) gb <- c(1e-3, 1 - 1e-3)      # defensive: fits lacking $bound
    .tb <- function(x) pmin(pmax(x, min(gb)), max(gb))
    g1den <- .tb(fit$g$g1W)                 # denominator for the A=1 score
    g0den <- .tb(1 - fit$g$g1W)             # denominator for the A=0 score
  } else {                                  # LEGACY: explicit symmetric truncation
    g1 <- pmin(pmax(fit$g$g1W, gbound), 1 - gbound)   # (kept for R14/R16/R18 arc_runs
    g1den <- g1; g0den <- 1 - g1                      #  that pass gbound explicitly)
  }
  Q0s <- Qst[, 1]; Q1s <- Qst[, 2]
  QAs <- ifelse(A == 1, Q1s, Q0s)
  H   <- A / g1den - (1 - A) / g0den
  psi <- weighted.mean(Q1s - Q0s, w)
  # ---- per-arm decomposition (for RR/OR + continuous; additive, non-breaking) --
  # psi_a = weighted mean of Q*(a,.) ; D_a = the clever-covariate score for the
  # marginal mean E[Y^a]. By construction D1 - D0 equals the contrast EIF returned
  # below and psi1 - psi0 equals `psi`, so the scalar `psi`/`eif` are byte-identical
  # to before -- callers that read only $psi/$eif are unaffected. The RR/OR post-
  # processor (codes/estimands_rr_or.R) consumes psi1/psi0/eif1/eif0.
  psi1 <- weighted.mean(Q1s, w); psi0 <- weighted.mean(Q0s, w)
  eif1 <- (A / g1den)       * (Y - Q1s) + Q1s - psi1
  eif0 <- ((1 - A) / g0den) * (Y - Q0s) + Q0s - psi0
  list(psi = psi, eif = H * (Y - QAs) + (Q1s - Q0s) - psi,
       psi1 = psi1, psi0 = psi0, eif1 = eif1, eif0 = eif0)
}

# ---- design SE of an EIF: Eq 8 (clustered) or weighted-no-cluster ------------
# Returns the SE and the design degrees of freedom (for a t-reference CI; the
# clustered design has df ~ #PSU - #strata, the unclustered design df ~ n-1).
.se_des <- function(eif, strata, cluster, weight, clustered = TRUE, nest = FALSE, inpop = NULL) {
  if (is.null(inpop)) {                       # simulation path -- unchanged
    dat <- data.frame(eif = eif, strata = strata, cluster = cluster, weight = weight)
    des <- if (clustered)
      svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dat)
    else
      svydesign(ids = ~1, weights = ~weight, data = dat)
    return(list(se = as.numeric(SE(svymean(~eif, des))), df = survey::degf(des)))
  }
  # sub-population (domain) path: build the design on the FULL sample, then
  # subset() it to the domain (EpiMethods surveydata8: subset the DESIGN, not the
  # data) so partially-covered PSUs contribute correctly to the domain variance.
  # NOTE: survey::degf drops PSUs with no in-domain rows, so wholly-absent PSUs
  # do NOT add df (none occur in the four NHANES example domains).
  # `eif` is over the in-pop rows; strata/cluster/weight/inpop are full-length.
  eif_full <- numeric(length(inpop)); eif_full[inpop] <- eif
  dat <- data.frame(eif = eif_full, strata = strata, cluster = cluster,
                    weight = weight, .inpop = inpop)
  des <- if (clustered)
    svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dat, nest = nest)
  else
    svydesign(ids = ~1, weights = ~weight, data = dat)
  des <- subset(des, .inpop)
  list(se = as.numeric(SE(svymean(~eif, des))), df = survey::degf(des))
}

# ---- cross-fit folds: whole PSUs within strata; V reduced to min PSUs/stratum
# so EVERY training fold still spans EVERY stratum (no cross-stratum
# extrapolation). For R1 (2 PSUs/stratum) this is leave-one-PSU-per-stratum-out.
make_cf_folds <- function(strata, cluster, V) {
  # key folds on the (stratum, PSU) PAIR, so this is correct whether the PSU id
  # is globally unique (simulation) or repeats across strata (NHANES SDMVPSU).
  # When the id is already globally unique this is identical to the original.
  ucl <- paste(strata, cluster, sep = "\r")
  key <- unique(data.frame(strata = strata, cluster = cluster, ucl = ucl,
                           stringsAsFactors = FALSE))
  key <- key[order(key$strata), ]
  min_psu <- min(tapply(key$ucl, key$strata, length))
  V_eff <- max(2L, min(as.integer(V), as.integer(min_psu)))
  key$fold <- ave(seq_len(nrow(key)), key$strata,
                  FUN = function(ix) sample(rep_len(1:V_eff, length(ix))))
  fold <- key$fold[match(ucl, key$ucl)]
  attr(fold, "unit") <- "psu"; attr(fold, "V_eff") <- V_eff
  fold
}

# ---- helper: SuperLearner fit returning a predictor closure -----------------
# For a SINGLE-learner library the internal CV only weights one learner (weight=1
# regardless of V) and the final fit is on all data either way, so V=2 is
# statistically identical to V=10 but ~5x cheaper -- a free speedup on the
# single-learner rungs (L1, L4). Multi-learner libraries keep the default CV.
.sl <- function(Y, X, weights = NULL, learners, family = binomial()) {
  cvc <- if (length(learners) == 1L) list(V = 2L) else list()
  fit <- SuperLearner(Y = Y, X = X, family = family, SL.library = learners,
                      obsWeights = if (is.null(weights)) rep(1, length(Y)) else weights,
                      cvControl = cvc)
  function(newX) as.numeric(predict(fit, newdata = newX)$pred)
}

# =============================== driver =====================================
# learners: SuperLearner library vector (e.g. "SL.glm" or c("SL.glm","SL.gam",...)).
run_estimators <- function(obs, learners = "SL.glm", V_cf = 5L, inner_cv_folds = 5L,
                           g_oof_bound = 0.05, W_cols = c("L1", "L2", "L3", "L4"),
                           nest = FALSE, inpop = NULL, family = "binomial") {
  # inpop = NULL -> whole sample (simulation default, behaviour unchanged).
  # inpop = logical over rows of `obs` -> a sub-population/domain: nuisances and
  # targeting are fit on the domain `d`, while the design variance is computed on
  # the FULL design subset to the domain (see .se_des / header).
  sub <- if (is.null(inpop)) seq_len(nrow(obs)) else which(inpop)
  d   <- obs[sub, , drop = FALSE]
  n   <- nrow(d); w <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)
  # Outcome family: 'binomial' (binary, default) or 'gaussian' (continuous, the
  # bounded-[0,1] device of as:bounded -- tmle() scales Y internally and applies
  # the SAME logistic fluctuation). The g (treatment) model is ALWAYS binomial.
  # Qfam = the Q-side SuperLearner family; qbound = the manual Q-truncation used
  # by the CV/CF arms, skipped for gaussian (tmle bounds Q to the Y range).
  Qfam   <- if (identical(family, "gaussian")) gaussian() else binomial()
  qbound <- if (identical(family, "gaussian")) function(Q) Q else function(Q) pmin(pmax(Q, 1e-3), 1 - 1e-3)
  # variance closure: simulation path uses the domain rows directly; domain path
  # builds the FULL design and subset()s it to inpop.
  se_des <- function(eif, clustered)
    if (is.null(inpop)) .se_des(eif, d$strata, d$cluster, w, clustered)
    else .se_des(eif, obs$strata, obs$cluster, obs$weight, clustered, nest = nest, inpop = inpop)
  res <- list(); diag <- list(); .tim <- numeric(0)

  ## ---- Fully-Aware + Partially-Aware (shared weighted tmle fit) ----
  .t0 <- proc.time()[3]
  fa <- tmle(Y = d$Y, A = d$A, W = W, family = family, obsWeights = w,
             Q.SL.library = learners, g.SL.library = learners,
             cvQinit = FALSE)   # TRUE single fit: tmle >= 2.x defaults cvQinit=TRUE
                                # (cross-validated initial Q), which would make this
                                # arm CV on the outcome side -- not the single-fit
                                # comparator the paper describes.
  .tim["single_fit_weighted"] <- proc.time()[3] - .t0
  e_fa <- .eif_from_tmle(fa, d$Y, d$A, w)
  s_fa <- se_des(e_fa$eif, TRUE)
  s_pa <- se_des(e_fa$eif, FALSE)
  res$FullyAware     <- data.frame(method = "Fully-Aware",     b = e_fa$psi, se = s_fa$se, df = s_fa$df)
  res$PartiallyAware <- data.frame(method = "Partially-Aware", b = e_fa$psi, se = s_pa$se, df = s_pa$df)

  ## ---- Non-Aware (unweighted) ----
  .t0 <- proc.time()[3]
  na <- tmle(Y = d$Y, A = d$A, W = W, family = family,
             Q.SL.library = learners, g.SL.library = learners,
             cvQinit = FALSE)   # TRUE single fit (see Fully-Aware note above)
  .tim["single_fit_unweighted"] <- proc.time()[3] - .t0
  e_na <- .eif_from_tmle(na, d$Y, d$A, rep(1, n))
  res$NonAware <- data.frame(method = "Non-Aware", b = e_na$psi,
                             se = sqrt(var(e_na$eif) / n), df = n - 1)

  ## ---- Fully-Aware-CV (cluster-aware internal CV; SL tier only) ----
  if (length(learners) > 1) {
    vr <- tryCatch(split(seq_len(n),
                         surveyCV::folds.svy(d, nfolds = inner_cv_folds, clusterID = "cluster")),
                   error = function(e) NULL)
    cvc <- if (is.null(vr)) list() else list(V = length(vr), validRows = vr)
    # Real-data domain (inpop supplied): fit CV-arm nuisances UNWEIGHTED, exactly
    # like the primary CF arm (valid under S _||_ (A,Y) | C), so CV and CF differ
    # only in internal-CV vs cross-fitting and the weighted full-sample fit cannot
    # separate (cf. the E2/L2 sign-flip). Simulation path keeps weighted CV as a foil.
    w_cv <- if (is.null(inpop)) w else rep(1, n)
    .t0 <- proc.time()[3]
    qg <- SuperLearner(d$Y, XA, family = Qfam, SL.library = learners,
                       obsWeights = w_cv, cvControl = cvc)
    gg <- SuperLearner(d$A, W,  family = binomial(), SL.library = learners,
                       obsWeights = w_cv, cvControl = cvc)
    Qcv <- cbind(as.numeric(predict(qg, Xa0)$pred), as.numeric(predict(qg, Xa1)$pred))
    gcv <- pmin(pmax(as.numeric(predict(gg, W)$pred), g_oof_bound), 1 - g_oof_bound)
    cv <- tmle(Y = d$Y, A = d$A, W = W, family = family, obsWeights = w,
               Q = qbound(Qcv), g1W = gcv)
    e_cv <- .eif_from_tmle(cv, d$Y, d$A, w)
    s_cv <- se_des(e_cv$eif, TRUE)
    res$FullyAwareCV <- data.frame(method = "Fully-Aware-CV", b = e_cv$psi, se = s_cv$se, df = s_cv$df)
    .tim["cv"] <- proc.time()[3] - .t0
  }

  ## ---- Fully-Aware-CF (cross-fitted; UNWEIGHTED per-fold fits) ----
  .t0 <- proc.time()[3]
  fold <- make_cf_folds(d$strata, d$cluster, V_cf)
  Q0o <- Q1o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = NULL, learners, family = Qfam)        # UNWEIGHTED (ignorability)
    gf <- .sl(d$A[tr], W[tr, ],  weights = NULL, learners, family = binomial())
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g1o <- pmin(pmax(g1o, g_oof_bound), 1 - g_oof_bound)
  Qoo <- qbound(cbind(Q0o, Q1o))
  cf <- tmle(Y = d$Y, A = d$A, W = W, family = family, obsWeights = w,
             Q = Qoo, g1W = g1o)                                # weighted POOLED targeting
  e_cf <- .eif_from_tmle(cf, d$Y, d$A, w)
  s_cf <- se_des(e_cf$eif, TRUE)
  res$FullyAwareCF <- data.frame(method = "Fully-Aware-CF", b = e_cf$psi, se = s_cf$se, df = s_cf$df)
  .tim["cf"] <- proc.time()[3] - .t0

  # ---- per-rep diagnostics (kept for reviewer-proofing) ----
  g_fa <- fa$g$g1W
  drow <- data.frame(
    eps_fa   = max(abs(fa$epsilon)),         # FA targeting fluctuation
    g_fa_min = min(g_fa), g_fa_max = max(g_fa),
    eps_cf   = max(abs(cf$epsilon)),         # CF targeting fluctuation
    g_cf_min = min(g1o),  g_cf_max = max(g1o),   # out-of-fold propensity range
    cf_V_eff = attr(fold, "V_eff")
  )
  # g_fa / g_cf added (2026-06-07) for NHANES overlap diagnostics; the simulation
  # driver ignores them, so this is non-breaking.
  # `arms` surfaces the per-arm psi_hat + EIF components for the single-fit (FA)
  # and cross-fitted (CF) arms so the RR/OR delta-method post-processor can run;
  # additive, the simulation/NHANES drivers ignore it unless RR/OR is requested.
  diag <- list(eif_fa = e_fa$eif, g_fa = fa$g$g1W, g_cf = g1o, drow = drow, timing = .tim,
               arms = list(
                 "Fully-Aware"    = list(psi1 = e_fa$psi1, psi0 = e_fa$psi0,
                                         eif1 = e_fa$eif1, eif0 = e_fa$eif0),
                 "Fully-Aware-CF" = list(psi1 = e_cf$psi1, psi0 = e_cf$psi0,
                                         eif1 = e_cf$eif1, eif0 = e_cf$eif0)))
  list(results = do.call(rbind, res), diagnostics = diag)
}
