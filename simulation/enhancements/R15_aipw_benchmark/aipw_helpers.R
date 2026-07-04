# =====================================================================
# aipw_helpers.R  —  R15_aipw_benchmark helper: survey-weighted AIPW arms
# (spec item A7 "External survey-causal benchmark" in
#  Writing/comments/phase4-arc-sim-specs.md)
#
# PURPOSE. The paper's five arms are the SAME estimator (survey-TMLE)
# design-ablated; A7 asks "compared to WHAT?". This helper implements an
# independent, literature-standard competitor: survey-weighted AIPW with a
# survey JKn (delete-one-PSU-within-stratum jackknife) replicate SE, run as
# a FULL arm across the ladder (L1-L4) x both designs, and on NHANES E1.
#
# TWO nuisance protocols (mirroring the paper's single-fit-vs-cross-fit axis):
#   AIPW-SF  WEIGHTED single-fit nuisances: .sl(weights = w/mean(w)) for Q on
#            (A,W) and g on W, predicted IN-SAMPLE -> the AIPW analogue of the
#            Fully-Aware single-fit arm. Weights are normalized to mean 1 for
#            the FITS only (scale-invariant estimand; tmle() normalizes its
#            obsWeights identically; raw-scale weights spuriously diverge
#            glm's IRLS -- see inline note).
#   AIPW-CF  UNWEIGHTED out-of-fold nuisances over make_cf_folds() -- the
#            SAME PSU-within-stratum cross-fit protocol as the engine's
#            Fully-Aware-CF arm (unweighted per-fold fits, valid under the
#            sampling-ignorability condition S _||_ (A,Y) | C).
# Both arms: propensity floored to [g_floor, 1 - g_floor] (default 0.05,
# matching the engine's OOF floor) and Q truncated to [1e-3, 1-1e-3]
# (the engine's Q clamp).
#
# POINT ESTIMATE (each arm): Hajek-weighted AIPW over the pseudo-outcome
#   D_i  = A_i (Y_i - Q1_i)/g_i - (1-A_i)(Y_i - Q0_i)/(1-g_i) + Q1_i - Q0_i
#   psi^ = sum(w_i D_i) / sum(w_i)
#
# VARIANCE, TWO WAYS per arm:
#   se_jkn  survey replicate variance: svydesign(ids=~cluster, strata=~strata,
#           weights=~weight) -> as.svrepdesign(type = "JKn") ->
#           withReplicates() of the Hajek ratio of D, with the NUISANCES
#           FROZEN (Q, g computed once on the full sample; only the weights
#           are perturbed across replicates). Freezing the nuisances is the
#           standard practice for replicate-variance AIPW: the replicate SE
#           is conditional on the fitted nuisances, exactly as the paper's
#           Eq-8 SE plugs the fitted EIF into the design linearization.
#   se_lin  the engine's Eq-8 design linearization .se_des() applied to the
#           CENTERED pseudo-outcome (D - psi^) -- apples-to-apples with the
#           five paper arms' EIF-based SE.
#
# DOMAIN (NHANES) path, inpop supplied: nuisances + point estimate are fit on
# the sub-population, while BOTH variances are computed on the FULL design:
#   se_lin  via .se_des(..., nest, inpop)  (the engine's subset-the-DESIGN
#           domain logic, byte-identical to run_estimators()).
#   se_jkn  the replicate design is built on the FULL design, then
#           subset(rep_des, inpop) -> withReplicates (the replicate-weight
#           analogue of subset-the-DESIGN). If the subset-of-svrepdesign path
#           errors, we fall back to the mathematically identical domain-zeroed
#           ratio on the full replicate design:
#           sum(wts * D * 1{inpop}) / sum(wts * 1{inpop}); the realized mode
#           is recorded in diagnostics$drow$jkn_mode.
#
# REUSES the canonical building blocks exposed by sourcing codes/estimators.R
# (read-only; this file NEVER edits engine code):
#   .sl(Y, X, weights, learners)         -> SuperLearner predictor closure
#   make_cf_folds(strata, cluster, V)    -> whole-PSU-within-stratum folds
#   .se_des(eif, strata, cluster, weight, clustered, nest, inpop) -> Eq-8 SE
# =====================================================================

# Q clamp used by the engine before targeting (estimators.R); keep identical.
.AIPW_Q_CLAMP <- 1e-3

# ---- JKn replicate SE of the Hajek AIPW mean (nuisances frozen) -------------
# D over the analysis rows (full sample, or the in-domain rows when inpop is
# given); strata/cluster/weight (and inpop) are FULL-length in the domain case,
# mirroring the .se_des() contract.
.aipw_jkn <- function(D, strata, cluster, weight, nest = FALSE, inpop = NULL) {
  if (is.null(inpop)) {                       # simulation / whole-sample path
    dat <- data.frame(D = D, strata = strata, cluster = cluster, weight = weight)
    des <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight,
                     data = dat, nest = nest)
    rep_des <- as.svrepdesign(des, type = "JKn")
    est <- withReplicates(rep_des, function(wts, dat) sum(wts * dat$D) / sum(wts))
    return(list(se = as.numeric(SE(est)), df = survey::degf(des), mode = "full"))
  }
  # domain path: replicate design on the FULL design, then subset to the domain
  D_full <- numeric(length(inpop)); D_full[inpop] <- D
  dat <- data.frame(D = D_full, strata = strata, cluster = cluster,
                    weight = weight, .inpop = inpop)
  des <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight,
                   data = dat, nest = nest)
  rep_des <- as.svrepdesign(des, type = "JKn")
  df_dom  <- survey::degf(subset(des, .inpop))   # engine df convention (.se_des)
  out <- tryCatch({
    sub_rep <- subset(rep_des, .inpop)
    est <- withReplicates(sub_rep, function(wts, dat) sum(wts * dat$D) / sum(wts))
    list(se = as.numeric(SE(est)), df = df_dom, mode = "subset_repdes")
  }, error = function(e) NULL)
  if (is.null(out)) {
    # fallback (documented in NOTES.md): domain-zeroed Hajek ratio on the FULL
    # replicate design -- algebraically the same domain estimator.
    est <- withReplicates(rep_des, function(wts, dat)
      sum(wts * dat$D * dat$.inpop) / sum(wts * dat$.inpop))
    out <- list(se = as.numeric(SE(est)), df = df_dom, mode = "domain_zeroed_full")
  }
  out
}

# =====================================================================
# aipw_arms(obs, learners, g_floor = 0.05, V_cf = 5L, ...)
#   Returns list(results, diagnostics) shaped like run_estimators():
#   results: data.frame(method in {AIPW-SF, AIPW-CF}, b, se_jkn, se_lin, df)
#     (df is the engine's design df from .se_des: #PSU - #strata, domain-aware)
#   diagnostics: drow (g ranges post-floor, pre-floor outside-floor share,
#     cf_V_eff, df_jkn/jkn_mode [SF arm], df_jkn_cf/jkn_mode_cf [CF arm],
#     g_floor) + the CENTERED pseudo-outcomes
#     D_sf / D_cf (for deff_clust audits, like run_estimators' eif_fa).
#   W_cols / nest / inpop generalize to the NHANES domain exactly as in
#   run_estimators() (defaults preserve the simulation behaviour).
# =====================================================================
aipw_arms <- function(obs, learners = "SL.glm", g_floor = 0.05, V_cf = 5L,
                      W_cols = c("L1", "L2", "L3", "L4"),
                      nest = FALSE, inpop = NULL) {
  sub <- if (is.null(inpop)) seq_len(nrow(obs)) else which(inpop)
  d   <- obs[sub, , drop = FALSE]
  n   <- nrow(d); w <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)

  ## ---- protocol (a): AIPW-SF -- WEIGHTED single-fit, in-sample predictions --
  ## Weights are NORMALIZED to mean 1 for the nuisance fits: the weighted MLE
  ## is invariant to the weight SCALE, but raw survey weights (1e1-1e5) make
  ## glm's IRLS diverge spuriously (coefs -> 1e15, fitted -> {0,1}; verified on
  ## E1 and on sim reps 3/20 during the smoke). tmle() does the SAME internal
  ## normalization (obsWeights/sum(obsWeights)*n), so this makes AIPW-SF the
  ## honest analogue of the locked weighted Fully-Aware arm. The Hajek point
  ## estimate and both survey variances keep the RAW weights.
  .t0 <- proc.time()[3]
  w_fit <- w / mean(w)
  qf_sf <- .sl(d$Y, XA, weights = w_fit, learners)
  gf_sf <- .sl(d$A, W,  weights = w_fit, learners)
  Q1_sf_raw <- qf_sf(Xa1); Q0_sf_raw <- qf_sf(Xa0); g_sf_raw <- gf_sf(W)
  .t_sf <- as.numeric(proc.time()[3] - .t0)

  ## ---- protocol (b): AIPW-CF -- UNWEIGHTED out-of-fold (engine CF protocol) -
  .t0 <- proc.time()[3]
  fold <- make_cf_folds(d$strata, d$cluster, V_cf)
  Q1_cf_raw <- Q0_cf_raw <- g_cf_raw <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = NULL, learners)   # UNWEIGHTED (ignorability)
    gf <- .sl(d$A[tr], W[tr, ],  weights = NULL, learners)
    Q1_cf_raw[ho] <- qf(Xa1[ho, ]); Q0_cf_raw[ho] <- qf(Xa0[ho, ])
    g_cf_raw[ho]  <- gf(W[ho, ])
  }
  .t_cf <- as.numeric(proc.time()[3] - .t0)

  clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  g_sf <- clamp(g_sf_raw, g_floor, 1 - g_floor)
  g_cf <- clamp(g_cf_raw, g_floor, 1 - g_floor)
  Q1s  <- clamp(Q1_sf_raw, .AIPW_Q_CLAMP, 1 - .AIPW_Q_CLAMP)
  Q0s  <- clamp(Q0_sf_raw, .AIPW_Q_CLAMP, 1 - .AIPW_Q_CLAMP)
  Q1c  <- clamp(Q1_cf_raw, .AIPW_Q_CLAMP, 1 - .AIPW_Q_CLAMP)
  Q0c  <- clamp(Q0_cf_raw, .AIPW_Q_CLAMP, 1 - .AIPW_Q_CLAMP)

  ## ---- one arm: Hajek AIPW point + JKn replicate SE + Eq-8 linearized SE ----
  one_arm <- function(Q1, Q0, g1, label) {
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
    list(row = data.frame(method = label, b = psi, se_jkn = jk$se,
                          se_lin = sl$se, df = sl$df, stringsAsFactors = FALSE),
         D_centered = D - psi, jkn_mode = jk$mode, df_jkn = jk$df)
  }

  sf <- one_arm(Q1s, Q0s, g_sf, "AIPW-SF")
  cf <- one_arm(Q1c, Q0c, g_cf, "AIPW-CF")

  drow <- data.frame(
    g_sf_min = min(g_sf), g_sf_max = max(g_sf),      # post-floor ranges
    g_cf_min = min(g_cf), g_cf_max = max(g_cf),
    g_sf_outside_floor = mean(g_sf_raw < g_floor | g_sf_raw > 1 - g_floor),
    g_cf_outside_floor = mean(g_cf_raw < g_floor | g_cf_raw > 1 - g_floor),
    cf_V_eff = attr(fold, "V_eff"),
    df_jkn = sf$df_jkn, jkn_mode = sf$jkn_mode,           # AIPW-SF arm
    df_jkn_cf = cf$df_jkn, jkn_mode_cf = cf$jkn_mode,     # AIPW-CF arm (a divergent
    # subset-vs-fallback choice between the arms must be visible, not masked)
    g_floor = g_floor, stringsAsFactors = FALSE)

  list(results = rbind(sf$row, cf$row),
       diagnostics = list(drow = drow, timing = c(sf = .t_sf, cf = .t_cf),
                          D_sf = sf$D_centered, D_cf = cf$D_centered))
}

# =====================================================================
# one_rep_aipw(i, pop, learners, g_floor, V_cf)  --  one simulation replication
#   i is the GLOBAL rep index: draw_sample(pop, sample_seed = SAMPLE_SEED_BASE+i)
#   makes the sample BYTE-IDENTICAL to the locked headline run (run_sim.R uses
#   the same formula), so the AIPW rows are paired rep-for-rep with the five
#   locked arms on the same scenario. Needs config.R/dgp.R/diagnostics.R
#   (SAMPLE_SEED_BASE, draw_sample, deff_clust) sourced -- run.R does this.
# =====================================================================
one_rep_aipw <- function(i, pop, learners, g_floor = 0.05, V_cf = 5L) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  r   <- aipw_arms(obs, learners = learners, g_floor = g_floor, V_cf = V_cf)
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(r$diagnostics$D_sf, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, r$results, deff = dd$deff_clust, icc_D = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design,
                    sumw_over_N = ch$sumw_over_N, w_cv = ch$w_cv),
    diag    = cbind(rep = i, r$diagnostics$drow,
                    deff = dd$deff_clust, icc_D = dd$icc_eif)
  )
}
