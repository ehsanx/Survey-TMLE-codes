# =====================================================================
# helpers_R11.R  —  NEW logic for run R11_resampling_eff (Web Appendix D)
#
# This file holds ALL run-specific logic so the canonical engine
# (R/*.R) is sourced READ-ONLY and never edited. It reuses the
# exposed building blocks from estimators.R:
#   run_estimators(), .se_des(), .eif_from_tmle(), make_cf_folds(), .sl()
# and from dgp.R / diagnostics.R: draw_sample(), deff_clust().
#
# Two additions:
#   (a) jk_se_on_eif()  — replication-based (survey JKn jackknife) SE of the
#       SAME cross-fitted / single-fit influence function whose Taylor
#       (linearization) SE the engine already returns. Corroborates Eq 8:
#       jackknife SE should ~ linearization SE (ratio ~ 1).
#   (b) ipw_svyglm_ate() — a survey-weighted IPW / g-computation (AIPW)
#       baseline built ONLY from svyglm nuisances, scored with the same
#       design-EIF machinery (.se_des), so its coverage and CI width are
#       directly comparable to the cross-fitted TMLE arms (efficiency
#       benchmark). Expectation: TMLE >= IPW efficiency, both valid coverage.
# =====================================================================

# ---------------------------------------------------------------------
# (a) JACKKNIFE (JKn) SE of an already-computed influence function -----
# ---------------------------------------------------------------------
# Given a per-unit EIF (centered: psi already subtracted by .eif_from_tmle)
# and the design (strata, cluster=PSU, weight), build the linearization
# design, convert it to a stratified delete-one-PSU jackknife replicate
# design (survey::as.svrepdesign, type='JKn'), and take the replicate SE of
# the weighted (Hajek) mean of the EIF. Because the EIF is centered, this is
# numerically the jackknife SE of the estimator itself (the EIF is the
# first-order delta-method linearization of psi-hat), so it is the correct
# apples-to-apples partner for the engine's Taylor SE.
#
# Returns list(se_jk, df_jk) or NA on failure (kept non-fatal: a single
# rep should never abort the chunk).
jk_se_on_eif <- function(eif, strata, cluster, weight) {
  out <- tryCatch({
    dat <- data.frame(eif = eif, strata = strata, cluster = cluster, weight = weight)
    # linearization design (same object .se_des uses internally)
    des_lin <- survey::svydesign(ids = ~cluster, strata = ~strata,
                                 weights = ~weight, data = dat)
    # stratified delete-one-PSU jackknife replicate design
    rep_des <- survey::as.svrepdesign(des_lin, type = "JKn")
    # replicate SE of the weighted mean of the (centered) EIF.
    # withReplicates recomputes weighted.mean(eif, w) on each jackknife
    # replicate-weight column and returns the jackknife variance.
    est <- survey::withReplicates(
      rep_des, quote(sum(.weights * eif) / sum(.weights)))
    list(se_jk = as.numeric(survey::SE(est)),
         df_jk = survey::degf(rep_des))
  }, error = function(e) list(se_jk = NA_real_, df_jk = NA_real_))
  out
}

# ---------------------------------------------------------------------
# (b) IPW / g-computation (AIPW) baseline from svyglm nuisances --------
# ---------------------------------------------------------------------
# Survey-weighted, parametric-only competitor to TMLE. We:
#   1. fit a survey-weighted propensity g(W)=P(A=1|W) via svyglm
#      (quasibinomial) on the SAME design;
#   2. fit a survey-weighted outcome Q(A,W) via svyglm (quasibinomial);
#   3. form the AIPW (doubly-robust) point estimate
#         psi = mean_w[ Q1 - Q0 + A/g (Y-Q1) - (1-A)/(1-g) (Y-Q0) ];
#   4. build its per-unit influence function (the AIPW EIF) and score the
#      design SE with the engine's .se_des() — identical variance machinery
#      to the TMLE arms, so coverage / CI width are directly comparable.
#
# g is floored at g_bound on BOTH tails (same spirit as the engine's
# g_oof_bound) to keep the IPW weights finite. This is the standard,
# defensible survey-IPW baseline; it is NOT cross-fitted (parametric, so it
# does not need to be), matching how an applied analyst would run svyglm.
#
# Args mirror run_estimators: obs (one draw_sample() output), W_cols, and
# g_bound. Returns a data.frame(method,b,se,df) row matching the engine.
ipw_svyglm_ate <- function(obs, W_cols = c("L1", "L2", "L3", "L4"),
                           g_bound = 0.05) {
  d   <- obs
  n   <- nrow(d)
  w   <- d$weight
  W   <- d[, W_cols, drop = FALSE]
  A   <- d$A
  Y   <- d$Y

  # design used to fit the weighted nuisances (clustered, stratified)
  des <- survey::svydesign(ids = ~cluster, strata = ~strata,
                           weights = ~weight, data = d)

  # 1. survey-weighted propensity g(W) = P(A=1|W)
  f_g  <- stats::as.formula(paste("A ~", paste(W_cols, collapse = " + ")))
  gfit <- survey::svyglm(f_g, design = des, family = stats::quasibinomial())
  g1   <- as.numeric(stats::predict(gfit, newdata = W, type = "response"))
  g1   <- pmin(pmax(g1, g_bound), 1 - g_bound)

  # 2. survey-weighted outcome Q(A,W) = E[Y|A,W]
  f_q  <- stats::as.formula(paste("Y ~ A +", paste(W_cols, collapse = " + ")))
  qfit <- survey::svyglm(f_q, design = des, family = stats::quasibinomial())
  Q1   <- as.numeric(stats::predict(qfit, newdata = data.frame(A = 1, W),
                                    type = "response"))
  Q0   <- as.numeric(stats::predict(qfit, newdata = data.frame(A = 0, W),
                                    type = "response"))
  Q1   <- pmin(pmax(Q1, 1e-3), 1 - 1e-3)
  Q0   <- pmin(pmax(Q0, 1e-3), 1 - 1e-3)
  QA   <- ifelse(A == 1, Q1, Q0)

  # 3. AIPW (doubly robust) point estimate, survey-weighted Hajek mean
  Hc  <- A / g1 - (1 - A) / (1 - g1)          # clever covariate
  blip_aug <- (Q1 - Q0) + Hc * (Y - QA)
  psi <- weighted.mean(blip_aug, w)

  # 4. AIPW influence function (same structure as the TMLE EIF, with
  #    parametric, non-targeted Q/g) -> design SE via the engine's .se_des
  eif <- blip_aug - psi
  se  <- .se_des(eif, d$strata, d$cluster, w, clustered = TRUE)

  list(row = data.frame(method = "IPW-svyglm", b = psi,
                        se = se$se, df = se$df),
       eif = eif, g1 = g1)
}
