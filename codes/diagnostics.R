# =====================================================================
# diagnostics.R  —  ICC / DEFF auditing for the survey-TMLE DGP
# The key gate (plan-phase2.md RC-2): the design effect that drives coverage
# lives on the EIF, not on Y. Measure it directly and clustering-only.
# =====================================================================

suppressMessages(library(survey))

# ---- population-level audit (realized ICCs on the fixed finite population) ---
audit_population <- function(pop_obj, n_sub = 40000L) {
  pop <- pop_obj$pop
  # subsample for the (slowish) ICC moment calc on the huge population
  idx <- if (nrow(pop) > n_sub) sort(sample.int(nrow(pop), n_sub)) else seq_len(nrow(pop))
  s <- pop[idx, ]
  list(
    icc_Y  = icc_anova(s$Y,  s$cluster),
    icc_A  = icc_anova(s$A,  s$cluster),
    icc_C1 = icc_anova(s$C1, s$cluster),
    icc_C3 = icc_anova(s$C3, s$cluster),
    mean_A = mean(pop$A),
    mean_Y = mean(pop$Y)
  )
}

# ---- clustering-only DEFF on a supplied EIF (the RC-2 gate) ------------------
# Both numerator and denominator are design SEs of the SAME weighted Hajek mean
# of the SAME eif; only clustering+strata are toggled -> divides out 1+CV(w)^2.
deff_clust <- function(eif, strata, cluster, weight) {
  dat <- data.frame(eif = eif, strata = strata, cluster = cluster, weight = weight)
  d_full <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dat)
  d_iid  <- svydesign(ids = ~1,                         weights = ~weight, data = dat)
  se_full <- as.numeric(SE(svymean(~eif, d_full)))
  se_iid  <- as.numeric(SE(svymean(~eif, d_iid)))
  list(deff_clust = (se_full / se_iid)^2,
       se_full = se_full, se_iid = se_iid,
       icc_eif = icc_anova(eif, cluster))
}

# ---- minimal GLM-TMLE used ONLY for early DGP de-risking (not a sim arm) -----
# Returns the targeted EIF + point estimate. weighted=TRUE -> Fully-Aware-style.
glm_tmle_eif <- function(obs, weighted = TRUE, gbound = 0.01) {
  W   <- obs[, c("L1", "L2", "L3", "L4")]
  wts <- if (weighted) obs$weight else rep(1, nrow(obs))
  qfit <- suppressWarnings(glm(Y ~ A + L1 + L2 + L3 + L4, data = obs,
                               family = binomial(), weights = wts))
  QA <- predict(qfit, type = "response")
  Q1 <- predict(qfit, newdata = data.frame(A = 1, W), type = "response")
  Q0 <- predict(qfit, newdata = data.frame(A = 0, W), type = "response")
  QA <- pmin(pmax(QA, 1e-4), 1 - 1e-4); Q1 <- pmin(pmax(Q1, 1e-4), 1 - 1e-4)
  Q0 <- pmin(pmax(Q0, 1e-4), 1 - 1e-4)
  gfit <- suppressWarnings(glm(A ~ L1 + L2 + L3 + L4, data = obs,
                               family = binomial(), weights = wts))
  g1 <- pmin(pmax(predict(gfit, type = "response"), gbound), 1 - gbound)
  H  <- obs$A / g1 - (1 - obs$A) / (1 - g1)
  fluc <- suppressWarnings(glm(obs$Y ~ -1 + H, family = binomial(),
                               offset = qlogis(QA), weights = wts))
  eps <- coef(fluc)
  Q1s <- plogis(qlogis(Q1) + eps * (1 / g1))
  Q0s <- plogis(qlogis(Q0) + eps * (-1 / (1 - g1)))
  QAs <- ifelse(obs$A == 1, Q1s, Q0s)
  psi <- weighted.mean(Q1s, wts) - weighted.mean(Q0s, wts)
  eif <- H * (obs$Y - QAs) + (Q1s - Q0s) - psi
  list(psi = psi, eif = eif, g_range = range(g1))
}
