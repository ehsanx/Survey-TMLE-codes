# =====================================================================
# estimands_rr_or.R  —  risk ratio (RR), odds ratio (OR), and the
# continuous ratio-of-means by the DELTA METHOD on the per-arm
# survey-TMLE influence functions.
#
# Theory (Web-Appendix RR/OR corollary): the per-arm targeted EIFs D1, D0
# (exposed by .eif_from_tmle as eif1/eif0) give a JOINT design-CLT for
# (psi1_hat, psi0_hat) via Cramer-Wold applied to c1*D1 + c0*D0 (lem:S).
# Any smooth contrast g(psi1, psi0) is then asymptotically normal by the
# delta method with design variance grad' Sigma grad, where Sigma is the
# design covariance of (psi1_hat, psi0_hat) = vcov(svymean(~D1 + D0))
# (Eq 8 / thm:var; the seminorm argument is generic so the plug-in
# variance carries to the transformed EIF). Gradients are Lipschitz
# because the Q-truncation in as:bounded keeps psi_a in [delta', 1-delta'].
#
# NON-DESTRUCTIVE post-processor: reads the per-arm components that
# run_estimators() now surfaces in diagnostics$arms. Used by the RR/OR
# demonstration simulation, the NHANES continuous example, and the
# svytmle Tier-1 feature.
# =====================================================================

suppressMessages({ library(survey) })

# Sigma = design covariance of (psi1_hat, psi0_hat) from the per-arm EIFs.
# Mirrors estimators.R::.se_des exactly: the simulation/whole-sample path
# (inpop = NULL) builds the design on the supplied rows; the domain path
# embeds the EIFs into full-length vectors and subset()s the FULL design
# (EpiMethods surveydata8), so partially-covered PSUs contribute correctly.
.cov2_des <- function(eif1, eif0, strata, cluster, weight,
                      clustered = TRUE, nest = FALSE, inpop = NULL) {
  if (is.null(inpop)) {
    dat <- data.frame(e1 = eif1, e0 = eif0, strata = strata,
                      cluster = cluster, weight = weight)
    des <- if (clustered)
      svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dat)
    else
      svydesign(ids = ~1, weights = ~weight, data = dat)
  } else {
    e1f <- numeric(length(inpop)); e1f[inpop] <- eif1
    e0f <- numeric(length(inpop)); e0f[inpop] <- eif0
    dat <- data.frame(e1 = e1f, e0 = e0f, strata = strata,
                      cluster = cluster, weight = weight, .inpop = inpop)
    des <- if (clustered)
      svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dat, nest = nest)
    else
      svydesign(ids = ~1, weights = ~weight, data = dat)
    des <- subset(des, .inpop)
  }
  m <- svymean(~e1 + e0, des)
  list(Sigma = vcov(m), df = survey::degf(des))
}

# Delta-method RR / OR (binary outcome) or ratio-of-means RR (continuous).
#   binary = TRUE  -> RR + OR on a [0,1] mean (probability scale)
#   binary = FALSE -> ratio of (positive) means only (continuous outcome; OR n/a)
# Returns one tidy row per estimand: point estimate, log-scale SE, t-df CI.
delta_rr_or <- function(psi1, psi0, eif1, eif0, strata, cluster, weight,
                        clustered = TRUE, level = 0.95, binary = TRUE,
                        eps = 1e-6, nest = FALSE, inpop = NULL) {
  cc   <- .cov2_des(eif1, eif0, strata, cluster, weight,
                    clustered = clustered, nest = nest, inpop = inpop)
  Sig  <- cc$Sigma; df <- cc$df
  crit <- qt(1 - (1 - level) / 2, max(1, df))
  mk <- function(estimand, log_est, grad) {
    v  <- as.numeric(crossprod(grad, Sig %*% grad))   # grad' Sigma grad
    se <- sqrt(max(v, 0))
    data.frame(estimand = estimand, est = exp(log_est), log_est = log_est,
               se_log = se, lo = exp(log_est - crit * se),
               hi = exp(log_est + crit * se), df = df,
               psi1 = psi1, psi0 = psi0, row.names = NULL)
  }
  if (binary) {
    p1 <- min(max(psi1, eps), 1 - eps); p0 <- min(max(psi0, eps), 1 - eps)
    rbind(
      mk("RR", log(p1) - log(p0),       c(1 / p1, -1 / p0)),
      mk("OR", qlogis(p1) - qlogis(p0), c(1 / (p1 * (1 - p1)), -1 / (p0 * (1 - p0))))
    )
  } else {
    if (psi1 <= 0 || psi0 <= 0)
      stop("ratio-of-means RR needs positive per-arm means; got psi1=", psi1, " psi0=", psi0)
    mk("RR", log(psi1) - log(psi0), c(1 / psi1, -1 / psi0))
  }
}
