# =====================================================================
# dgp_infsamp.R  —  R07_informative helper: INFORMATIVE-SAMPLING-BEYOND-C DGP
#
# Lives in THIS run folder. Reuses (does NOT edit) R/dgp.R. It EXTENDS the
# canonical two-stage sampling mechanism in draw_sample() by replacing the
# stage-1 SRS-of-PSUs with a PROBABILITY-PROPORTIONAL-TO-uS without-replacement
# selection, where uS is a per-PSU latent SELECTION driver that is correlated
# with the OUTCOME / TREATMENT PSU random effects (uY, uA) but is NOT a function
# of the confounders C. This deliberately breaks the sampling-ignorability
# assumption  S _||_ (A,Y) | C  that justifies the de-weighted out-of-fold (OOF)
# nuisance fits in the Fully-Aware-CF arm, so we can MEASURE the cost of
# de-weighting as a function of the informativeness knob rho.
#
# ---------------------------------------------------------------------
# SELECTION MECHANISM (math; see NOTES.md for the full write-up)
# ---------------------------------------------------------------------
# The canonical population (R/dgp.R make_population) already carries, per
# unit i, the PSU-level random intercepts:
#     uY_i  (outcome logit PSU effect; integrated OUT of the estimand)
#     uA_i  (treatment logit PSU effect)
# These are CONSTANT within a PSU. We aggregate them to one value per PSU.
#
# For each PSU g (globally-unique `cluster` id) we build a latent selection
# driver, standardized to mean 0 / unit variance across PSUs:
#     r_g    = z( uY_g + uA_g )                  # the "outcome/treatment" part
#     eps_g  ~ N(0,1)  i.i.d. across PSUs        # pure idiosyncratic noise
#     uS_g   = rho * r_g  +  sqrt(1 - rho^2) * eps_g
# so Var(uS_g) = 1 for every rho, and
#     Corr(uS_g, r_g) = rho.
#   * rho = 0 : uS is pure independent noise. Selection is informative on a
#               variable independent of (A,Y,C) -> sampling is STILL ignorable
#               given C (indeed unconditionally). BOTH arms should be unbiased:
#               this is the canonical "non-informative-given-C" reduction and
#               the negative-control end of the sweep.
#   * rho > 0 : uS correlates with uY,uA, which drive Y and A but are NOT
#               captured by C (C shares only uC, drawn independently of uY/uA in
#               R/dgp.R). Hence  S _||_ (A,Y) | C  is VIOLATED. The
#               de-weighted OOF fits (Fully-Aware-CF, unweighted) lose
#               consistency; the weighted-OOF variant should stay closer to
#               truth. As rho->1, uS -> r_g exactly.
#
# Stage-1 PSU selection (per stratum h): size measure
#     s_g = exp(lambda * uS_g)
# Select m_h PSUs WITHOUT replacement by SEQUENTIAL (systematic) PPS with target
# first-order inclusion probabilities
#     pi1_g = min(1, m_h * s_g / sum_{g' in h} s_{g'}),
# rescaled so sum_g pi1_g = m_h exactly (standard inclusion-probability
# normalization; capped values are held at 1 and the remainder is shared among
# the uncapped). Stage-2 is the canonical SRS-WOR of units within a chosen PSU
# (unchanged), giving pi2 = n0 / M. The Horvitz-Thompson weight is
#     w_i = 1 / (pi1_g * pi2),    so  E[sum_i in sample w_i] = N_pop,
# i.e. weights are design-unbiased for population totals (the Hajek mean the
# estimators target is consistent for Psi by construction). lambda controls the
# STRENGTH of the informative selection (weight spread); rho controls whether
# that selection is on something the C-only OOF fit can see.
#
# rho = 0 with lambda > 0 still produces UNEQUAL WEIGHTS (selection genuinely
# bites: w_cv > 0), but on a C-independent nuisance variable -> the de-weighting
# is harmless. This is exactly the design that lets the decision-rule check
# below be meaningful: the WEIGHTED estimator's bias must MOVE with rho (the
# selection truly perturbs the Y/A law in-sample), otherwise the curve is flat
# because the DGP is mis-implemented, not because de-weighting is safe.
#
# ---------------------------------------------------------------------
# This file defines:
#   attach_uS(pop, rho, seed)          -> per-PSU uS table (deterministic in seed)
#   draw_sample_infsamp(pop, sample_seed, rho, lambda, ...) -> obs (same contract
#       as draw_sample: L1..L4,A,Y,strata,cluster,weight,...; attr 'checks')
# No estimator/aggregation code is duplicated; run.R reuses run_estimators() and
# the .se_des/.sl building blocks from R/estimators.R.
# =====================================================================

# null-coalescing (R/dgp.R defines this too, but keep self-contained)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- standardize to mean 0 / unit variance (population PSUs) -----------------
.zscore <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) return(x - mean(x))
  (x - mean(x)) / s
}

# ---- build the per-PSU latent selection driver uS ---------------------------
# Deterministic given (population, rho, seed). Returns a data.frame keyed by the
# globally-unique PSU id `cluster`, with columns: strata, cluster, uS.
attach_uS <- function(population, rho, seed = 20260606L) {
  stopifnot(inherits(population, "svytmle_pop"))
  stopifnot(rho >= 0, rho <= 1)
  pop <- population$pop
  # one row per PSU; uY/uA/uC/b_h are constant within PSU so `unique` collapses
  psu <- unique(pop[, c("strata", "cluster", "uY", "uA", "uC", "b_h")])
  psu <- psu[order(psu$strata, psu$cluster), ]
  # "outcome/treatment" part of the driver, standardized across PSUs
  r <- .zscore(psu$uY + psu$uA)
  # idiosyncratic independent part (own RNG stream; does not touch SL CV RNG)
  eps_raw <- local({
    has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_old) old <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit({
      if (has_old) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    })
    stats::rnorm(nrow(psu))
  })
  # ORTHOGONALIZE the noise against the FULL span of ALL PSU-level drivers of the
  # data law -- {uY, uA, uC, b_h} plus an intercept -- NOT merely against r=z(uY+uA).
  # This is the crucial negative-control fix. Two leakage channels exist if we
  # project out too little:
  #   * projecting out only r leaves a residual eps-uY correlation (eps may load on
  #     uY as long as it loads oppositely on uA) -> at rho=0 the de-weighting biases;
  #   * leaving eps correlated with uC / b_h makes selection informative on the
  #     CONFOUNDERS at rho=0, so the raw in-sample E[Y] drifts (a "within-C" channel
  #     the estimator should ADJUST away, but which muddies a clean control).
  # Regressing eps on [1, uY, uA, uC, b_h] and keeping the residual makes eps _||_
  # every PSU effect individually -> at rho=0, uS = z(eps) is independent of
  # (A, Y, C): cor(uS, uY)=cor(uS, uA)=cor(uS, uC)=0, the unweighted PSU-mean of Y
  # does NOT drift, and BOTH arms are unbiased -- the proper negative control. As
  # rho grows, uS rotates toward r = z(uY+uA), reintroducing the BEYOND-C
  # informativeness (uY/uA are independent of uC, so C cannot see it).
  Bmat <- cbind(1, psu$uY, psu$uA, psu$uC, psu$b_h)
  # use .lm.fit() for the OLS residual: it is the projection onto the ORTHOGONAL
  # complement of span(Bmat). (Do NOT use qr.solve(Bmat, eps_raw) for an
  # overdetermined Bmat -- it returns a pivoted basic solution, not the
  # least-squares fit, leaving the residual NON-orthogonal to the basis.)
  eps_res <- as.numeric(.lm.fit(Bmat, eps_raw)$residuals)           # residual _|_ span(Bmat)
  eps <- .zscore(eps_res)                                            # unit variance
  uS <- rho * r + sqrt(max(0, 1 - rho^2)) * eps                      # cor(uS, r) == rho
  data.frame(strata = psu$strata, cluster = psu$cluster, uS = uS,
             r_psu = r, uY = psu$uY, uA = psu$uA, uC = psu$uC,
             stringsAsFactors = FALSE)
}

# ---- inclusion probabilities for SEQUENTIAL PPS-WOR within one stratum -------
# size s (length J), sample size m -> pi1 (length J), sum(pi1) == m, each in (0,1].
# Standard inclusion-probability construction with capping at 1 (Tille 2006):
# units whose proportional pi would exceed 1 are taken with certainty (pi=1) and
# the remaining sample size is reallocated proportionally among the rest.
.pps_incprob <- function(s, m) {
  J <- length(s)
  if (m >= J) return(rep(1, J))            # take all
  s <- pmax(s, .Machine$double.eps)
  pi <- rep(0, J)
  remaining <- seq_len(J)
  m_rem <- m
  repeat {
    p <- m_rem * s[remaining] / sum(s[remaining])
    over <- remaining[p >= 1]
    if (!length(over)) { pi[remaining] <- p; break }
    pi[over] <- 1
    m_rem <- m_rem - length(over)
    remaining <- setdiff(remaining, over)
    if (m_rem <= 0 || !length(remaining)) break
  }
  pi
}

# ---- systematic (sequential) PPS-WOR draw given inclusion probs pi ----------
# Returns the indices selected. Uses one uniform start + cumulative pi crossing
# integer boundaries (Madow systematic selection); exact first-order incl. probs.
.pps_systematic <- function(pi) {
  J <- length(pi)
  ord <- sample.int(J)                      # random permutation removes order bias
  cs  <- cumsum(pi[ord])
  m   <- round(sum(pi))
  if (m <= 0) return(integer(0))
  u   <- stats::runif(1)
  hits <- u + seq.int(0L, m - 1L)           # m equally spaced points in (0, m)
  sel_ord <- findInterval(hits, c(0, cs)) # which cumulative interval each lands in
  sel_ord <- pmin(pmax(sel_ord, 1L), J)
  unique(ord[sel_ord])
}

# ---- the informative-sampling sampler ---------------------------------------
# Same output contract as R/dgp.R draw_sample(): a data.frame `obs` with
# L1..L4, A, Y, strata, cluster, weight, psu_within, pi_i, .Q1_int, .Q0_int and
# attr(obs,'checks'). Differences vs canonical:
#   * stage-1 PSU selection is PPS-WOR on s_g = exp(lambda * uS_g)  (informative)
#   * stage-2 units-per-PSU n0 is BALANCED across strata (base_n0), so the ONLY
#     informativeness is the new PSU-level uS channel (clean rho sweep). We do
#     NOT also oversample high-b_h strata here; that canonical channel is C-based
#     and would muddy the "beyond C" interpretation.
draw_sample_infsamp <- function(population, sample_seed = NULL, model_type = NULL,
                                rho = 0, lambda = 0.8, uS_tab = NULL,
                                base_m = NULL, base_n0 = NULL,
                                uS_seed = 20260606L) {
  stopifnot(inherits(population, "svytmle_pop"))
  p   <- population$params
  pop <- population$pop
  mt  <- model_type %||% p$model_type
  scenario <- p$scenario
  J <- p$J_per_stratum; M <- p$M_per_psu; H <- p$H

  # canonical per-scenario stage sizes (mirror draw_sample defaults)
  if (scenario == "standard") { base_m  <- base_m  %||% 6L;  base_n0 <- base_n0 %||% 25L }
  if (scenario == "R1")       { base_m  <- base_m  %||% 2L;  base_n0 <- base_n0 %||% 20L }

  # per-PSU latent selection driver (deterministic; cache via uS_tab if provided)
  if (is.null(uS_tab)) uS_tab <- attach_uS(population, rho = rho, seed = uS_seed)

  m_h  <- rep(base_m,  H); names(m_h)  <- sort(unique(pop$strata))
  n0_h <- rep(base_n0, H); names(n0_h) <- names(m_h)

  # ---- RNG-scoped sampling (mirrors R/dgp.R .scoped_seed behaviour) ----
  samp <- local({
    has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_old) old <- get(".Random.seed", envir = .GlobalEnv)
    if (!is.null(sample_seed)) set.seed(sample_seed)
    on.exit({
      if (has_old) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    })

    keep <- integer(0); wt <- numeric(0); pii <- numeric(0)
    psu_tab <- unique(pop[, c("strata", "cluster")])
    for (hh in sort(unique(pop$strata))) {
      hk    <- as.character(hh)
      psus  <- psu_tab$cluster[psu_tab$strata == hh]
      # size measure for THIS stratum's PSUs, in the same order as `psus`
      uS_h  <- uS_tab$uS[match(psus, uS_tab$cluster)]
      s_h   <- exp(lambda * uS_h)
      pi1_h <- .pps_incprob(s_h, m_h[[hk]])             # first-order incl. probs
      sel_local <- .pps_systematic(pi1_h)                # STAGE 1: PPS-WOR of PSUs
      # guard: systematic selection can return != m due to rounding; top up/trim
      if (length(sel_local) > m_h[[hk]])
        sel_local <- sel_local[seq_len(m_h[[hk]])]
      if (length(sel_local) < m_h[[hk]]) {
        extra <- setdiff(order(pi1_h, decreasing = TRUE), sel_local)
        sel_local <- c(sel_local, extra[seq_len(m_h[[hk]] - length(sel_local))])
      }
      for (li in sel_local) {
        pg   <- psus[li]
        pi1  <- pi1_h[li]
        idx  <- which(pop$cluster == pg)
        take <- idx[sample.int(length(idx), min(n0_h[[hk]], length(idx)))]  # STAGE 2 SRS-WOR
        pi2  <- length(take) / M
        keep <- c(keep, take)
        pii  <- c(pii, rep(pi1 * pi2, length(take)))
        wt   <- c(wt,  rep(1 / (pi1 * pi2), length(take)))
      }
    }
    list(keep = keep, wt = wt, pii = pii)
  })

  s  <- pop[samp$keep, , drop = FALSE]
  Ls <- apply_L(s$C1, s$C2, s$C3, s$C4, mt)
  obs <- data.frame(
    Ls, A = s$A, Y = s$Y,
    strata = s$strata, cluster = s$cluster, weight = samp$wt,
    psu_within = s$psu_within, pi_i = samp$pii,
    .Q1_int = s$Q1_int, .Q0_int = s$Q0_int
  )

  # verification identities (cheap; same set as canonical draw_sample)
  ppst <- tapply(obs$cluster, obs$strata, function(x) length(unique(x)))
  attr(obs, "checks") <- list(
    n            = nrow(obs),
    sumw_over_N  = sum(obs$weight) / p$N_pop,            # ~1 in expectation (PPS-HT)
    min_psu_str  = min(ppst),                            # should be >= 2
    df_design    = length(unique(obs$cluster)) - length(unique(obs$strata)),
    w_cv         = sd(obs$weight) / mean(obs$weight),    # weight spread (rises with lambda)
    n_distinct_w = length(unique(round(obs$weight, 6))),
    rho          = rho,
    lambda       = lambda
  )
  obs
}
