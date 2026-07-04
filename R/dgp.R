# =====================================================================
# dgp.R  —  canonical data-generating process for the survey-TMLE sims
# Replaces the 6-way-duplicated create.data() block. See plan-phase2.md.
#
# Two functions:
#   make_population(scenario, ...) -> fixed finite population + truth (run ONCE)
#   draw_sample(population, ...)    -> one multistage design sample (run per rep)
#
# Design corrections baked in (plan-phase2.md RC-1/2/3):
#   RC-1  outcome PSU random intercept u^Y is INDEPENDENT of u^C and of A;
#         the truth integrates u^Y OUT (Gauss-Hermite), so it equals the
#         estimator's target Q0(a,c)=E[Y|A=a,C=c].
#   RC-2  an INDEPENDENT treatment PSU random intercept u^A is the strongest
#         lever for the design effect ON THE EIF; estimand-safe.
#   RC-3  weights are unequal BY DESIGN via differential stratum sampling
#         fractions (oversample high-b_h strata); SRS at both stages keeps
#         inclusion probs exact and sum(w) = N exactly.
# =====================================================================

# null-coalescing helper (defined first; used throughout)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Gauss-Hermite nodes/weights in BASE R (Golub-Welsch) so we need no `statmod`
# dependency (it is not in the ARC R library). Physicists' convention:
# integral f(x) exp(-x^2) dx  ~=  sum_i weights[i] * f(nodes[i]).
.gauss_hermite <- function(n) {
  i <- seq_len(n - 1L)
  b <- sqrt(i / 2)
  J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- b; J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(nodes = e$values[o], weights = sqrt(pi) * (e$vectors[1, o])^2)
}

# ---- RNG scoping: keep DGP/sample RNG from leaking into SuperLearner CV -----
.scoped_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_old) old <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(seed)
  on.exit(if (has_old) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                 rm(".Random.seed", envir = .GlobalEnv))
  force(expr)
}

# ---- Gauss-Hermite: E_{u~N(0,s2)}[ expit(m + u) ] for a vector of linear preds
# integral expit(m+u) phi(u;0,s2) du  ~=  sum_q (w_q/sqrt(pi)) expit(m + sqrt(2)*s*z_q)
gh_expit <- function(m_lin, sigma2, nGH = 30L) {
  if (sigma2 <= 0) return(plogis(m_lin))
  gh    <- .gauss_hermite(nGH)
  nodes <- sqrt(2 * sigma2) * gh$nodes
  wts   <- gh$weights / sqrt(pi)
  out <- numeric(length(m_lin))
  for (q in seq_along(nodes)) out <- out + wts[q] * plogis(m_lin + nodes[q])
  out
}

# ---- one-way ANOVA / method-of-moments ICC (continuous OR 0/1 x) ------------
icc_anova <- function(x, cluster) {
  grp <- split(x, cluster)
  nj  <- vapply(grp, length, integer(1))
  M   <- length(grp); N <- length(x); gm <- mean(x)
  if (M < 2L) return(NA_real_)
  MSB <- sum(nj * (vapply(grp, mean, numeric(1)) - gm)^2) / (M - 1)
  MSW <- sum(vapply(grp, function(v) sum((v - mean(v))^2), numeric(1))) / max(1, (N - M))
  n0  <- (N - sum(nj^2) / N) / (M - 1)
  s2b <- (MSB - MSW) / n0
  max(0, s2b / (s2b + MSW))
}

# ---- Kang-Schafer-style covariate transforms (B5: CODE/canonical form) ------
# Decision (2026-06-06): keep the numerically-sane CODE form; the APPENDIX is
# the artifact to edit, NOT the code. (Appendix L3=(25*C1*C3+0.6)^3 explodes.)
apply_L <- function(C1, C2, C3, C4, model_type) {
  if (model_type == "complex") {
    data.frame(
      L1 = exp(C1 / 2),
      L2 = C2 / (1 + exp(C1)) + 10,
      L3 = (C1 * C3 / 25 + 0.6)^3,
      L4 = (C2 + C4 + 20)^2
    )
  } else {
    data.frame(L1 = C1, L2 = C2, L3 = C3, L4 = C4)
  }
}

# ---- scenario geometry (frozen table, plan-phase2.md sec 1) -----------------
.scenario_defaults <- function(scenario) {
  switch(scenario,
    standard = list(H = 10L,  J_per_stratum = 60L, M_per_psu = 200L),   # N=120k, 600 pop PSUs
    R1       = list(H = 50L,  J_per_stratum = 8L,  M_per_psu = 300L),   # N=120k, 400 pop PSUs
    stop("scenario must be 'standard' or 'R1'")
  )
}

# ============================================================================
# make_population(): build ONE fixed finite population + compute the truth ONCE
# ============================================================================
make_population <- function(scenario      = c("standard", "R1"),
                            model_type    = c("complex", "simple"),
                            outcome_type  = c("binary", "continuous"),
                            te_log_odds   = 1.5,
                            pop_seed      = 20260606L,
                            H = NULL, J_per_stratum = NULL, M_per_psu = NULL,
                            beta_strat    = 1.0,     # half-range of fixed stratum grid (informativeness)
                            sigma2_C      = 0.3,     # PSU effect in confounders  (ICC_C ~ 0.23)
                            sigma2_Y      = 1.5,     # PSU effect in OUTCOME logit (drives EIF clustering/DEFF)
                            sigma2_A      = 0.8,     # PSU effect in TREATMENT logit (clusters the clever cov.)
                            sigma2_eps    = 1.0,     # residual variance of the CONTINUOUS outcome (ignored if binary)
                            alpha_g       = log(1.3),# confounding strength (moderate per RC-2, keeps overlap)
                            p_treat_target= 0.35,    # marginal P(A=1) (relaxed from rare GDM; see plan)
                            gamma0        = -2,      # outcome intercept
                            gamma_C       = 0.5,     # confounder -> outcome
                            truth_M       = 2e6L,    # huge-N MC draws for the super-pop truth
                            nGH           = 30L) {
  scenario     <- match.arg(scenario)
  model_type   <- match.arg(model_type)
  outcome_type <- match.arg(outcome_type)
  d <- .scenario_defaults(scenario)
  H             <- H %||% d$H
  J_per_stratum <- J_per_stratum %||% d$J_per_stratum
  M_per_psu     <- M_per_psu %||% d$M_per_psu
  theta <- te_log_odds
  N_pop <- H * J_per_stratum * M_per_psu

  built <- .scoped_seed(pop_seed, {
    # ---- hierarchy indices ----
    n_psu      <- H * J_per_stratum
    strata     <- rep(seq_len(H), each = J_per_stratum * M_per_psu)
    psu_within <- rep(rep(seq_len(J_per_stratum), each = M_per_psu), times = H)
    psu_global <- (strata - 1L) * J_per_stratum + psu_within          # globally unique PSU id

    # ---- fixed stratum grid + INDEPENDENT PSU random effects (drawn ONCE) ----
    b_grid <- if (H > 1) beta_strat * ((seq_len(H) - (H + 1) / 2) / ((H - 1) / 2)) else 0
    uC_psu <- rnorm(n_psu, 0, sqrt(sigma2_C))
    uY_psu <- rnorm(n_psu, 0, sqrt(sigma2_Y))   # independent of uC
    uA_psu <- rnorm(n_psu, 0, sqrt(sigma2_A))   # independent of uY (RC-2)

    b_h  <- b_grid[strata]
    uC_i <- uC_psu[psu_global]
    uY_i <- uY_psu[psu_global]
    uA_i <- uA_psu[psu_global]

    # ---- confounders: 2 clustered (share b_h + uC), 2 idiosyncratic ----
    C1 <- b_h + uC_i + rnorm(N_pop)
    C2 <- b_h + uC_i + rnorm(N_pop)
    C3 <- rnorm(N_pop)
    C4 <- rnorm(N_pop)
    Csum <- C1 + C2 + C3 + C4

    # ---- treatment: confounding via C + INDEPENDENT cluster effect uA ----
    # solve alpha0 so the realized marginal P(A=1) hits the target
    eta_noint <- alpha_g * Csum + uA_i
    f_root <- function(a0) mean(plogis(a0 + eta_noint)) - p_treat_target
    alpha0 <- uniroot(f_root, lower = -12, upper = 12)$root
    pscore <- plogis(alpha0 + eta_noint)
    A <- rbinom(N_pop, 1, pscore)

    # ---- outcome: m(a,C) + INDEPENDENT cluster effect uY ----
    m1 <- gamma0 + theta * 1 + gamma_C * Csum
    m0 <- gamma0 + theta * 0 + gamma_C * Csum
    if (outcome_type == "binary") {
      Q1_real <- plogis(m1 + uY_i)               # realized-u potential-outcome means
      Q0_real <- plogis(m0 + uY_i)
      QA <- ifelse(A == 1, Q1_real, Q0_real)
      Y  <- rbinom(N_pop, 1, QA)
      # ---- u-INTEGRATED conditional means (what the TMLE targets): Q0(a,c) ----
      Q1_int <- gh_expit(m1, sigma2_Y, nGH)
      Q0_int <- gh_expit(m0, sigma2_Y, nGH)
    } else {
      # CONTINUOUS (Gaussian, identity link). E[uY]=0 and uY _||_ C, so the
      # u-integrated target is E[Y|A=a,C=c]=m_a (no Gauss-Hermite needed).
      Q1_real <- m1 + uY_i                       # realized-u conditional means
      Q0_real <- m0 + uY_i
      QA <- ifelse(A == 1, Q1_real, Q0_real)
      Y  <- QA + rnorm(N_pop, 0, sqrt(sigma2_eps))
      Q1_int <- m1
      Q0_int <- m0
    }

    list(
      pop = data.frame(strata, psu_within, cluster = psu_global,
                       C1, C2, C3, C4, pscore, A, Y,
                       uC = uC_i, uY = uY_i, uA = uA_i, b_h = b_h,
                       Q1_real, Q0_real, Q1_int, Q0_int),
      alpha0 = alpha0
    )
  })
  pop    <- built$pop
  alpha0 <- built$alpha0

  # ---- TRUTH (computed ONCE) ----
  # (a) super-population Psi: fresh huge-M MC over (C,u) law, u^Y integrated out
  truth <- .scoped_seed(pop_seed + 1L, {
    b_grid <- if (H > 1) beta_strat * ((seq_len(H) - (H + 1) / 2) / ((H - 1) / 2)) else 0
    h   <- sample.int(H, truth_M, replace = TRUE)
    bm  <- b_grid[h]
    uCm <- rnorm(truth_M, 0, sqrt(sigma2_C))
    Csum_m <- (bm + uCm) + rnorm(truth_M) +     # C1
              (bm + uCm) + rnorm(truth_M) +     # C2
              rnorm(truth_M) + rnorm(truth_M)   # C3, C4
    m1m <- gamma0 + theta + gamma_C * Csum_m
    m0m <- gamma0 +         gamma_C * Csum_m
    if (outcome_type == "binary") {
      mu1 <- gh_expit(m1m, sigma2_Y, nGH); mu0 <- gh_expit(m0m, sigma2_Y, nGH)
    } else {
      mu1 <- m1m; mu0 <- m0m                    # identity link, uY integrates out
    }
    blip <- mu1 - mu0
    list(psi = mean(blip), se_mc = sd(blip) / sqrt(truth_M),
         psi1 = mean(mu1), psi0 = mean(mu0))    # per-arm means -> RR/OR truth
  })
  # (b) finite-pop census parameters on THIS realized population
  psi_N_Q0   <- mean(pop$Q1_int  - pop$Q0_int)    # u-integrated  (what estimator targets here)
  psi_N_real <- mean(pop$Q1_real - pop$Q0_real)   # realized-u    (actual realized effect)
  psi1_N     <- mean(pop$Q1_int); psi0_N <- mean(pop$Q0_int)  # per-arm census means (RR/OR)

  structure(
    list(
      pop   = pop,
      truth = list(psi        = truth$psi,       # <- coverage target (super-pop ATE, RD)
                   se_mc      = truth$se_mc,
                   psi1       = truth$psi1,       # super-pop per-arm means (RR/OR coverage target)
                   psi0       = truth$psi0,
                   psi1_N     = psi1_N,           # finite-pop census per-arm means
                   psi0_N     = psi0_N,
                   psi_N_Q0   = psi_N_Q0,
                   psi_N_real = psi_N_real,
                   gap_super_census = truth$psi - psi_N_Q0),
      params = list(scenario = scenario, model_type = model_type, outcome_type = outcome_type,
                    te_log_odds = te_log_odds,
                    H = H, J_per_stratum = J_per_stratum, M_per_psu = M_per_psu, N_pop = N_pop,
                    beta_strat = beta_strat, sigma2_C = sigma2_C, sigma2_Y = sigma2_Y,
                    sigma2_A = sigma2_A, sigma2_eps = sigma2_eps, alpha_g = alpha_g,
                    p_treat_target = p_treat_target,
                    alpha0 = alpha0, gamma0 = gamma0, gamma_C = gamma_C,
                    pop_seed = pop_seed, nGH = nGH)
    ),
    class = "svytmle_pop"
  )
}

print.svytmle_pop <- function(x, ...) {
  p <- x$params; t <- x$truth; cat(sprintf(
    "svytmle_pop [%s/%s]  N=%d  H=%d  J/str=%d  M/PSU=%d\n  Psi(super-pop)=%.5f  (mc se %.6f)\n  Psi_N(u-int)=%.5f  Psi_N(realized-u)=%.5f  |super-census gap|=%.6f\n  alpha0=%.3f (target P(A=1)=%.2f)\n",
    p$scenario, p$model_type, p$N_pop, p$H, p$J_per_stratum, p$M_per_psu,
    t$psi, t$se_mc, t$psi_N_Q0, t$psi_N_real, abs(t$gap_super_census),
    p$alpha0, p$p_treat_target)); invisible(x)
}

# ============================================================================
# draw_sample(): one two-stage stratified cluster sample with unequal weights
# ============================================================================
draw_sample <- function(population, sample_seed = NULL, model_type = NULL,
                        base_m = NULL, base_n0 = NULL,
                        alpha_strat = 2.0, oversample = TRUE) {
  stopifnot(inherits(population, "svytmle_pop"))
  p   <- population$params
  pop <- population$pop
  mt  <- model_type %||% p$model_type
  scenario <- p$scenario
  J <- p$J_per_stratum; M <- p$M_per_psu; H <- p$H

  if (scenario == "standard") { base_m  <- base_m  %||% 6L;  base_n0 <- base_n0 %||% 25L }
  if (scenario == "R1")       { base_m  <- base_m  %||% 2L;  base_n0 <- base_n0 %||% 20L }

  # ---- informativeness via UNITS-per-PSU (n0_h), with m_h BALANCED ------------
  # n0_h increases in b_h so high-b_h units are oversampled -> unequal weights and
  # a sample C-distribution shifted from the population (the bias the weights fix).
  # Keeping m_h equal across strata avoids sparse (2-PSU) strata, which would
  # (a) thin the design df and (b) force the CF arm into a cross-stratum
  # extrapolating stratum-block fold split. (R1 is intrinsically m_h = 2.)
  strat_tab <- unique(pop[, c("strata", "b_h")])
  strat_tab <- strat_tab[order(strat_tab$strata), ]
  f_h <- if (oversample) exp(alpha_strat * strat_tab$b_h) else rep(1, H)
  f_h <- f_h / mean(f_h)
  m_h  <- rep(base_m, H)
  n0_h <- pmax(4L, round(base_n0 * f_h))
  names(m_h) <- names(n0_h) <- strat_tab$strata

  samp <- .scoped_seed(sample_seed, {
    keep <- integer(0); wt <- numeric(0); pii <- numeric(0)
    psu_tab <- unique(pop[, c("strata", "cluster")])
    for (hh in strat_tab$strata) {
      hk    <- as.character(hh)
      psus  <- psu_tab$cluster[psu_tab$strata == hh]
      sel   <- psus[sample.int(length(psus), m_h[[hk]])]      # STAGE 1: SRS-WOR of PSUs
      pi1   <- m_h[[hk]] / J
      for (pg in sel) {
        idx  <- which(pop$cluster == pg)
        take <- idx[sample.int(length(idx), min(n0_h[[hk]], length(idx)))]  # STAGE 2: SRS-WOR
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
    # carry latent truth pieces for diagnostics (NOT used by estimators):
    .Q1_int = s$Q1_int, .Q0_int = s$Q0_int
  )

  # ---- verification identities (cheap; assert design correctness) ----
  ppst <- tapply(obs$cluster, obs$strata, function(x) length(unique(x)))
  attr(obs, "checks") <- list(
    n            = nrow(obs),
    sumw_over_N  = sum(obs$weight) / p$N_pop,                 # should be ~1 (exact for SRS)
    min_psu_str  = min(ppst),                                  # should be >= 2
    df_design    = length(unique(obs$cluster)) - length(unique(obs$strata)),
    w_cv         = sd(obs$weight) / mean(obs$weight),          # weight variability (0 if equal)
    n_distinct_w = length(unique(round(obs$weight, 6)))
  )
  obs
}
