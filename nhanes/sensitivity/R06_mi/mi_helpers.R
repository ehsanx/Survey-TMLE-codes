# =====================================================================
# Nhanes/arc_runs/R06_mi/mi_helpers.R
#   Helper logic for the R06_mi enhancement run. NON-DESTRUCTIVE: this file
#   lives in the run's own folder and only REUSES the canonical engine
#   (estimators.R via run_estimators) plus mice. It edits NOTHING in codes/
#   or Nhanes/R/.
#
#   Two pieces of new logic:
#     (1) impute_m()   — mirrors Nhanes/R/02_impute.R's mice setup but with
#                        m = M completed datasets (02_impute.R hardcodes m = 1).
#                        Returns a list of M obs frames, each with the imputed
#                        covariate values written back into the FULL MEC frame's
#                        inpop rows (identical write-back contract to 02_impute.R)
#                        so each is a drop-in for run_estimators(... inpop=...).
#     (2) rubin_pool() — combines the M per-imputation point estimates + SEs of
#                        ONE arm via Rubin's rules: between + within variance,
#                        total variance, t reference with Barnard-Rubin df.
# =====================================================================

suppressMessages({ library(mice) })

# ---- (1) m-imputation, reusing 02_impute.R's exact mice configuration --------
# obs : the *_analytic.rds frame (has attr 'covs','example', column inpop, design
#       columns SDMVSTRA/SDMVPSU/WTMEC_POOLED, A, Y).
# M   : number of completed datasets.
# seed: base mice seed (02_impute.R used 20260607L). Same defaults as 02_impute:
#       maxit = 5, A and Y used as PREDICTORS but never imputed (observed by
#       construction). If no covariate is missing, returns M identical copies
#       (mirrors 02_impute.R's `else dc <- d`), so downstream pooling still works.
impute_m <- function(obs, M = 10L, seed = 20260607L, maxit = 5L) {
  covs <- attr(obs, "covs")
  sub  <- which(obs$inpop)
  d    <- obs[sub, c(covs, "A", "Y")]
  # character -> factor for mice (same as 02_impute.R)
  for (cc in names(d)) if (is.character(d[[cc]])) d[[cc]] <- factor(d[[cc]])
  miss_before <- colSums(is.na(d[, covs, drop = FALSE]))
  to_impute   <- names(miss_before)[miss_before > 0]

  if (length(to_impute)) {
    set.seed(seed)
    mi <- mice(d, m = M, maxit = maxit, printFlag = FALSE, seed = seed)
    completed <- lapply(seq_len(M), function(k) complete(mi, k))
  } else {
    completed <- replicate(M, d, simplify = FALSE)   # nothing to impute
  }

  # write each completed set's covariates back into a fresh copy of the FULL frame
  out <- lapply(completed, function(dc) {
    o <- obs
    for (cc in covs) o[[cc]][sub] <- dc[[cc]]
    attr(o, "covs")    <- covs
    attr(o, "example") <- attr(obs, "example")
    o
  })
  attr(out, "n_imputed_covs") <- length(to_impute)
  attr(out, "imputed_covs")   <- to_impute
  out
}

# ---- (2) Rubin's rules for a single arm --------------------------------------
# b  : length-M vector of per-imputation point estimates (Q_l).
# se : length-M vector of per-imputation standard errors    (sqrt(U_l)).
# df_complete : the complete-data design df for ONE imputation (all M share the
#       same design, so they are equal; we take df[1]). Used for the Barnard-
#       Rubin (1999) small-sample df correction so the pooled df never exceeds
#       the complete-data design df.
# Returns: pooled estimate qbar, within Ubar, between B, total T, pooled se,
#          fraction of missing information (fmi), relative increase in variance
#          (riv), and Barnard-Rubin df (df_br). Matches Rubin (1987) /
#          Barnard & Rubin (1999); this is exactly what mice::pool() computes.
rubin_pool <- function(b, se, df_complete) {
  M     <- length(b)
  qbar  <- mean(b)                      # pooled point estimate
  Ubar  <- mean(se^2)                   # within-imputation variance
  B     <- if (M > 1L) var(b) else 0    # between-imputation variance (n-1 denom)
  Tvar  <- Ubar + (1 + 1 / M) * B       # total variance
  se_p  <- sqrt(Tvar)

  if (M > 1L && B > 0) {
    # classic Rubin (1987) df
    riv     <- (1 + 1 / M) * B / Ubar           # relative increase in variance
    df_old  <- (M - 1) * (1 + 1 / riv)^2
    # Barnard-Rubin (1999) small-sample correction toward the complete-data df
    lambda  <- (1 + 1 / M) * B / Tvar           # fraction of missing information core
    df_obs  <- if (is.finite(df_complete) && df_complete > 0) {
                 (df_complete + 1) / (df_complete + 3) * df_complete * (1 - lambda)
               } else Inf
    df_br   <- if (is.finite(df_obs)) df_old * df_obs / (df_old + df_obs) else df_old
    fmi     <- (riv + 2 / (df_br + 3)) / (riv + 1)
  } else {
    riv <- 0; df_br <- df_complete; fmi <- 0; lambda <- 0
  }

  list(b = qbar, se = se_p, df = df_br,
       Ubar = Ubar, B = B, Tvar = Tvar, riv = riv, fmi = fmi, M = M)
}

# ---- one-hot encode a completed frame's domain covariates (mirrors 03_*.R) ---
# Returns the obs frame with numeric dummy columns added on inpop rows, plus the
# design columns set, plus the dummy column names (W_cols for run_estimators).
encode_domain <- function(obs) {
  covs <- attr(obs, "covs"); sub <- which(obs$inpop)
  covdf <- droplevels(obs[sub, covs, drop = FALSE])
  mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
  stopifnot(nrow(mm) == length(sub))
  wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
  for (j in seq_along(wnames)) {
    obs[[wnames[j]]] <- NA_real_
    obs[[wnames[j]]][sub] <- mm[, j]
  }
  obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED
  list(obs = obs, wnames = wnames, sub = sub)
}
