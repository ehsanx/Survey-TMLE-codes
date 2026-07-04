# =====================================================================
# nhanes/arc_runs/R09_sensitivity/sensitivity_helpers.R
#   Helper logic for the R09_sensitivity (causal-validity) enhancement run.
#
#   NON-DESTRUCTIVE: this file lives in the run's own folder and only REUSES the
#   canonical engine (R/estimators.R via run_estimators) plus the already
#   IMPUTED analytic frames (nhanes/analytic/<id>_imputed.rds). It edits NOTHING
#   in codes/ or nhanes/.
#
#   Pieces of new logic (all are covariate-set swaps or design tabulations; no
#   new estimator math):
#     (1) set_covs()        — return a copy of an imputed obs frame with attr
#                             'covs' overridden to a modified adjustment set
#                             (e.g. BMI omitted, PHQ-9 omitted). The dropped
#                             covariate's COLUMN is left in the frame untouched;
#                             only the 'covs' attribute (which 03/encode read to
#                             build the design matrix) changes. This is exactly
#                             the "alter the covariate set via attr 'covs'" hook
#                             that 01_build_analytic.R sets and 03 consumes.
#     (2) encode_domain()   — one-hot encode the domain covariates named in
#                             attr 'covs' to a numeric design matrix and set the
#                             survey design columns. BYTE-IDENTICAL to the
#                             encoder in 03_run_estimators.R / nhanes_arc.R /
#                             R06_mi (copied here so this run is self-contained).
#     (3) run_primary()     — run the PRIMARY (3-learner near-Donsker
#                             library, == 03_run_estimators.R LIB) on a given
#                             adjustment set and return the Fully-Aware-CF row
#                             plus all arms, with t-CI. Wraps run_estimators with
#                             inpop = domain, nest = TRUE (NHANES SDMVPSU repeats
#                             across strata).
#     (4) ci_delta()        — ATE/CI delta of a sensitivity vs its main analysis
#                             (same example, full adjustment set): difference in
#                             point estimate and in each CI bound, and whether the
#                             qualitative conclusion (CI excludes 0) flipped.
#     (5) strata_cycle_check() — the empirical SDMVSTRA cross-cycle NON-OVERLAP
#                             check: tabulate strata x cycle and confirm no
#                             SDMVSTRA value spans >1 NHANES cycle in the pooled
#                             design (which is what makes nest=TRUE / the pooled
#                             stratification valid). A quick tabulation, NOT an
#                             estimator run.
# =====================================================================

suppressMessages({ library(survey) })

# ---- (1) override the adjustment set via attr 'covs' -------------------------
# obs  : an imputed *_imputed.rds frame (has attr 'covs','example').
# drop : character vector of covariate names to REMOVE from the adjustment set.
# Returns a copy of obs with attr 'covs' = setdiff(old covs, drop). Errors if a
# name in `drop` is not actually in the current covs (guards against typos /
# stale specs). The underlying data columns are left in place; only the
# adjustment set (the attribute the encoder reads) shrinks.
set_covs <- function(obs, drop = character(0)) {
  covs <- attr(obs, "covs")
  miss <- setdiff(drop, covs)
  if (length(miss))
    stop("set_covs: covariate(s) not in adjustment set: ", paste(miss, collapse = ", "),
         " | available: ", paste(covs, collapse = ", "))
  new_covs <- setdiff(covs, drop)
  attr(obs, "covs") <- new_covs
  attr(obs, "dropped_covs") <- drop
  obs
}

# ---- (2) one-hot encode the domain covariates (mirrors 03_run_estimators.R) --
# Identical contract to nhanes_arc.R / 03_run_estimators.R / R06_mi encode_domain:
# build mm on the inpop rows, write numeric dummy columns back onto inpop rows,
# set strata = SDMVSTRA, cluster = SDMVPSU, weight = WTMEC_POOLED.
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

# ---- (3) run the primary on a given adjustment set -----------------
# obs      : an imputed frame whose attr 'covs' is the adjustment set to use
#            (already passed through set_covs() if dropping anything).
# LIB      : SuperLearner library (default the primary, == 03_*.R).
# seed     : run seed (the CF/SL split seed). 03_run_estimators.R used 20260607L.
# Returns a data.frame of all five arms with method,b,se,df,lcl,ucl, plus the
# adjustment-set label, n_domain, and the FA single-fit g-range (overlap quick
# look). nest = TRUE because NHANES SDMVPSU ids repeat across strata.
run_primary <- function(obs, LIB = c("SL.glm", "SL.earth", "SL.glmnet"),
                        seed = 20260607L) {
  enc <- encode_domain(obs)
  set.seed(seed)
  r <- run_estimators(enc$obs, learners = LIB, V_cf = 5L, inner_cv_folds = 5L,
                      W_cols = enc$wnames, nest = TRUE, inpop = enc$obs$inpop)
  rr <- r$results
  crit <- qt(0.975, pmax(1, rr$df))
  rr$lcl <- rr$b - crit * rr$se
  rr$ucl <- rr$b + crit * rr$se
  rr$signif <- (rr$lcl > 0) | (rr$ucl < 0)       # CI excludes 0
  attr(rr, "g_fa") <- range(r$diagnostics$g_fa)
  attr(rr, "g_cf") <- range(r$diagnostics$g_cf)
  attr(rr, "n_domain") <- length(enc$sub)
  attr(rr, "covs") <- attr(obs, "covs")
  rr
}

# ---- (4) ATE/CI delta of a sensitivity arm vs the main analysis --------------
# main, sens : the all-arms data.frames from run_primary() for, respectively, the
#              MAIN (full adjustment) and the SENSITIVITY (modified) fits of the
#              SAME example. We report the delta on a chosen arm (default the
#              primary, Fully-Aware-CF).
# Returns a one-row data.frame: the main and sensitivity b/CI, the deltas, and a
# flag for whether the qualitative conclusion (CI excludes 0) changed.
ci_delta <- function(main, sens, arm = "Fully-Aware-CF") {
  m <- main[main$method == arm, ]
  s <- sens[sens$method == arm, ]
  data.frame(
    arm        = arm,
    b_main     = m$b,   lcl_main = m$lcl, ucl_main = m$ucl, signif_main = m$signif,
    b_sens     = s$b,   lcl_sens = s$lcl, ucl_sens = s$ucl, signif_sens = s$signif,
    d_b        = s$b   - m$b,
    d_lcl      = s$lcl - m$lcl,
    d_ucl      = s$ucl - m$ucl,
    rel_d_b    = if (m$b != 0) (s$b - m$b) / abs(m$b) else NA_real_,
    conclusion_flip = (m$signif != s$signif),
    stringsAsFactors = FALSE)
}

# ---- (5) empirical SDMVSTRA cross-cycle non-overlap check --------------------
# obs : an analytic/imputed frame with SDMVSTRA + SDDSRVYR (cycle code) columns.
# domain_only : if TRUE, restrict to inpop rows (the domain actually analyzed);
#               default FALSE -> the full MEC design (what svydesign sees).
# Returns a list: n_strata, n_cycles, max #cycles ANY stratum spans, the number
# of strata that span >1 cycle (MUST be 0 for the pooled stratification to be
# valid -> nest=TRUE safe), and the offending strata if any. A pure tabulation.
strata_cycle_check <- function(obs, domain_only = FALSE) {
  d <- if (domain_only) obs[obs$inpop, ] else obs
  tab <- table(strata = d$SDMVSTRA, cycle = d$SDDSRVYR)
  n_cyc_per_stratum <- rowSums(tab > 0)
  multi <- names(n_cyc_per_stratum)[n_cyc_per_stratum > 1]
  list(
    domain_only      = domain_only,
    n_strata         = nrow(tab),
    n_cycles         = ncol(tab),
    max_cycles_spanned = as.integer(max(n_cyc_per_stratum)),
    n_strata_multicycle = length(multi),
    overlap_clean    = length(multi) == 0L,      # TRUE = no stratum spans >1 cycle
    offending_strata = multi,
    strata_range     = range(as.integer(rownames(tab))),
    cycles_present   = colnames(tab))
}
