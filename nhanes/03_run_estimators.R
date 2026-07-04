# =====================================================================
# nhanes/03_run_estimators.R  —  the five survey-TMLE arms per example
#
# For each example: load the imputed frame, one-hot encode the (imputed) domain
# covariates to a numeric design matrix, set the survey design columns
# (strata=SDMVSTRA, cluster=SDMVPSU, weight=WTMEC_POOLED), and call the generalized
# run_estimators() with inpop = the domain. run_estimators builds the design on the
# FULL MEC sample and subset()s it to the domain for the design-based variance
# (proper sub-population handling). Five arms: Non-Aware, Partially-Aware,
# Fully-Aware (single-fit), Fully-Aware-CV, Fully-Aware-CF (primary).
#
# Saves nhanes/results/<id>/estimators_all_arms.rds and a combined CSV.
# Usage: Rscript nhanes/03_run_estimators.R            (all four)
#        Rscript nhanes/03_run_estimators.R E3         (one example)
# =====================================================================

suppressMessages({ library(SuperLearner); library(survey); library(surveyCV); library(tmle) })
REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
source(file.path(REPO_ROOT, "R", "estimators.R"))
ana_dir <- file.path(REPO_ROOT, "nhanes", "analytic")
res_dir <- file.path(REPO_ROOT, "nhanes", "results"); dir.create(res_dir, showWarnings = FALSE)

LIB <- c("SL.glm", "SL.earth", "SL.glmnet")     # flexible ensemble (gam/ranger -> sensitivity)
SEED <- 20260607L
args <- commandArgs(trailingOnly = TRUE)
ids  <- if (length(args)) args else c("E1", "E2", "E3", "E4")

run_one <- function(id) {
  obs  <- readRDS(file.path(ana_dir, paste0(id, "_imputed.rds")))
  covs <- attr(obs, "covs"); sub <- which(obs$inpop)
  # one-hot encode the (complete) domain covariates -> numeric matrix aligned to `sub`
  covdf <- droplevels(obs[sub, covs, drop = FALSE])
  mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
  stopifnot(nrow(mm) == length(sub))
  wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
  for (j in seq_along(wnames)) { obs[[wnames[j]]] <- NA_real_; obs[[wnames[j]]][sub] <- mm[, j] }
  obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED
  set.seed(SEED)
  t0 <- proc.time()[3]
  r  <- run_estimators(obs, learners = LIB, V_cf = 5L, inner_cv_folds = 5L,
                       W_cols = wnames, nest = TRUE, inpop = obs$inpop)
  secs <- round(proc.time()[3] - t0)
  rr <- r$results
  crit <- qt(0.975, pmax(1, rr$df))
  rr$lcl <- rr$b - crit * rr$se; rr$ucl <- rr$b + crit * rr$se
  rr$example <- id; rr$label <- attr(obs, "example")
  drow <- r$diagnostics$drow
  out <- list(results = rr, diag = drow, n_domain = length(sub),
              A_prev = mean(obs$A[sub]), Y_prev = mean(obs$Y[sub]),
              n_psu = length(unique(paste(obs$SDMVSTRA[sub], obs$SDMVPSU[sub]))),
              n_strata = length(unique(obs$SDMVSTRA[sub])), secs = secs, lib = LIB)
  dir.create(file.path(res_dir, id), showWarnings = FALSE)
  saveRDS(out, file.path(res_dir, id, "estimators_all_arms.rds"))
  cat(sprintf("\n=== %s (%s)  [%ds, n=%d] ===\n", id, attr(obs,"example"), secs, length(sub)))
  print(rr[, c("method","b","se","lcl","ucl","df")], row.names = FALSE)
  cat(sprintf("  g (Fully-Aware): [%.3f, %.3f]   g (CF): [%.3f, %.3f]\n",
              drow$g_fa_min, drow$g_fa_max, drow$g_cf_min, drow$g_cf_max))
  out
}

allres <- lapply(ids, run_one)

## combined cross-example CSV (only when running >1 example, so parallel single-id
## runs do not clobber each other; otherwise 06_compare assembles it from the rds)
if (length(ids) > 1) {
  comb <- do.call(rbind, lapply(allres, function(o) o$results))
  write.csv(comb[, c("example","label","method","b","se","lcl","ucl","df")],
            file.path(res_dir, "cross_example_summary.csv"), row.names = FALSE)
  cat("\nDONE. Combined results ->", file.path(res_dir, "cross_example_summary.csv"), "\n")
} else cat("\nDONE (single example).\n")
