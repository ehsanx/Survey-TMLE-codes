# =====================================================================
# Nhanes/R/02_impute.R  —  single imputation of covariates within each domain
#
# For each example, impute missing COVARIATES among the analytic sub-population
# (inpop) with mice (m=1; exposure A and outcome Y are observed by construction
# and are used as predictors but never imputed). The imputed covariate values are
# written back into the FULL MEC frame's inpop rows (out-of-pop rows are not used
# downstream); design vars + inpop are preserved. Single imputation under-states
# variance (stated caveat); multiple imputation + Rubin's rules is the principled
# extension (Web Appendix).
# =====================================================================

suppressMessages({ library(mice) })
REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
ana_dir <- file.path(REPO_ROOT, "Nhanes", "analytic")
SEED <- 20260607L

ids <- c("E1", "E2", "E3", "E4")
for (id in ids) {
  obs  <- readRDS(file.path(ana_dir, paste0(id, "_analytic.rds")))
  covs <- attr(obs, "covs")
  sub  <- which(obs$inpop)
  d    <- obs[sub, c(covs, "A", "Y")]
  # character -> factor for mice
  for (c in names(d)) if (is.character(d[[c]])) d[[c]] <- factor(d[[c]])
  miss_before <- colSums(is.na(d[, covs, drop = FALSE]))
  to_impute <- names(miss_before)[miss_before > 0]
  if (length(to_impute)) {
    set.seed(SEED)
    mi <- mice(d, m = 1, maxit = 5, printFlag = FALSE, seed = SEED)
    dc <- complete(mi, 1)
  } else dc <- d
  # write imputed covariates back into the full frame's inpop rows
  for (c in covs) obs[[c]][sub] <- dc[[c]]
  attr(obs, "covs") <- covs; attr(obs, "example") <- attr(obs, "example")
  saveRDS(obs, file.path(ana_dir, paste0(id, "_imputed.rds")))
  still_na <- sapply(covs, function(c) sum(is.na(obs[[c]][sub])))
  cat(sprintf("%s: imputed %d/%d covariates; residual NA in domain covs: %s\n",
              id, length(to_impute), length(covs),
              if (all(still_na == 0)) "none" else paste(names(still_na)[still_na>0], collapse=",")))
}
cat("\nDONE. Imputed frames saved to", ana_dir, "\n")
