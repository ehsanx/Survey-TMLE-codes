# =====================================================================
# Nhanes/R/03b_run_continuous_E1.R  —  A8: a CONTINUOUS-outcome NHANES example.
#
# E1 with the outcome taken CONTINUOUS: Y = BMXBMI (body-mass index) instead of
# the obesity indicator 1{BMXBMI >= 30}. Same exposure (short sleep A), same
# covariates, same complex survey design. Demonstrates the Web-Appendix RR/OR
# corollary's bounded-continuous extension on real survey data: the design-based
# ATE (mean-BMI difference) + CI for the single-fit (Fully-Aware) and cross-fit
# (Fully-Aware-CF) arms, plus the ratio-of-means risk ratio via the delta method.
#
# BMXBMI was not retained in E1_analytic (replaced by the binary Y), so it is
# merged back by SEQN from the cached raw body-measures files (BMX_*). Mirrors
# 03_run_estimators.R for the one-hot/design/domain handling; only family and the
# outcome differ.   Rscript Nhanes/R/03b_run_continuous_E1.R
# =====================================================================
suppressMessages({ library(SuperLearner); library(survey); library(surveyCV); library(tmle) })
REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
source(file.path(REPO_ROOT, "R", "estimators.R"))
source(file.path(REPO_ROOT, "R", "estimands_rr_or.R"))
ana_dir <- file.path(REPO_ROOT, "nhanes", "analytic")
raw_dir <- file.path(REPO_ROOT, "nhanes", "raw")
res_dir <- file.path(REPO_ROOT, "nhanes", "results", "E1_continuous")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
LIB <- c("SL.glm", "SL.earth", "SL.glmnet"); SEED <- 20260613L

# ---- E1 imputed covariates + merge the continuous BMI outcome by SEQN ----------
obs  <- readRDS(file.path(ana_dir, "E1_imputed.rds"))
covs <- attr(obs, "covs")
bmx  <- do.call(rbind, lapply(list.files(raw_dir, "^BMX_.*\\.rds$", full.names = TRUE),
                              function(f) { d <- readRDS(f); d[, c("SEQN", "BMXBMI")] }))
obs$BMI   <- bmx$BMXBMI[match(obs$SEQN, bmx$SEQN)]
obs$Y     <- obs$BMI                               # CONTINUOUS outcome
obs$inpop <- obs$inpop & !is.na(obs$BMI)           # domain: E1 eligibility AND BMI observed
sub <- which(obs$inpop)
cat(sprintf("E1-continuous (short sleep -> BMI): n_domain=%d  BMI in [%.1f, %.1f] mean %.1f  A_prev=%.3f\n",
            length(sub), min(obs$BMI[sub]), max(obs$BMI[sub]), mean(obs$BMI[sub]), mean(obs$A[sub])))

# ---- one-hot the (imputed) domain covariates -> numeric design matrix ----------
covdf <- droplevels(obs[sub, covs, drop = FALSE])
mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
for (j in seq_along(wnames)) { obs[[wnames[j]]] <- NA_real_; obs[[wnames[j]]][sub] <- mm[, j] }
obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED

# ---- five arms with gaussian outcome family (bounded-[0,1] device inside tmle) --
set.seed(SEED); t0 <- proc.time()[3]
r <- run_estimators(obs, learners = LIB, V_cf = 5L, inner_cv_folds = 5L, family = "gaussian",
                    W_cols = wnames, nest = TRUE, inpop = obs$inpop)
secs <- round(proc.time()[3] - t0)
rr <- r$results; crit <- qt(0.975, pmax(1, rr$df))
rr$lcl <- rr$b - crit * rr$se; rr$ucl <- rr$b + crit * rr$se
cat(sprintf("\n=== E1-continuous arms (ATE = mean-BMI difference)  [%ds, n=%d] ===\n", secs, length(sub)))
print(rr[, c("method", "b", "se", "lcl", "ucl", "df")], row.names = FALSE, digits = 4)

# ---- ratio-of-means RR (BMI > 0) for the single-fit + cross-fit arms -----------
rrtab <- do.call(rbind, lapply(c("Fully-Aware", "Fully-Aware-CF"), function(m) {
  a <- r$diagnostics$arms[[m]]
  x <- delta_rr_or(a$psi1, a$psi0, a$eif1, a$eif0, obs$strata, obs$cluster, obs$weight,
                   binary = FALSE, nest = TRUE, inpop = obs$inpop); x$method <- m; x }))
cat("\n--- ratio-of-means RR (mean BMI: short sleep / normal sleep) ---\n")
print(rrtab[, c("method", "estimand", "est", "lo", "hi", "se_log", "psi1", "psi0")],
      row.names = FALSE, digits = 4)

saveRDS(list(results = rr, rr_means = rrtab, n_domain = length(sub), secs = secs,
             A_prev = mean(obs$A[sub]), Y_mean = mean(obs$Y[sub]), lib = LIB),
        file.path(res_dir, "estimators_continuous.rds"))
write.csv(rr[, c("method", "b", "se", "lcl", "ucl", "df")],
          file.path(res_dir, "E1_continuous_arms.csv"), row.names = FALSE)
cat("\nsaved ->", res_dir, "\n")
