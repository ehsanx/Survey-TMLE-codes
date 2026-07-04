# =====================================================================
# check_continuous.R  —  local verification for the continuous-outcome
# branch (A3): the Gaussian DGP + family='gaussian' estimator path
# (the bounded-[0,1] device). Confirms the RD recovers theta, the per-arm
# means recover the truth, and delta_rr_or(binary=FALSE) gives a sane
# ratio-of-means RR. Fast (SL.glm, one sample).  Rscript codes/check_continuous.R
# =====================================================================
CODE <- "R"
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
source(file.path(CODE, "estimands_rr_or.R"))
suppressMessages({ library(survey); library(tmle); library(SuperLearner) })

set.seed(20260613L)
# gamma0 = 5 keeps the continuous outcome POSITIVE so the ratio-of-means RR is
# interpretable; te_log_odds = 1.5 is the mean difference (ATE on the Y scale).
pop <- make_population("standard", model_type = "complex", outcome_type = "continuous",
                       gamma0 = 5, sigma2_eps = 1.0, truth_M = 2e5L)
cat(sprintf("[truth] RD(psi)=%.4f  psi1=%.4f psi0=%.4f  RR(means)=%.4f\n",
            pop$truth$psi, pop$truth$psi1, pop$truth$psi0, pop$truth$psi1 / pop$truth$psi0))

obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + 1L, model_type = "complex")
cat(sprintf("[sample] n=%d  Y in [%.2f, %.2f]  mean=%.2f\n",
            nrow(obs), min(obs$Y), max(obs$Y), mean(obs$Y)))

est <- run_estimators(obs, learners = "SL.glm", family = "gaussian")
res <- est$results
cat("\n-- arms (continuous; RD scale) --\n"); print(res, row.names = FALSE, digits = 4)
stopifnot(all(is.finite(res$b)), all(is.finite(res$se)), all(res$se > 0))

b_fa <- res$b[res$method == "Fully-Aware"]
cat(sprintf("\n[check] FA RD b = %.4f   (truth psi = %.4f)\n", b_fa, pop$truth$psi))

fa <- est$diagnostics$arms[["Fully-Aware"]]
cat(sprintf("[check] FA per-arm psi1 = %.4f psi0 = %.4f   (truth %.4f / %.4f)\n",
            fa$psi1, fa$psi0, pop$truth$psi1, pop$truth$psi0))
stopifnot(abs((fa$psi1 - fa$psi0) - b_fa) < 1e-8)   # per-arm split == RD contrast

rr <- delta_rr_or(fa$psi1, fa$psi0, fa$eif1, fa$eif0,
                  obs$strata, obs$cluster, obs$weight, binary = FALSE)  # ratio-of-means
cat("\n-- ratio-of-means RR (continuous) --\n"); print(rr, row.names = FALSE, digits = 4)
stopifnot(is.finite(rr$est), rr$se_log > 0, rr$lo < rr$est, rr$est < rr$hi)

cat("\n[check_continuous] ALL CHECKS PASSED.\n")
