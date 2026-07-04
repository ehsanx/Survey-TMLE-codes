# =====================================================================
# check_rr_or.R  —  local verification for the per-arm EIF exposure (A1)
# and the delta_rr_or() post-processor (A2). Confirms the estimator change
# is NON-BREAKING (the per-arm split reproduces the contrast scalar + EIF
# byte-for-byte) and that the delta-method RR/OR are sane on one sim cell.
# Fast (SL.glm, one sample).   Rscript codes/check_rr_or.R
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
pop <- make_population("standard", model_type = "complex", truth_M = 2e5L)
obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + 1L, model_type = "complex")
est <- run_estimators(obs, learners = "SL.glm")

# ---- (1) NON-BREAKING identity checks --------------------------------------
res  <- est$results
b_fa <- res$b[res$method == "Fully-Aware"]
arms <- est$diagnostics$arms
fa   <- arms[["Fully-Aware"]]
id_psi <- abs((fa$psi1 - fa$psi0) - b_fa)
id_eif <- max(abs((fa$eif1 - fa$eif0) - est$diagnostics$eif_fa))
cat(sprintf("[check] (psi1 - psi0) vs FA b   : |diff| = %.3e\n", id_psi))
cat(sprintf("[check] (eif1 - eif0) vs eif_fa : max|diff| = %.3e\n", id_eif))
stopifnot(id_psi < 1e-10, id_eif < 1e-10)
cat("[check] PASS: per-arm decomposition reproduces the contrast scalar + EIF exactly.\n\n")

# ---- (2) delta_rr_or sanity on FA + CF -------------------------------------
for (m in c("Fully-Aware", "Fully-Aware-CF")) {
  a  <- arms[[m]]
  rr <- delta_rr_or(a$psi1, a$psi0, a$eif1, a$eif0,
                    obs$strata, obs$cluster, obs$weight, clustered = TRUE)
  cat(sprintf("-- %s --  psi1=%.4f psi0=%.4f\n", m, a$psi1, a$psi0))
  print(rr, row.names = FALSE, digits = 4)
  stopifnot(all(is.finite(rr$est)), all(rr$se_log > 0),
            all(rr$lo < rr$est), all(rr$est < rr$hi))
  cat("\n")
}

# ---- (3) truth for context (per-arm population means -> RR/OR) -------------
p1 <- mean(pop$pop$Q1_int); p0 <- mean(pop$pop$Q0_int)
cat(sprintf("[truth] psi1=%.4f psi0=%.4f  RD=%.4f  RR=%.4f  OR=%.4f\n",
            p1, p0, p1 - p0, p1 / p0, (p1 / (1 - p1)) / (p0 / (1 - p0))))
cat("\n[check_rr_or] ALL CHECKS PASSED.\n")
