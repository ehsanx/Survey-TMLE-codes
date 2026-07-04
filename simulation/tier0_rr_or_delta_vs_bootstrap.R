# =====================================================================
# tier0_rr_or_delta_vs_bootstrap.R
#   contribution.md  §6 Tier-0 executable evidence for the RR/OR corollary
#   (cor:rrora): on ONE simulation cell, the delta-method design SE of
#   log-RR and log-OR must agree with an INDEPENDENT survey bootstrap of
#   the arms (resample whole PSUs within strata; recompute RR/OR
#   NON-linearly per resample). A SEEDED-ERROR CONTROL plants a wrong
#   gradient; the harness must FLAG the resulting disagreement.
#
#   PASS  = correct delta SE within TOL of the bootstrap SE for BOTH RR,OR
#           AND the planted-error control is detected (disagrees).
#   Rscript codes/tier0_rr_or_delta_vs_bootstrap.R
# =====================================================================
CODE <- "R"
source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R")); source(file.path(CODE, "estimands_rr_or.R"))
suppressMessages({ library(survey); library(tmle); library(SuperLearner) })

TOL <- 0.10          # relative agreement tolerance (delta vs bootstrap)
B   <- 4000          # bootstrap resamples
set.seed(20260613L)

# ---- one cell: one binary population, one sample, one fit --------------------
pop <- make_population("standard", model_type = "complex", outcome_type = "binary", truth_M = 2e5L)
obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + 7L, model_type = "complex")
est <- run_estimators(obs, learners = "SL.glm")
a   <- est$diagnostics$arms[["Fully-Aware"]]      # psi1, psi0, eif1 (D1), eif0 (D0)
psi1 <- a$psi1; psi0 <- a$psi0; D1 <- a$eif1; D0 <- a$eif0
cat(sprintf("[cell] standard/binary/L1, n=%d  psi1=%.4f psi0=%.4f  RR=%.4f OR=%.4f\n",
            nrow(obs), psi1, psi0, psi1/psi0, (psi1/(1-psi1))/(psi0/(1-psi0))))

# ---- (1) delta-method design SE (the quantity under test) -------------------
del <- delta_rr_or(psi1, psi0, D1, D0, obs$strata, obs$cluster, obs$weight, binary = TRUE)
se_delta_RR <- del$se_log[del$estimand == "RR"]
se_delta_OR <- del$se_log[del$estimand == "OR"]

# ---- (2) INDEPENDENT survey bootstrap of the arms ---------------------------
# Resample whole PSUs WITH REPLACEMENT within each stratum (m_h per stratum),
# perturb the arm means by the resampled per-PSU influence totals, and form
# RR/OR NON-linearly each resample. With-replacement resampling of m_h from m_h
# is variance-deflated by (m_h-1)/m_h; apply the Rao-Wu factor sqrt(m_h/(m_h-1)).
w   <- obs$weight; Nhat <- sum(w)
key <- paste(obs$strata, obs$cluster, sep = "\r")
s1  <- tapply(w * D1, key, sum); s0 <- tapply(w * D0, key, sum)   # per-PSU influence totals
str_of_key <- tapply(obs$strata, key, `[`, 1)
psu_by_str <- split(names(s1), str_of_key[names(s1)])
mh         <- vapply(psu_by_str, length, integer(1))
rw_factor  <- sqrt(mean(mh / (mh - 1)))                            # Rao-Wu rescaling (m_h=6 -> 1.095)

logRR <- numeric(B); logOR <- numeric(B)
for (b in seq_len(B)) {
  d1 <- 0; d0 <- 0
  for (h in names(psu_by_str)) {
    ks  <- psu_by_str[[h]]; sel <- sample(ks, length(ks), replace = TRUE)
    d1  <- d1 + sum(s1[sel]); d0 <- d0 + sum(s0[sel])
  }
  p1b <- min(max(psi1 + d1 / Nhat, 1e-6), 1 - 1e-6)
  p0b <- min(max(psi0 + d0 / Nhat, 1e-6), 1 - 1e-6)
  logRR[b] <- log(p1b) - log(p0b)
  logOR[b] <- qlogis(p1b) - qlogis(p0b)
}
se_boot_RR <- sd(logRR) * rw_factor
se_boot_OR <- sd(logOR) * rw_factor

# ---- (3) seeded-error CONTROL: a wrong gradient (treat psi0 as KNOWN) --------
# log-RR gradient should be (1/psi1, -1/psi0); the planted error drops the arm-0
# term -> ignores arm-0 variance -> SE too small -> must DISAGREE with bootstrap.
cc  <- vcov(svymean(~e1 + e0, svydesign(ids = ~cluster, strata = ~strata, weights = ~weight,
            data = data.frame(e1 = D1, e0 = D0, cluster = obs$cluster, strata = obs$strata, weight = w))))
g_bad <- c(1/psi1, 0)
se_delta_RR_bad <- sqrt(as.numeric(crossprod(g_bad, cc %*% g_bad)))

rel <- function(a, b) abs(a - b) / b
cat("\n--- delta-method vs bootstrap (log scale SE) ---\n")
cat(sprintf("  RR : delta=%.4f  boot=%.4f  rel.diff=%.1f%%\n", se_delta_RR, se_boot_RR, 100*rel(se_delta_RR, se_boot_RR)))
cat(sprintf("  OR : delta=%.4f  boot=%.4f  rel.diff=%.1f%%\n", se_delta_OR, se_boot_OR, 100*rel(se_delta_OR, se_boot_OR)))
cat(sprintf("  CONTROL (wrong RR gradient): delta=%.4f vs boot=%.4f  rel.diff=%.1f%%\n",
            se_delta_RR_bad, se_boot_RR, 100*rel(se_delta_RR_bad, se_boot_RR)))

pass_RR   <- rel(se_delta_RR, se_boot_RR) < TOL
pass_OR   <- rel(se_delta_OR, se_boot_OR) < TOL
ctrl_flag <- rel(se_delta_RR_bad, se_boot_RR) > TOL          # control must be DETECTED (disagree)
verdict   <- pass_RR && pass_OR && ctrl_flag

cat(sprintf("\n[TIER0-GATE] %s  (RR agree=%s, OR agree=%s, control detected=%s; TOL=%.0f%%, B=%d)\n",
            if (verdict) "PASS" else "FAIL", pass_RR, pass_OR, ctrl_flag, 100*TOL, B))

dir.create(file.path(RESULTS_DIR, "tier0"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(psi1 = psi1, psi0 = psi0,
             se_delta_RR = se_delta_RR, se_boot_RR = se_boot_RR,
             se_delta_OR = se_delta_OR, se_boot_OR = se_boot_OR,
             se_delta_RR_control = se_delta_RR_bad, rw_factor = rw_factor,
             B = B, TOL = TOL, pass = verdict),
        file.path(RESULTS_DIR, "tier0", "tier0_rr_or_delta_vs_bootstrap.rds"))
if (!verdict) quit(status = 1)
