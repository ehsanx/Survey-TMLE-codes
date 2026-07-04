# =====================================================================
# helpers.R  —  R02_largem_sweep run-local helpers (NON-DESTRUCTIVE)
#
# Purpose: sweep the number of SAMPLED PSUs-per-stratum (m = base_m) while
# holding the finite population + DGP + library FIXED, to test whether the
# L4 cross-fit over-coverage (se_ratio ~ 1.44) shrinks as m grows (finite-m,
# conservative-consistent per Theorem 2) or is structural.
#
# IMPORTANT: this file does NOT edit R/dgp.R or R/estimators.R. It
# wraps draw_sample()'s EXISTING `base_m` argument (PSUs sampled per stratum,
# stage-1 SRS-WOR). Increasing base_m raises the total sampled PSU count
#   m_total = H * base_m
# and the design degrees of freedom
#   df_design = H*base_m - H = H*(base_m - 1),
# WITHOUT touching the population geometry (H, J_per_stratum, M_per_psu) or
# the units-per-PSU (n0_h, controlled by base_n0) -- so the per-PSU cluster
# size and weight structure are held fixed and only the # of independent PSUs
# (the asymptotic dimension in Theorem 2) increases.
#
# Constraint: base_m must be <= J_per_stratum (the # of population PSUs per
# stratum). For the 'standard' scenario J_per_stratum = 60, so base_m in
# {6, 12, 20, 30} (= baseline, 2x, ~3.3x, 5x) is valid (stage-1 SRS-WOR).
#
# Reuses the engine building blocks via run_estimators() + deff_clust(); the
# ONLY new logic here is (a) the per-rep wrapper that forwards base_m, and
# (b) a per-arm summary that matches aggregate_sim.R column conventions.
# =====================================================================

# ---- one replication at a chosen sampled-PSU-count m = base_m ---------------
# Mirrors run_sim.R::one_rep but forwards base_m to draw_sample (the m knob)
# and keeps base_n0 fixed (default for the scenario) so units-per-PSU is held.
one_rep_m <- function(i, pop, learners, base_m) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i,
                     model_type = "complex", base_m = base_m)
  est <- run_estimators(obs, learners = learners)
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results,
                    deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design,
                    sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow,
                    deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}

# ---- per-(m, rung, method) summary matching aggregate_sim.R conventions ------
# bias = mean(b) - Psi ; emp_sd = sd(b) ; se_ratio = mean(se)/emp_sd ;
# coverage with t df ; mcse_cov = sqrt(cov(1-cov)/nreps). Adds the swept knobs
# (base_m, m_total = H*base_m, mean df_design) and the rung label as columns so
# the m-sweep is plottable directly.
summarise_sweep <- function(rows) {
  z_or_t <- function(df) qt(0.975, pmax(1, df))
  do.call(rbind, by(rows, list(rows$base_m, rows$rung_label, rows$method),
    function(d) {
      if (is.null(d) || !nrow(d)) return(NULL)
      Psi  <- d$Psi[1]
      crit <- z_or_t(d$df)
      cov  <- mean(abs(d$b - Psi) <= crit * d$se)
      data.frame(
        scenario   = d$scenario[1],
        rung       = d$rung_label[1],
        base_m     = d$base_m[1],            # sampled PSUs per stratum (the m knob)
        m_total    = d$m_total[1],           # total sampled PSUs = H * base_m
        method     = d$method[1],
        n_reps     = nrow(d),
        Psi        = Psi,
        bias       = mean(d$b) - Psi,
        emp_sd     = sd(d$b),
        mean_se    = mean(d$se),
        se_ratio   = mean(d$se) / sd(d$b),
        coverage   = cov,
        mcse_cov   = sqrt(cov * (1 - cov) / nrow(d)),
        mean_df    = mean(d$df),
        deff_clust = mean(d$deff),
        icc_eif    = mean(d$icc_eif),
        stringsAsFactors = FALSE
      )
    }))
}
