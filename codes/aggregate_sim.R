# =====================================================================
# aggregate_sim.R  —  combine per-task RDS from run_sim.R into summaries
# Usage:  Rscript aggregate_sim.R   (reads SIM_OUT, writes to RESULTS_DIR/DATA_RESULTS)
# =====================================================================
# locate config.R robustly: SIM_CODE env, else repo-relative, else Windows default
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("codes", "codes")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from the repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))
# OUT defaults to DATA_INTERMEDIATE, which (post config.R fix) points at
# .../survey-tmle2-data/sim_output/intermediate on Windows and the ARC scratch path on ARC.
OUT <- Sys.getenv("SIM_OUT", DATA_INTERMEDIATE)

files <- list.files(OUT, pattern = "^sim_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no sim_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))
diags <- do.call(rbind, lapply(objs, function(o)
  if (!is.null(o$diagnostics)) cbind(scenario = o$scenario, rung = o$rung, o$diagnostics)))

z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows, list(rows$scenario, rows$rung, rows$method), function(d) {
  if (is.null(d)) return(NULL)
  Psi <- d$Psi[1]; crit <- z_or_t(d$df)
  cov <- mean(abs(d$b - Psi) <= crit * d$se)
  data.frame(scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
             n_reps = nrow(d), Psi = Psi,
             bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
             se_ratio = mean(d$se) / sd(d$b),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
             deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif))
}))
rownames(agg) <- NULL
ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
agg <- agg[order(agg$scenario, agg$rung, match(agg$method, ord_m)), ]  # L1_/L2_... sort in ladder order

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
print(agg_print, row.names = FALSE)

# ---- diagnostic summary per (scenario,rung): targeting + positivity -----------
diag_sum <- NULL
if (!is.null(diags)) {
  diag_sum <- do.call(rbind, by(diags, list(diags$scenario, diags$rung), function(d) {
    if (is.null(d)) return(NULL)
    data.frame(scenario = d$scenario[1], rung = d$rung[1], n_reps = nrow(d),
               eps_fa = mean(d$eps_fa), g_fa_min = mean(d$g_fa_min), g_fa_max = mean(d$g_fa_max),
               eps_cf = mean(d$eps_cf), g_cf_min = mean(d$g_cf_min), g_cf_max = mean(d$g_cf_max),
               cf_V_eff = d$cf_V_eff[1], deff = mean(d$deff), icc_eif = mean(d$icc_eif))
  }))
  rownames(diag_sum) <- NULL
  diag_sum <- diag_sum[order(diag_sum$scenario, diag_sum$rung), ]
}

saveRDS(list(per_rep = rows, summary = agg, per_rep_diag = diags, diag_summary = diag_sum),
        file.path(DATA_RESULTS, "sim_full_combined.rds"))
write.csv(agg_print, file.path(RESULTS_DIR, "sim_full_summary.csv"), row.names = FALSE)
if (!is.null(diag_sum)) {
  ds <- diag_sum; dn <- sapply(ds, is.numeric); ds[dn] <- round(ds[dn], 4)
  write.csv(ds, file.path(RESULTS_DIR, "sim_full_diagnostics.csv"), row.names = FALSE)
}
cat("\nSaved:\n -", file.path(RESULTS_DIR, "sim_full_summary.csv"),
    "\n -", file.path(RESULTS_DIR, "sim_full_diagnostics.csv"),
    "\n -", file.path(DATA_RESULTS, "sim_full_combined.rds"), "\n")
