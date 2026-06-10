# =====================================================================
# aggregate.R  —  combine R02_largem_sweep per-task RDS into the m-sweep summary
# Usage:  Rscript aggregate.R        (after all array tasks finish)
# Reads:  sim_output/arc_runs/R02_largem_sweep/sim_m*_chunk*.rds
# Writes: results/arc/R02_largem_sweep_summary.csv   (per m x rung x method)
#         results/arc/R02_largem_sweep_diag.csv      (targeting/positivity per cell)
#         sim_output/arc_runs/R02_largem_sweep/R02_combined.rds  (full detail)
# NON-DESTRUCTIVE: writes only under results/arc and the run's own output dir.
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("codes", "codes")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

RUN_ID <- "R02_largem_sweep"
OUT    <- file.path(DATA_ROOT, "arc_runs", RUN_ID)
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RESULTS)) dir.create(ARC_RESULTS, recursive = TRUE, showWarnings = FALSE)

# run-local summary helper (matches aggregate_sim.R column conventions)
RUN_DIR <- file.path(REPO_ROOT, "codes", "arc_runs", RUN_ID)
source(file.path(RUN_DIR, "helpers.R"))

# read real (non-smoke) task files
files <- list.files(OUT, pattern = "^sim_m\\d+_chunk\\d+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no sim_m*.rds found in ", OUT, " (run the array first)")
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o) o$results))
diag <- do.call(rbind, lapply(objs, function(o) o$diagnostics))

# ---- per (m, rung, method) summary ------------------------------------------
agg <- summarise_sweep(rows)
# sort: rung (ladder/contrast order) -> m -> method order from aggregate_sim.R
ord_r <- c("L1_param", "L3_adaptive", "RF_shallow", "L4_aggressive")
ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
agg <- agg[order(match(agg$rung, ord_r), agg$base_m, match(agg$method, ord_m)), ]
rownames(agg) <- NULL

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
cat("\n==== R02 m-sweep summary (focus: Fully-Aware-CF se_ratio vs m) ====\n")
print(agg_print, row.names = FALSE)

# ---- headline: CF se_ratio / coverage vs m for L4 and the shallow-RF contrast -
cf <- agg[agg$method == "Fully-Aware-CF" &
          agg$rung %in% c("L4_aggressive", "RF_shallow"), ]
cf <- cf[order(cf$rung, cf$base_m), ]
cat("\n==== Fully-Aware-CF: se_ratio & coverage vs m (decision view) ====\n")
print(cf[, c("rung", "base_m", "m_total", "se_ratio", "coverage", "mcse_cov", "bias")],
      row.names = FALSE)

# ---- bias track at L1/L3 vs m (O(1/m) Hajek-ratio check) ---------------------
bt <- agg[agg$method == "Fully-Aware" & agg$rung %in% c("L1_param", "L3_adaptive"), ]
bt <- bt[order(bt$rung, bt$base_m), ]
cat("\n==== Fully-Aware bias vs m at L1/L3 (O(1/m) Hajek-ratio check) ====\n")
print(bt[, c("rung", "base_m", "m_total", "bias", "emp_sd", "coverage")], row.names = FALSE)

# ---- diagnostic summary per (m, rung): targeting + OOF positivity ------------
diag_sum <- NULL
if (!is.null(diag)) {
  diag_sum <- do.call(rbind, by(diag, list(diag$base_m, diag$rung_label), function(d) {
    if (is.null(d) || !nrow(d)) return(NULL)
    data.frame(base_m = d$base_m[1], rung = d$rung_label[1], n_reps = nrow(d),
               eps_fa = mean(d$eps_fa), g_fa_min = mean(d$g_fa_min), g_fa_max = mean(d$g_fa_max),
               eps_cf = mean(d$eps_cf), g_cf_min = mean(d$g_cf_min), g_cf_max = mean(d$g_cf_max),
               cf_V_eff = d$cf_V_eff[1], deff = mean(d$deff), icc_eif = mean(d$icc_eif),
               stringsAsFactors = FALSE)
  }))
  rownames(diag_sum) <- NULL
  diag_sum <- diag_sum[order(match(diag_sum$rung, ord_r), diag_sum$base_m), ]
}

# ---- write (NON-DESTRUCTIVE; results/arc + run dir only) ---------------------
saveRDS(list(per_rep = rows, summary = agg, per_rep_diag = diag, diag_summary = diag_sum),
        file.path(OUT, "R02_combined.rds"))
write.csv(agg_print, file.path(ARC_RESULTS, "R02_largem_sweep_summary.csv"), row.names = FALSE)
if (!is.null(diag_sum)) {
  ds <- diag_sum; dn <- sapply(ds, is.numeric); ds[dn] <- round(ds[dn], 4)
  write.csv(ds, file.path(ARC_RESULTS, "R02_largem_sweep_diag.csv"), row.names = FALSE)
}
cat("\nSaved:\n -", file.path(ARC_RESULTS, "R02_largem_sweep_summary.csv"),
    "\n -", file.path(ARC_RESULTS, "R02_largem_sweep_diag.csv"),
    "\n -", file.path(OUT, "R02_combined.rds"), "\n")
