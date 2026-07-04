# =====================================================================
# Nhanes/R/aggregate_nhanes.R  —  combine per-task nh_*.rds into result tables
# Usage: Rscript Nhanes/R/aggregate_nhanes.R   (reads NH_OUT, writes NH_RESULTS)
# =====================================================================

OUT <- Sys.getenv("NH_OUT", "Nhanes/nhanes_output/intermediate")
RES <- Sys.getenv("NH_RESULTS", "Nhanes/nhanes_output/results")
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^nh_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no nh_*.rds found in ", OUT, " -- run nhanes_arc.R (the array) first")
cat("aggregating", length(files), "task files from", OUT, "\n")
objs <- lapply(files, readRDS)

## (1) the results table: example x rung x arm (estimate, SE, 95% CI, df, split-SD)
summary_tbl <- do.call(rbind, lapply(objs, `[[`, "summary"))
ord_m <- c("Non-Aware","Partially-Aware","Fully-Aware","Fully-Aware-CV","Fully-Aware-CF")
ord_r <- c("L1_param","L2_smooth","L3_adaptive","L4_aggressive")
summary_tbl <- summary_tbl[order(summary_tbl$example,
                                 match(summary_tbl$rung, ord_r),
                                 match(summary_tbl$method, ord_m)), ]
num <- sapply(summary_tbl, is.numeric)
st_print <- summary_tbl; st_print[num] <- round(st_print[num], 4)
write.csv(st_print, file.path(RES, "nhanes_results_summary.csv"), row.names = FALSE)

## (2) the diagnostics table: example x rung
diag_tbl <- do.call(rbind, lapply(objs, function(o) with(o$diagnostics, data.frame(
  example = o$example, label = o$label, rung = o$rung,
  n_domain = n_domain, n_psu = n_psu, n_strata = n_strata, design_df = design_df,
  A_prev = A_prev, Y_prev = Y_prev, deff_eif = deff_clust, icc_eif = icc_eif,
  weight_cv = weight_cv, min_psu_stratum = min_psu_stratum,
  g_fa_min = min(g_fa), g_fa_max = max(g_fa), g_fa_near_bound = g_fa_near_bound,
  g_cf_min = min(g_cf), g_cf_max = max(g_cf),
  eps_fa = drow$eps_fa, eps_cf = drow$eps_cf, cf_V_eff = drow$cf_V_eff))))
diag_tbl <- diag_tbl[order(diag_tbl$example, match(diag_tbl$rung, ord_r)), ]
dn <- sapply(diag_tbl, is.numeric); dp <- diag_tbl; dp[dn] <- round(dp[dn], 4)
write.csv(dp, file.path(RES, "nhanes_diagnostics.csv"), row.names = FALSE)

## (3) combined object (everything incl per-split + g vectors for figures)
saveRDS(list(summary = summary_tbl, diagnostics = diag_tbl,
             per_split = do.call(rbind, lapply(objs, `[[`, "per_split")),
             cells = objs),
        file.path(RES, "nhanes_combined.rds"))

cat("\n================ NHANES RESULTS (example x rung x arm) ================\n")
print(st_print, row.names = FALSE)
cat("\n================ DIAGNOSTICS (example x rung) ================\n")
print(dp[, c("example","rung","n_domain","design_df","A_prev","Y_prev","deff_eif",
             "g_fa_min","g_fa_near_bound","g_cf_min")], row.names = FALSE)
cat("\nSaved:\n -", file.path(RES, "nhanes_results_summary.csv"),
    "\n -", file.path(RES, "nhanes_diagnostics.csv"),
    "\n -", file.path(RES, "nhanes_combined.rds"), "\n")
