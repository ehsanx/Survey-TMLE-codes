# =====================================================================
# Nhanes/arc_runs/R05_harmonized_floor/aggregate_R05.R
#   Combine R05 per-cell nh_*.rds -> results/arc/R05_summary.csv (+ diagnostics).
# Mirrors Nhanes/R/aggregate_nhanes.R but for the 6-arm R05 output and writes to
# the NON-DESTRUCTIVE results/arc/ location (never the locked nhanes_output/results).
# Usage:  Rscript Nhanes/arc_runs/R05_harmonized_floor/aggregate_R05.R
# =====================================================================

REPO <- Sys.getenv("REPO_ROOT", ".")
OUT  <- Sys.getenv("NH_OUT", file.path(REPO, "Nhanes", "nhanes_output", "arc_runs", "R05"))
RES  <- Sys.getenv("R05_RESULTS", file.path(REPO, "results", "arc"))
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^nh_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no nh_*.rds found in ", OUT, " -- run run.R (the array) first")
cat("aggregating", length(files), "R05 cell files from", OUT, "\n")
objs <- lapply(files, readRDS)

## (1) results table: example x rung x arm (incl the 2 NEW arms) ---------------
summary_tbl <- do.call(rbind, lapply(objs, `[[`, "summary"))
ord_m <- c("Fully-Aware", "Fully-Aware-h05",
           "Fully-Aware-CF", "Fully-Aware-CF-wOOF")          # FA, harmonized, CF(unwt), CF(wt)
ord_r <- c("L1_param", "L2_smooth", "L3_adaptive", "L4_aggressive")
summary_tbl <- summary_tbl[order(summary_tbl$example,
                                 match(summary_tbl$rung, ord_r),
                                 match(summary_tbl$method, ord_m)), ]
num <- sapply(summary_tbl, is.numeric)
st_print <- summary_tbl; st_print[num] <- round(st_print[num], 4)
write.csv(st_print, file.path(RES, "R05_summary.csv"), row.names = FALSE)

## (2) point-divergence table: the FLOOR-ROBUST decision quantities ------------
# Per example x rung: headline FA-vs-CF divergence, the (should-be-zero) floor
# effect on the FA point, and where weighted-OOF-CF lands (closer to FA or CF).
b_at <- function(s, m) { v <- s$b[s$method == m]; if (length(v)) v[1] else NA_real_ }
div_tbl <- do.call(rbind, lapply(objs, function(o) {
  s <- o$summary
  bFA <- b_at(s, "Fully-Aware"); bFAh <- b_at(s, "Fully-Aware-h05")
  bCF <- b_at(s, "Fully-Aware-CF"); bCFw <- b_at(s, "Fully-Aware-CF-wOOF")
  data.frame(example = o$example, label = o$label, rung = o$rung,
             b_FA = bFA, b_FA_h05 = bFAh, b_CF = bCF, b_CF_wOOF = bCFw,
             d_FA_CF = abs(bFA - bCF),            # HEADLINE divergence (floor-robust)
             d_FA_FAh05 = abs(bFA - bFAh),        # floor effect on the POINT (~0 by construction)
             d_CFwOOF_CF = abs(bCFw - bCF),       # de-weighting effect within cross-fit
             d_CFwOOF_FA = abs(bCFw - bFA),       # weighted-OOF vs single-fit FA
             # which side does weighted-OOF-CF sit closer to?
             cfwoof_closer_to = ifelse(abs(bCFw - bCF) <= abs(bCFw - bFA), "CF", "FA"))
}))
div_tbl <- div_tbl[order(div_tbl$example, match(div_tbl$rung, ord_r)), ]
dn <- sapply(div_tbl, is.numeric); dvp <- div_tbl; dvp[dn] <- round(dvp[dn], 4)
write.csv(dvp, file.path(RES, "R05_point_divergence.csv"), row.names = FALSE)

## (3) diagnostics table: example x rung (floor sensitivity + overlap) ---------
diag_tbl <- do.call(rbind, lapply(objs, function(o) with(o$diagnostics, data.frame(
  example = o$example, label = o$label, rung = o$rung,
  fa_gbound = fa_gbound, g_oof_bound = g_oof_bound,
  n_domain = n_domain, n_psu = n_psu, n_strata = n_strata, design_df = design_df,
  A_prev = A_prev, Y_prev = Y_prev,
  deff_eif = deff_clust, deff_eif_h05 = deff_clust_h05,        # floor-sensitivity on DEFF
  icc_eif = icc_eif, icc_eif_h05 = icc_eif_h05,
  weight_cv = weight_cv, min_psu_stratum = min_psu_stratum,
  g_fa_min = min(g_fa), g_fa_max = max(g_fa), g_fa_near_bound = g_fa_near_bound,
  g_cf_min = min(g_cf), g_cf_max = max(g_cf),                  # UNWEIGHTED-OOF overlap
  g_cfw_min = min(g_cfw), g_cfw_max = max(g_cfw),              # WEIGHTED-OOF overlap
  eps_fa = drow$eps_fa, eps_cf = drow$eps_cf, eps_cfw = drow$eps_cfw,
  cf_V_eff = drow$cf_V_eff))))
diag_tbl <- diag_tbl[order(diag_tbl$example, match(diag_tbl$rung, ord_r)), ]
dgn <- sapply(diag_tbl, is.numeric); dp <- diag_tbl; dp[dgn] <- round(dp[dgn], 4)
write.csv(dp, file.path(RES, "R05_diagnostics.csv"), row.names = FALSE)

## (4) combined object (per-split + g vectors for figures) ---------------------
saveRDS(list(summary = summary_tbl, point_divergence = div_tbl, diagnostics = diag_tbl,
             per_split = do.call(rbind, lapply(objs, `[[`, "per_split")),
             cells = objs),
        file.path(RES, "R05_combined.rds"))

cat("\n========== R05 RESULTS (example x rung x arm) ==========\n")
print(st_print, row.names = FALSE)
cat("\n========== R05 POINT DIVERGENCE (floor-robust decision quantities) ==========\n")
print(dvp[, c("example","rung","b_FA","b_FA_h05","b_CF","b_CF_wOOF",
              "d_FA_CF","d_FA_FAh05","cfwoof_closer_to")], row.names = FALSE)
cat("\nSaved:\n -", file.path(RES, "R05_summary.csv"),
    "\n -", file.path(RES, "R05_point_divergence.csv"),
    "\n -", file.path(RES, "R05_diagnostics.csv"),
    "\n -", file.path(RES, "R05_combined.rds"), "\n")
