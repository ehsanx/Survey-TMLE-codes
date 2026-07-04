# =====================================================================
# aggregate.R  —  combine R07_informative per-task RDS into the rho-sweep summary
# Usage:  Rscript simulation/enhancements/R07_informative/aggregate.R
#   reads R07_OUT (default sim_output/arc_runs/R07_informative)
#   writes results/R07_informative_summary.csv  (+ combined RDS alongside)
#
# Column conventions match aggregate_sim.R EXACTLY where shared:
#   method,b,se,df,bias,emp_sd,mean_se,se_ratio,coverage,mcse_cov,deff_clust,icc_eif
# plus the run-specific key column `rho`. Coverage uses a t-reference CI with the
# design df (qt), identical to aggregate_sim.R's z_or_t().
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R07_OUT", file.path(DATA_ROOT, "arc_runs", "R07_informative"))
files <- list.files(OUT, pattern = "^r07_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no r07_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(rho = o$rho, rung = o$rung, Psi = o$Psi, o$results[, setdiff(names(o$results), c("rho","rung"))])))

z_or_t <- function(df) qt(0.975, pmax(1, df))     # matches aggregate_sim.R
agg <- do.call(rbind, by(rows, list(rows$rho, rows$rung, rows$method), function(d) {
  if (is.null(d)) return(NULL)
  Psi <- d$Psi[1]; crit <- z_or_t(d$df)
  cov <- mean(abs(d$b - Psi) <= crit * d$se)
  data.frame(rho = d$rho[1], rung = d$rung[1], method = d$method[1],
             n_reps = nrow(d), Psi = Psi,
             bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
             se_ratio = mean(d$se) / sd(d$b),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
             deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif),
             w_cv = mean(d$w_cv), df_design = mean(d$df_design))
}))
rownames(agg) <- NULL
ord_m <- c("Fully-Aware-CF-unwt", "Fully-Aware-CF-wt")
agg <- agg[order(agg$rung, agg$rho, match(agg$method, ord_m)), ]

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
print(agg_print, row.names = FALSE)

ARC_RES <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RES)) dir.create(ARC_RES, recursive = TRUE, showWarnings = FALSE)
write.csv(agg_print, file.path(ARC_RES, "R07_informative_summary.csv"), row.names = FALSE)
saveRDS(list(per_rep = rows, summary = agg), file.path(ARC_RES, "R07_informative_combined.rds"))

# ---- decision-rule readout: does WEIGHTED bias ALSO move with rho? -----------
# de-risk: |bias_wt(rho=max)| - |bias_wt(rho=0)| should be NON-trivially > 0,
# otherwise the curve is flat from a mis-implemented DGP (selection not biting).
cat("\n--- DECISION-RULE READOUT (selection-bites check) ---\n")
for (rg in unique(agg$rung)) {
  a <- agg[agg$rung == rg, ]
  for (mth in ord_m) {
    am <- a[a$method == mth, ]; am <- am[order(am$rho), ]
    if (nrow(am) < 2) next
    b0 <- abs(am$bias[am$rho == min(am$rho)])
    bM <- abs(am$bias[am$rho == max(am$rho)])
    cat(sprintf("  [%s] %-22s |bias|: rho=%.2f -> %.4f  rho=%.2f -> %.4f  (delta=%+.4f)\n",
                rg, mth, min(am$rho), b0, max(am$rho), bM, bM - b0))
  }
}
cat(sprintf("\nSaved:\n - %s\n - %s\n",
            file.path(ARC_RES, "R07_informative_summary.csv"),
            file.path(ARC_RES, "R07_informative_combined.rds")))
