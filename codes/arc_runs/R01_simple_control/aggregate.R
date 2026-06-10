# =====================================================================
# R01_simple_control / aggregate.R
#   Combine the per-task RDS from run.R into ONE summary CSV that mirrors
#   aggregate_sim.R's column conventions EXACTLY, with one extra leading
#   column: model_type. Writes to results/arc/R01_simple_control_summary.csv
#   (never to the locked results/sim_full_summary.csv).
#
# Coverage / bias / se_ratio / MCSE conventions are copied verbatim from
# aggregate_sim.R so the numbers are directly comparable:
#   crit     = qt(0.975, max(1, df))               (t reference)
#   coverage = mean(|b - Psi| <= crit*se)
#   bias     = mean(b) - Psi ; emp_sd = sd(b)
#   se_ratio = mean(se)/sd(b) ; mcse_cov = sqrt(cov*(1-cov)/n_reps)
#
# Usage:
#   Rscript codes/arc_runs/R01_simple_control/aggregate.R
#   (reads ARC_OUT; writes RESULTS_DIR/arc/R01_simple_control_summary.csv)
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("codes", "codes")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from the repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

RUN_ID  <- "R01_simple_control"
ARC_OUT <- Sys.getenv("ARC_OUT", file.path(dirname(DATA_INTERMEDIATE), "arc_runs", RUN_ID))
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RESULTS)) dir.create(ARC_RESULTS, recursive = TRUE, showWarnings = FALSE)

files <- list.files(ARC_OUT, pattern = "^r01_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no r01_*.rds found in ", ARC_OUT)
cat("aggregating", length(files), "task files from", ARC_OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, model_type = o$model_type, rung = o$rung,
        Psi = o$Psi, o$results)))

z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows,
  list(rows$scenario, rows$model_type, rows$rung, rows$method), function(d) {
    if (is.null(d)) return(NULL)
    Psi <- d$Psi[1]; crit <- z_or_t(d$df)
    cov <- mean(abs(d$b - Psi) <= crit * d$se)
    data.frame(scenario = d$scenario[1], model_type = d$model_type[1], rung = d$rung[1],
               method = d$method[1], n_reps = nrow(d), Psi = Psi,
               bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
               se_ratio = mean(d$se) / sd(d$b),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
               deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif))
  }))
rownames(agg) <- NULL
ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
agg <- agg[order(agg$scenario, agg$model_type, agg$rung,
                 match(agg$method, ord_m)), ]

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
print(agg_print, row.names = FALSE)

# ---- smoke-gate / decision helper: flag any 'simple' L1 row that violates ----
# the control's expectation (coverage in [0.94,0.95] +/- 2 MCSE; |bias| within
# ~2 MCSE-of-the-mean of 0). Printed so the user sees a PASS/FLAG verdict.
ctrl <- agg[agg$model_type == "simple" & agg$rung == "L1_param", ]
if (nrow(ctrl)) {
  cat("\n--- simple/L1 pipeline-control check (Fully-Aware arm) ---\n")
  fa <- ctrl[ctrl$method == "Fully-Aware", ]
  for (i in seq_len(nrow(fa))) {
    r <- fa[i, ]
    bias_mcse <- r$emp_sd / sqrt(r$n_reps)          # MCSE of the mean estimate
    bias_ok   <- abs(r$bias) <= 2 * bias_mcse
    cov_ok    <- (r$coverage >= 0.93) && (r$coverage <= 0.96)  # ~.94-.95 +/- slack
    cat(sprintf("  %-9s FA: bias=%+.4f (2*MCSE=%.4f) [%s]  coverage=%.3f [%s]\n",
                r$scenario, r$bias, 2 * bias_mcse, if (bias_ok) "OK" else "FLAG",
                r$coverage, if (cov_ok) "OK" else "FLAG"))
  }
  cat("  Rule: any FLAG on a 'simple'/L1 Fully-Aware row => STOP-and-report\n",
      "  (real pipeline finding) BEFORE writing Theorem 1/2 prose.\n", sep = "")
}

out_csv <- file.path(ARC_RESULTS, "R01_simple_control_summary.csv")
write.csv(agg_print, out_csv, row.names = FALSE)
saveRDS(list(per_rep = rows, summary = agg),
        file.path(ARC_RESULTS, "R01_simple_control_combined.rds"))
cat("\nSaved:\n -", out_csv,
    "\n -", file.path(ARC_RESULTS, "R01_simple_control_combined.rds"), "\n")
