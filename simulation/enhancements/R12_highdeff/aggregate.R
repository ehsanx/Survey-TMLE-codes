# =====================================================================
# R12_highdeff / aggregate.R  —  combine per-task RDS into one summary
#
# Reads the private per-task RDS written by run.R (sim_output/arc_runs/
# R12_highdeff/) and writes ONE summary CSV to results/. Column
# conventions (bias/emp_sd/mean_se/se_ratio/coverage/mcse_cov/
# deff_clust/icc_eif) match aggregate_sim.R so figures/tables stay
# consistent with results/sim_full_summary.csv.
#
# Usage:  Rscript aggregate.R         (uses SIM_OUT or the default dir)
# =====================================================================
RUN_ID <- "R12_highdeff"

CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("SIM_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
files <- list.files(OUT, pattern = "^r12_DesignC_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_SMOKE\\.rds$", files)]   # exclude smoke artifacts by default
if (Sys.getenv("INCLUDE_SMOKE") == "1")
  files <- list.files(OUT, pattern = "^r12_DesignC_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no r12_DesignC_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))
diags <- do.call(rbind, lapply(objs, function(o)
  if (!is.null(o$diagnostics)) cbind(scenario = o$scenario, rung = o$rung, o$diagnostics)))

z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows, list(rows$rung, rows$method), function(d) {
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
ord_r <- c("L1_param", "L3_adaptive", "L4_aggressive")
agg <- agg[order(match(agg$rung, ord_r), match(agg$method, ord_m)), ]

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
cat("\n========== R12_highdeff (Design C) summary ==========\n")
print(agg_print, row.names = FALSE)

# realized-DEFF headline (the number reviewers asked for): per rung, FA arm
dh <- agg[agg$method == "Fully-Aware", c("rung", "deff_clust", "icc_eif")]
cat("\nREALIZED deff_clust (Fully-Aware EIF) per rung:\n")
print(transform(dh, deff_clust = round(deff_clust, 2), icc_eif = round(icc_eif, 3)), row.names = FALSE)

RUN_SUMMARY <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(RUN_SUMMARY)) dir.create(RUN_SUMMARY, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(RUN_SUMMARY, paste0(RUN_ID, "_summary.csv"))
write.csv(agg_print, csv, row.names = FALSE)

DATA_RESULTS_RUN <- file.path(DATA_RESULTS, "arc_runs", RUN_ID)
if (!dir.exists(DATA_RESULTS_RUN)) dir.create(DATA_RESULTS_RUN, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(per_rep = rows, summary = agg, per_rep_diag = diags),
        file.path(DATA_RESULTS_RUN, "r12_combined.rds"))
cat("\nSaved:\n -", csv, "\n -", file.path(DATA_RESULTS_RUN, "r12_combined.rds"), "\n")
