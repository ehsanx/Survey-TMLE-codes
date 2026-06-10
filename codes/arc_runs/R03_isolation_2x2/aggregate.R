# =====================================================================
# aggregate.R  —  combine R03_isolation_2x2 per-task RDS into a summary CSV
# Usage:  Rscript aggregate.R     (after the full array job completes)
# Mirrors codes/aggregate_sim.R conventions (same coverage/se_ratio/MCSE math
# and the same column names) so the result drops straight into the figures.
# Writes results/arc/R03_isolation_2x2_summary.csv -- it does NOT clobber the
# locked results/sim_full_summary.csv.
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("codes", "codes")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R03_OUT",
                  file.path(REPO_ROOT, "sim_output", "arc_runs", "R03_isolation_2x2"))
# exclude SMOKE files from the full aggregation
files <- list.files(OUT, pattern = "^R03_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_SMOKE", files)]
if (!length(files)) stop("no full R03_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, g_floor = o$g_floor, Psi = o$Psi, o$results)))

z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows, list(rows$scenario, rows$rung, rows$method), function(d) {
  if (is.null(d)) return(NULL)
  Psi <- d$Psi[1]
  ok <- is.finite(d$b) & is.finite(d$se)        # drop non-convergent (diverged) reps
  nd <- sum(!ok); nrok <- sum(ok)
  bo <- d$b[ok]; so <- d$se[ok]; crit <- z_or_t(d$df[ok])
  cov <- mean(abs(bo - Psi) <= crit * so)
  data.frame(scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
             n_reps = nrok, n_diverged = nd, Psi = Psi,
             bias = mean(bo) - Psi, emp_sd = sd(bo), mean_se = mean(so),
             se_ratio = mean(so) / sd(bo),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / max(1, nrok)),
             deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif))
}))
rownames(agg) <- NULL
ord_m <- c("SF-W", "SF-U", "CF-W", "CF-U")
ord_r <- c("L4_aggressive", "L1_param")
agg <- agg[order(match(agg$rung, ord_r), match(agg$method, ord_m)), ]

num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
print(ap, row.names = FALSE)

dir.create(file.path(RESULTS_DIR, "arc"), recursive = TRUE, showWarnings = FALSE)
csv <- file.path(RESULTS_DIR, "arc", "R03_isolation_2x2_summary.csv")
write.csv(ap, csv, row.names = FALSE)
saveRDS(list(per_rep = rows, summary = agg),
        file.path(RESULTS_DIR, "arc", "R03_isolation_2x2_combined.rds"))
cat("\nSaved:\n -", csv, "\n -",
    file.path(RESULTS_DIR, "arc", "R03_isolation_2x2_combined.rds"), "\n")

# ---- restate the decision rule on the FULL run -------------------------------
sfu <- ap$coverage[ap$rung == "L4_aggressive" & ap$method == "SF-U"]
cf  <- ap$coverage[ap$rung == "L4_aggressive" & ap$method %in% c("CF-W","CF-U")]
cat("\n==== DECISION (L4_aggressive) ====\n")
cat(sprintf("  SF-U (single-fit/UNWEIGHTED) coverage = %.3f\n", sfu))
cat(sprintf("  CF-W/CF-U (cross-fit) coverage = %s\n", paste(round(cf,3), collapse = ", ")))
if (length(sfu) && sfu >= 0.90) {
  cat("  >>> STOP-AND-REPORT: de-weighting, not cross-fitting, may rescue coverage.\n")
} else {
  cat("  >>> CONFIRMS: single-fit under-covers, cross-fit recovers -> cross-fitting is the active ingredient.\n")
}
