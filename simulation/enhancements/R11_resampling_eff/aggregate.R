# =====================================================================
# aggregate.R  —  combine R11 per-task RDS into the Web Appendix D summary
# Usage:  Rscript simulation/enhancements/R11_resampling_eff/aggregate.R
#         (reads R11_OUT; writes results/R11_resampling_eff_summary.csv)
#
# Produces ONE long summary CSV with two blocks distinguished by `kind`:
#   kind="taylor"  -> standard simulation summary (matches aggregate_sim.R
#                     column conventions: bias,emp_sd,mean_se,se_ratio,
#                     coverage,mcse_cov) for every method incl. IPW-svyglm.
#   kind="jack"    -> resampling check on Fully-Aware & Fully-Aware-CF:
#                     mean_se_taylor, mean_se_jack, jk_lin_ratio (mean per-rep
#                     se_jk/se_taylor; ~1 corroborates Eq 8), plus jackknife
#                     coverage (CI from se_jk, t df_jk).
# A separate CI-width view (efficiency) is included as columns mean_ci_width /
# ci_width_ratio_vs_TMLE in the taylor block.
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R11_OUT",
                  file.path(REPO_ROOT, "sim_output", "arc_runs", "R11_resampling_eff"))
files <- list.files(OUT, pattern = "^R11_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("manifest", files)]
if (!length(files)) stop("no R11_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))

tcrit <- function(df) qt(0.975, pmax(1, df))

# ---- block 1: Taylor (design-linearization) summary, all methods ------------
# reference TMLE arm for the efficiency (CI-width) ratio = Fully-Aware-CF
agg_taylor <- do.call(rbind, by(rows, list(rows$scenario, rows$rung, rows$method),
  function(d) {
    if (is.null(d)) return(NULL)
    Psi <- d$Psi[1]; crit <- tcrit(d$df)
    cov <- mean(abs(d$b - Psi) <= crit * d$se)
    data.frame(
      kind = "taylor",
      scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
      n_reps = nrow(d), Psi = Psi,
      bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
      se_ratio = mean(d$se) / sd(d$b),
      coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
      mean_ci_width = mean(2 * crit * d$se),
      stringsAsFactors = FALSE)
  }))
rownames(agg_taylor) <- NULL
# CI-width ratio vs Fully-Aware-CF within each (scenario,rung): <1 => more efficient
agg_taylor$ci_width_ratio_vs_CF <- NA_real_
for (k in unique(paste(agg_taylor$scenario, agg_taylor$rung))) {
  sel <- paste(agg_taylor$scenario, agg_taylor$rung) == k
  ref <- agg_taylor$mean_ci_width[sel & agg_taylor$method == "Fully-Aware-CF"]
  if (length(ref) == 1L && is.finite(ref))
    agg_taylor$ci_width_ratio_vs_CF[sel] <- agg_taylor$mean_ci_width[sel] / ref
}

# ---- block 2: jackknife resampling check (only arms with se_jk) -------------
jrows <- rows[is.finite(rows$se_jk), , drop = FALSE]
agg_jack <- NULL
if (nrow(jrows)) {
  agg_jack <- do.call(rbind, by(jrows, list(jrows$scenario, jrows$rung, jrows$method),
    function(d) {
      if (is.null(d)) return(NULL)
      Psi <- d$Psi[1]
      crit_t <- tcrit(d$df)                 # Taylor CI
      crit_j <- tcrit(d$df_jk)              # jackknife CI
      cov_t  <- mean(abs(d$b - Psi) <= crit_t * d$se)
      cov_j  <- mean(abs(d$b - Psi) <= crit_j * d$se_jk)
      data.frame(
        kind = "jack",
        scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
        n_reps = nrow(d), Psi = Psi,
        bias = mean(d$b) - Psi, emp_sd = sd(d$b),
        mean_se_taylor = mean(d$se), mean_se_jack = mean(d$se_jk),
        # per-rep ratio of jackknife to linearization SE; ~1 corroborates Eq 8
        jk_lin_ratio = mean(d$se_jk / d$se),
        se_ratio_taylor = mean(d$se) / sd(d$b),
        se_ratio_jack   = mean(d$se_jk) / sd(d$b),
        coverage_taylor = cov_t, coverage_jack = cov_j,
        mcse_cov = sqrt(cov_t * (1 - cov_t) / nrow(d)),
        stringsAsFactors = FALSE)
    }))
  rownames(agg_jack) <- NULL
}

# ---- write a single long CSV (rbind.fill-style: union of columns) -----------
all_cols <- union(names(agg_taylor), if (!is.null(agg_jack)) names(agg_jack) else character(0))
pad <- function(df) { for (c in setdiff(all_cols, names(df))) df[[c]] <- NA; df[all_cols] }
combined <- agg_taylor; if (!is.null(agg_jack)) combined <- rbind(pad(agg_taylor), pad(agg_jack))

ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV",
           "Partially-Aware", "Non-Aware", "IPW-svyglm")
combined <- combined[order(combined$kind, combined$scenario, combined$rung,
                           match(combined$method, ord_m)), ]
num <- sapply(combined, is.numeric); combined[num] <- round(combined[num], 4)

out_dir <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(out_dir, "R11_resampling_eff_summary.csv")
write.csv(combined, csv, row.names = FALSE)
saveRDS(list(per_rep = rows, summary = combined),
        file.path(out_dir, "R11_resampling_eff_combined.rds"))

print(combined, row.names = FALSE)
cat("\nSaved:\n -", csv,
    "\n -", file.path(out_dir, "R11_resampling_eff_combined.rds"), "\n")
