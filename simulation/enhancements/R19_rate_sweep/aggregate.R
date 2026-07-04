# =====================================================================
# aggregate.R  —  combine R19_rate_sweep per-task RDS into the rate-TREND table
# Usage:  Rscript codes/arc_runs/R19_rate_sweep/aggregate.R   (after the array)
# Reads:  R19_OUT (default sim_output/arc_runs/R19_rate_sweep)/r19_*.rds
#         (SMOKE_* files and the manifest/ subdir are skipped: the pattern
#          anchors on ^r19_ and list.files is non-recursive)
# Writes: results/arc/R19_rate_sweep_summary.csv
#         <R19_OUT>/R19_rate_sweep_combined.rds
#
# Two views are printed:
#   1. per (base_m, rung): mean_eQ, mean_eg, mean_prod, mean_prod_sqrtm,
#      sd_prod_sqrtm, prod_sqrtm_vs_L1 (vs L1 at the SAME base_m), n_reps,
#      truth_join  — joins the existing R02/R04 tables on (rung, base_m).
#   2. HEADLINE TREND: per rung, mean_prod_sqrtm(base_m)/mean_prod_sqrtm(base_m=6)
#      across the sweep. Flat/declining at L2-L3 = the assumed product rate (C1)
#      plausibly met; GROWING at L4 = rate fails (the expected demonstration);
#      GROWING at L1 is ALSO expected (SL.glm misspecified for both nuisances
#      -> approximation floor; matches R02's locked L1 bias floor + coverage
#      decay). Only L4 flat/declining or L2 clearly growing signals a bug
#      (see NOTES.md "Interpreting the full-run table").
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R19_OUT", file.path(DATA_ROOT, "arc_runs", "R19_rate_sweep"))
files <- list.files(OUT, pattern = "^r19_.*\\.rds$", full.names = TRUE)  # skips SMOKE_*/manifest/
if (!length(files)) stop("no r19_*.rds found in ", OUT, " (run the array first)")
cat("aggregating", length(files), "task files from", OUT, "\n")

per_rep <- do.call(rbind, lapply(files, function(f) readRDS(f)$per_rep))

# ladder order for sorting
rung_lv <- c("L1_param", "L2_smooth", "L3_adaptive", "L4_aggressive")
per_rep$rung <- factor(per_rep$rung, levels = rung_lv)

# ---- view 1: per (base_m, rung) rate table -----------------------------------
agg <- do.call(rbind, by(per_rep, list(per_rep$base_m, per_rep$rung), function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  data.frame(
    scenario = d$scenario[1], rung = as.character(d$rung[1]),
    base_m = d$base_m[1], m_total = d$m_total[1],          # join keys to R02/R04
    n_reps = nrow(d), mean_n = mean(d$n), mean_m_psu = mean(d$m_psu),
    # PRIMARY (u-integrated Q, uA-integrated g) -- the rate object
    mean_eQ = mean(d$eQ_int), mean_eg = mean(d$eg_int),
    mean_prod = mean(d$prod_int),
    mean_prod_sqrtm = mean(d$prod_int_sqrtm),
    sd_prod_sqrtm   = sd(d$prod_int_sqrtm),
    # SECONDARY (realized-u) for completeness, mirrors R04's table
    mean_eQ_real = mean(d$eQ_real), mean_eg_real = mean(d$eg_real),
    mean_prod_real_sqrtm = mean(d$prod_real_sqrtm),
    truth_join = paste(unique(d$truth_join), collapse = ","),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$base_m, match(agg$rung, rung_lv)), ]
rownames(agg) <- NULL

# prod_sqrtm_vs_L1: ratio vs the L1 rung at the SAME base_m (R04's decision col)
l1_tab <- agg[agg$rung == "L1_param", c("base_m", "mean_prod_sqrtm")]
agg$prod_sqrtm_vs_L1 <- agg$mean_prod_sqrtm /
  l1_tab$mean_prod_sqrtm[match(agg$base_m, l1_tab$base_m)]

# trend_vs_m6: per rung, ratio vs that rung's own base_m=6 value (the RATE probe)
m6_tab <- agg[agg$base_m == 6L, c("rung", "mean_prod_sqrtm")]
agg$trend_vs_m6 <- agg$mean_prod_sqrtm /
  m6_tab$mean_prod_sqrtm[match(agg$rung, m6_tab$rung)]

num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 5)
cat("\n==== R19 view 1: per (base_m, rung) sqrt(m)-scaled product-error table ====\n")
print(ap, row.names = FALSE)

# ---- view 2: HEADLINE TREND — per rung, prod_sqrtm(base_m)/prod_sqrtm(m=6) ----
m_lv <- sort(unique(agg$base_m))
trend <- data.frame(rung = rung_lv[rung_lv %in% agg$rung], stringsAsFactors = FALSE)
for (m in m_lv) {
  v <- agg$trend_vs_m6[match(paste(trend$rung, m), paste(agg$rung, agg$base_m))]
  trend[[sprintf("m%02d_ratio", m)]] <- round(v, 4)
}
cat("\n==== R19 view 2 (HEADLINE): mean_prod_sqrtm(base_m) / mean_prod_sqrtm(base_m=6) ====\n")
cat("(columns = base_m 6/12/20/30 -> m_total 60/120/200/300; per-rung trend)\n")
print(trend, row.names = FALSE)
cat("\nINTERPRETATION: flat/DECLINING ratios across m at L2-L3 -> the assumed\n",
    "product rate (C1) is plausibly met at that rung; GROWING ratios at L4 ->\n",
    "the rate FAILS there (expected: the interpolating deep RF). GROWING at L1\n",
    "is ALSO expected (SL.glm misspecified for both nuisances -> approximation\n",
    "floor; matches R02's locked L1 bias floor 0.0183/0.0169/0.0172/0.0188 and\n",
    "Fully-Aware-CF coverage decay 0.949->0.848). Only L4 flat/declining or L2\n",
    "clearly growing contradicts the locked table -> STOP and investigate.\n", sep = "")

# ---- write (NON-DESTRUCTIVE; results/arc + the run's own dir only) -----------
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RESULTS)) dir.create(ARC_RESULTS, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(ARC_RESULTS, "R19_rate_sweep_summary.csv")
write.csv(ap, csv, row.names = FALSE)
saveRDS(list(per_rep = per_rep, summary = agg, trend = trend),
        file.path(OUT, "R19_rate_sweep_combined.rds"))
cat("\nSaved:\n -", csv, "\n -", file.path(OUT, "R19_rate_sweep_combined.rds"), "\n")
