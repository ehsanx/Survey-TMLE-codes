# =====================================================================
# aggregate.R  —  R20_cv_vs_cf_nondonsker
# Pools the per-task RDS into a per-(scenario, method) summary at the
# non-Donsker multi-learner rung L5, with the headline decision view:
# Fully-Aware-CV vs Fully-Aware-CF vs Fully-Aware (single-fit) coverage
# and se_ratio. Uses aggregate_sim.R's EXACT formulas.
#
# Env: SIM_CODE, REPO_ROOT, DATA_ROOT, R20_OUT (per-task RDS dir),
#      RESULTS_DIR (summary CSV dir; via config.R = REPO_ROOT/results).
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) CODE <- "R"
source(file.path(CODE, "config.R"))

RUN_ID <- "R20_cv_vs_cf_nondonsker"
OUT  <- Sys.getenv("R20_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
ARCR <- file.path(RESULTS_DIR, "arc"); if (!dir.exists(ARCR)) dir.create(ARCR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^r20_.*\\.rds$", full.names = TRUE)  # skips SMOKE_, manifest/
if (!length(files)) stop("no r20_*.rds in ", OUT, " -- run the array first")

per_rep <- do.call(rbind, lapply(files, function(f) {
  o <- readRDS(f); cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)
}))

# exact aggregate_sim.R math
summ <- do.call(rbind, by(per_rep, list(per_rep$scenario, per_rep$method), function(d) {
  if (!nrow(d)) return(NULL)
  Psi <- d$Psi[1]; df <- pmax(1, d$df); crit <- qt(0.975, df)
  cov <- mean(abs(d$b - Psi) <= crit * d$se); n <- nrow(d)
  data.frame(scenario = d$scenario[1], rung = d$rung[1], method = d$method[1], n_reps = n,
             Psi = round(Psi, 4), bias = round(mean(d$b) - Psi, 4),
             emp_sd = round(sd(d$b), 4), mean_se = round(mean(d$se), 4),
             se_ratio = round(mean(d$se) / sd(d$b), 4), coverage = round(cov, 4),
             mcse_cov = round(sqrt(cov * (1 - cov) / n), 4),
             deff_clust = round(mean(d$deff), 4), icc_eif = round(mean(d$icc_eif), 4),
             stringsAsFactors = FALSE)
}))
rownames(summ) <- NULL
ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
summ <- summ[order(summ$scenario, match(summ$method, ord_m)), ]

csv <- file.path(ARCR, paste0(RUN_ID, "_summary.csv"))
write.csv(summ, csv, row.names = FALSE)
saveRDS(list(summary = summ, per_rep = per_rep), file.path(OUT, paste0(RUN_ID, "_combined.rds")))
cat("wrote", csv, "\n\n=== FULL SUMMARY (L5 = glm + earth + deep RF) ===\n")
print(summ, row.names = FALSE)

# ---- HEADLINE: CV vs CF vs single-fit at the non-Donsker library --------
cat("\n=== DECISION VIEW: internal CV vs cross-fitting at a NON-Donsker library ===\n")
cat("Question: at L5 (a multi-learner library that INCLUDES the deep RF),\n",
    "does internal CV's down-weighting of the overfitter recover coverage,\n",
    "or does only PSU-level cross-fitting (CF) hold the nominal level?\n\n", sep = "")
key <- c("Fully-Aware-CF", "Fully-Aware-CV", "Fully-Aware")
lab <- c("Fully-Aware-CF" = "cross-fitted (CF)",
         "Fully-Aware-CV" = "internal CV     ",
         "Fully-Aware"     = "single-fit (FA) ")
for (sc in unique(summ$scenario)) {
  cat(sprintf("-- %s --\n", sc))
  for (m in key) {
    r <- summ[summ$scenario == sc & summ$method == m, ]
    if (nrow(r)) cat(sprintf("   %s  coverage=%.3f  se_ratio=%.3f  bias=%+.4f\n",
                             lab[[m]], r$coverage, r$se_ratio, r$bias))
  }
  cf <- summ[summ$scenario == sc & summ$method == "Fully-Aware-CF", "coverage"]
  cv <- summ[summ$scenario == sc & summ$method == "Fully-Aware-CV", "coverage"]
  if (length(cf) && length(cv)) {
    verdict <- if (cf >= 0.93 && cv < 0.90)
      "CF holds, CV undercovers -> internal CV is NOT a substitute for cross-fitting (expected)."
    else if (cf >= 0.93 && cv >= 0.93)
      "BOTH hold -> internal CV's down-weighting recovered coverage here; report as a nuance (CF still guaranteed; CV depends on the ensemble retreating to smooth members)."
    else "unexpected pattern -> inspect bias/se_ratio jointly before quoting."
    cat(sprintf("   => %s\n", verdict))
  }
}

# ---- optional ladder context from the locked headline run (graceful skip) ----
locked <- file.path(RESULTS_DIR, "sim_full_summary.csv")
if (file.exists(locked)) {
  L <- read.csv(locked, stringsAsFactors = FALSE)
  ctx <- L[L$method %in% c("Fully-Aware-CV", "Fully-Aware-CF") &
           L$rung %in% c("L3_adaptive", "L4_aggressive"),
           c("scenario", "rung", "method", "coverage", "se_ratio")]
  cat("\n=== LADDER CONTEXT (locked run): CV at L3, CF at L4 ===\n")
  print(ctx[order(ctx$scenario, ctx$rung, ctx$method), ], row.names = FALSE)
} else {
  cat("\n[note] locked results/sim_full_summary.csv not found under RESULTS_DIR; ladder-context view skipped.\n")
}
