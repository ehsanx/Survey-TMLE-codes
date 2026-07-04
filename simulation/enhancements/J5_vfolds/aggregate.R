# Aggregate J5_vfolds per-task RDS -> output/newruns/J5_vfolds_summary.csv
IN  <- Sys.getenv("J_OUT", "sim_output/enhancements/J5_vfolds")
fls <- list.files(IN, pattern = "^j5_.*_chunk\\d+\\.rds$", full.names = TRUE)
stopifnot(length(fls) > 0)
allr <- do.call(rbind, lapply(fls, function(f) {
  x <- readRDS(f); cbind(rung = x$rung, Psi = x$Psi, x$results)
}))
allr$covered <- with(allr, abs(b - Psi) <= qt(0.975, pmax(1, df)) * se)
agg <- do.call(rbind, lapply(split(allr, list(allr$rung, allr$V_cf, allr$method), drop = TRUE),
  function(g) data.frame(
    rung = g$rung[1], V_cf = g$V_cf[1], method = g$method[1],
    n_reps = length(unique(g$rep)), Psi = g$Psi[1],
    bias = mean(g$b) - g$Psi[1], emp_sd = sd(g$b), mean_se = mean(g$se),
    se_ratio = mean(g$se) / sd(g$b), coverage = mean(g$covered),
    mcse_cov = sqrt(mean(g$covered) * (1 - mean(g$covered)) / length(g$covered)))))
agg <- agg[order(agg$rung, agg$V_cf, agg$method), ]
out <- "results/J5_vfolds_summary.csv"
write.csv(agg, out, row.names = FALSE)
cat("wrote", out, "\n"); print(agg, digits = 3)
cat("\nANCHOR: the V=5 Fully-Aware-CF rows should reproduce the (new-engine)\n",
    "headline CF cells within MC error. GATE: CF coverage stable across V.\n",
    "Only the Fully-Aware-CF rows vary with V; other arms are V-invariant\n",
    "replicates (a consistency check, not new information).\n")
