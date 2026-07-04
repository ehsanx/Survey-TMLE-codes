# Aggregate J3_calib per-task RDS -> output/newruns/J3_calib_summary.csv
# Canonical coverage math (aggregate_sim.R): crit = qt(.975, max(1, df));
# covered = |b - Psi| <= crit * se.
IN  <- Sys.getenv("J_OUT", "sim_output/enhancements/J3_calib")
fls <- list.files(IN, pattern = "^j3_.*_chunk\\d+\\.rds$", full.names = TRUE)
stopifnot(length(fls) > 0)
allr <- do.call(rbind, lapply(fls, function(f) {
  x <- readRDS(f); cbind(rung = x$rung, Psi = x$Psi, x$results)
}))
allr$covered <- with(allr, abs(b - Psi) <= qt(0.975, pmax(1, df)) * se)
agg <- do.call(rbind, lapply(split(allr, list(allr$rung, allr$weight_scheme, allr$method), drop = TRUE),
  function(g) data.frame(
    rung = g$rung[1], weight_scheme = g$weight_scheme[1], method = g$method[1],
    n_reps = length(unique(g$rep)), Psi = g$Psi[1],
    bias = mean(g$b) - g$Psi[1], emp_sd = sd(g$b), mean_se = mean(g$se),
    se_ratio = mean(g$se) / sd(g$b), coverage = mean(g$covered),
    mcse_cov = sqrt(mean(g$covered) * (1 - mean(g$covered)) / length(g$covered)),
    mean_w_cv = mean(g$w_cv))))
agg <- agg[order(agg$rung, agg$method, agg$weight_scheme), ]
out <- "results/J3_calib_summary.csv"
write.csv(agg, out, row.names = FALSE)
cat("wrote", out, "\n"); print(agg, digits = 3)
cat("\nGATE: compare Fully-Aware-CF coverage design vs raked at each rung;\n",
    "PASS if raked is within ~2 MCSE of design (and both near 0.95 at L1/L2).\n")
