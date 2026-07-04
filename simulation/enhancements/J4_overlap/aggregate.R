# Aggregate J4_overlap per-task RDS -> output/newruns/J4_overlap_summary.csv
IN  <- Sys.getenv("J_OUT", "sim_output/enhancements/J4_overlap")
fls <- list.files(IN, pattern = "^j4_.*_chunk\\d+\\.rds$", full.names = TRUE)
stopifnot(length(fls) > 0)
allr <- do.call(rbind, lapply(fls, function(f) {
  x <- readRDS(f)
  cbind(rung = x$rung, alpha_tag = x$alpha_tag, alpha_g = x$alpha_g,
        Psi = x$Psi, x$results)
}))
allr$covered <- with(allr, abs(b - Psi) <= qt(0.975, pmax(1, df)) * se)
agg <- do.call(rbind, lapply(split(allr, list(allr$rung, allr$alpha_tag, allr$method), drop = TRUE),
  function(g) data.frame(
    rung = g$rung[1], alpha_tag = g$alpha_tag[1], alpha_g = g$alpha_g[1],
    method = g$method[1], n_reps = length(unique(g$rep)), Psi = g$Psi[1],
    bias = mean(g$b) - g$Psi[1], emp_sd = sd(g$b), mean_se = mean(g$se),
    se_ratio = mean(g$se) / sd(g$b), coverage = mean(g$covered),
    mcse_cov = sqrt(mean(g$covered) * (1 - mean(g$covered)) / length(g$covered)))))
agg <- agg[order(agg$rung, agg$alpha_g, agg$method), ]
out <- "results/J4_overlap_summary.csv"
write.csv(agg, out, row.names = FALSE)
cat("wrote", out, "\n"); print(agg, digits = 3)
cat("\nNOTE: Psi should be ~invariant across alpha_tag (the knob only enters\n",
    "treatment assignment); a drifting Psi means an aggregation bug, stop.\n",
    "GATE: report CF coverage vs alpha_g; degradation at 'high' is a FINDING to\n",
    "report honestly (bounded-relative-weights/positivity edge), not a bug.\n")
