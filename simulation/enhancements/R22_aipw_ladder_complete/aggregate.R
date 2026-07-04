# =====================================================================
# aggregate.R  —  R22_aipw_ladder_complete (minimal AIPW completion)
# Pools the per-task RDS into a per-(scenario, rung, method) summary of the
# genuinely-new AIPW cells (AIPW-CV @ L2/L3 ; AIPW-{SF,CV,CF} @ L5), and prints
# them next to their LOCKED neighbours (TMLE-CV @ L2/L3 from sim_full; AIPW @ L4
# from R15; AIPW @ L6 from R21) for a plausibility eyeball. No within-run
# cross-check is possible because R22 recomputes NO locked cell -- trust comes
# from CODE IDENTITY with the already-validated R15/R21 AIPW machinery.
#
# Uses aggregate_sim.R's EXACT formulas (qt(.975, df)).
# Env: SIM_CODE, REPO_ROOT, DATA_ROOT, R22_OUT, RESULTS_DIR.
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) CODE <- "R"
source(file.path(CODE, "config.R"))

RUN_ID <- "R22_aipw_ladder_complete"
OUT  <- Sys.getenv("R22_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
ARCR <- file.path(RESULTS_DIR, "arc"); if (!dir.exists(ARCR)) dir.create(ARCR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^r22_.*\\.rds$", full.names = TRUE)  # skips SMOKE_, manifest/, *_combined
files <- files[!grepl("_combined\\.rds$", files)]
if (!length(files)) stop("no r22_*.rds in ", OUT, " -- run the array first")

per_rep <- do.call(rbind, lapply(files, function(f) {
  o <- readRDS(f)
  cbind(scenario = o$scenario, Psi = o$Psi, o$results)   # carries rung + base_m
}))

# ---- integrity guard: no duplicated (scenario, rung, method, rep) ------------
key <- per_rep[, c("scenario", "rung", "method", "rep")]
stopifnot(!anyDuplicated(key))

summ <- do.call(rbind, by(per_rep, list(per_rep$scenario, per_rep$rung, per_rep$method),
                          function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  Psi <- d$Psi[1]; df <- pmax(1, d$df); crit <- qt(0.975, df)
  cov <- mean(abs(d$b - Psi) <= crit * d$se); n <- nrow(d)
  data.frame(scenario = d$scenario[1], rung = d$rung[1],
             estimator = d$estimator[1], method = d$method[1],
             n_reps = n, Psi = round(Psi, 4), bias = round(mean(d$b) - Psi, 4),
             emp_sd = round(sd(d$b), 4), mean_se = round(mean(d$se), 4),
             se_ratio = round(mean(d$se) / sd(d$b), 4), coverage = round(cov, 4),
             mcse_cov = round(sqrt(cov * (1 - cov) / n), 4),
             deff_clust = round(mean(d$deff), 4), icc_eif = round(mean(d$icc_eif), 4),
             stringsAsFactors = FALSE)
}))
rownames(summ) <- NULL
ord_m <- c("AIPW-SF", "AIPW-CV", "AIPW-CF")
ord_r <- c("L2_smooth", "L3_adaptive", "L5_nondonsker")
summ <- summ[order(match(summ$rung, ord_r), summ$scenario, match(summ$method, ord_m)), ]

csv <- file.path(ARCR, paste0(RUN_ID, "_summary.csv"))
write.csv(summ, csv, row.names = FALSE)
saveRDS(list(summary = summ, per_rep = per_rep), file.path(OUT, paste0(RUN_ID, "_combined.rds")))
cat("wrote", csv, "\n\n=== DELIVERABLE cells (the new AIPW numbers that feed the appendix figure) ===\n")
print(summ[, c("rung","scenario","method","coverage","se_ratio","mcse_cov","bias","n_reps")], row.names = FALSE)

# ---- plausibility context: print LOCKED neighbours (eyeball, not a check) ----
cat("\n=== LOCKED-NEIGHBOUR CONTEXT (eyeball plausibility; NOT a cross-check) ===\n")
say <- function(lbl, path, filt) {
  if (!file.exists(path)) { cat(sprintf("  [skip] %s (%s not found)\n", lbl, basename(path))); return(invisible()) }
  L <- read.csv(path, stringsAsFactors = FALSE); r <- filt(L)
  if (nrow(r)) { cat(sprintf("  -- %s --\n", lbl)); print(r, row.names = FALSE) }
}
SF <- file.path(RESULTS_DIR, "sim_full_summary.csv")
say("TMLE-CV @ L2/L3 (sim_full) -- compare to AIPW-CV @ L2/L3", SF, function(L)
  L[L$method=="Fully-Aware-CV" & L$rung %in% c("L2_smooth","L3_adaptive"),
    c("scenario","rung","method","coverage","se_ratio")])
R15 <- file.path(ARCR, "R15_aipw_benchmark_summary.csv")
say("AIPW @ L4 (R15, lin) -- the collapse below L5", R15, function(L)
  L[L$se_type=="lin" & L$rung=="L4_aggressive" & L$method %in% c("AIPW-SF","AIPW-CF"),
    c("scenario","rung","method","coverage","se_ratio")])
R21 <- file.path(ARCR, "R21_deployable_cvcf_summary.csv")
say("AIPW @ L6 (R21) -- the deployable ensemble above L5", R21, function(L)
  L[L$rung=="L6_deployable" & L$method %in% c("AIPW-SF","AIPW-CV","AIPW-CF"),
    c("scenario","rung","method","coverage","se_ratio")])

cat("\nSaved:\n -", csv, "\n -", file.path(OUT, paste0(RUN_ID, "_combined.rds")), "\n")
cat("\nNOTE: AIPW-CV @ L2/L3 should sit near nominal (~0.90-0.94, se_ratio<1) like TMLE-CV there;\n")
cat("      AIPW @ L5 should sit BETWEEN the L4 collapse and the L6 deployable values.\n")
