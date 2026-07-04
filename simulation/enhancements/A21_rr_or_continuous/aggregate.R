# =====================================================================
# aggregate.R  —  A21_rr_or_continuous
# Pools the per-chunk RDS into coverage / bias / width summaries for the
# delta-method estimands (RR, OR on binary; RD, ratio-of-means RR on
# continuous) for the Fully-Aware (single-fit) and Fully-Aware-CF arms.
# Writes results/arc/A21_rr_or_continuous_{rr,rd}_summary.csv.
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) CODE <- "R"
source(file.path(CODE, "config.R"))

RUN_ID <- "A21_rr_or_continuous"
OUT  <- Sys.getenv("A21_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
ARCR <- file.path(RESULTS_DIR, "arc"); if (!dir.exists(ARCR)) dir.create(ARCR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^a21_.*\\.rds$", full.names = TRUE)  # skips SMOKE_, manifest/
if (!length(files)) stop("no a21_*.rds in ", OUT, " -- run run.R first")

rr <- do.call(rbind, lapply(files, function(f) { o <- readRDS(f); cbind(outcome = o$outcome_type, rung = o$rung, o$rr) }))
rd <- do.call(rbind, lapply(files, function(f) { o <- readRDS(f); cbind(outcome = o$outcome_type, rung = o$rung, o$rd) }))

summ_rr <- do.call(rbind, by(rr, list(rr$outcome, rr$rung, rr$method, rr$estimand), function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  tru <- d$truth[1]; n <- nrow(d); cov <- mean(d$cover)
  data.frame(outcome = d$outcome[1], rung = d$rung[1], method = d$method[1], estimand = d$estimand[1],
             n_reps = n, truth = round(tru, 4), est_med = round(median(d$est), 4),
             bias_log = round(mean(d$log_est) - log(tru), 4),
             coverage = round(cov, 4), mcse_cov = round(sqrt(cov * (1 - cov) / n), 4),
             mean_width = round(mean(d$width), 4),
             se_ratio = round(mean(d$se_log) / sd(d$log_est), 4), stringsAsFactors = FALSE) }))
rownames(summ_rr) <- NULL

summ_rd <- do.call(rbind, by(rd, list(rd$outcome, rd$rung, rd$method), function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  tru <- d$truth[1]; n <- nrow(d); cov <- mean(d$cover)
  data.frame(outcome = d$outcome[1], rung = d$rung[1], method = d$method[1], estimand = "RD",
             n_reps = n, truth = round(tru, 4), est_med = round(median(d$b), 4),
             bias = round(mean(d$b) - tru, 4),
             coverage = round(cov, 4), mcse_cov = round(sqrt(cov * (1 - cov) / n), 4),
             mean_width = round(mean(d$width), 4),
             se_ratio = round(mean(d$se) / sd(d$b), 4), stringsAsFactors = FALSE) }))
rownames(summ_rd) <- NULL

ord_m <- c("Fully-Aware", "Fully-Aware-CF")
summ_rr <- summ_rr[order(summ_rr$outcome, summ_rr$rung, match(summ_rr$method, ord_m), summ_rr$estimand), ]
summ_rd <- summ_rd[order(summ_rd$outcome, summ_rd$rung, match(summ_rd$method, ord_m)), ]

write.csv(summ_rr, file.path(ARCR, paste0(RUN_ID, "_rr_summary.csv")), row.names = FALSE)
write.csv(summ_rd, file.path(ARCR, paste0(RUN_ID, "_rd_summary.csv")), row.names = FALSE)
saveRDS(list(rr = summ_rr, rd = summ_rd, n_files = length(files)),
        file.path(OUT, paste0(RUN_ID, "_combined.rds")))
cat(sprintf("pooled %d chunk files\n\n=== RR / OR (delta method) ===\n", length(files)))
print(summ_rr, row.names = FALSE)
cat("\n=== RD (context) ===\n"); print(summ_rd, row.names = FALSE)
cat(sprintf("\nwrote %s/%s_{rr,rd}_summary.csv\n", ARCR, RUN_ID))
