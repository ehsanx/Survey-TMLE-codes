# =====================================================================
# R10_fpc/aggregate.R  —  combine per-task RDS into the two final summary CSVs.
# Run AFTER the array job finishes:  Rscript codes/arc_runs/R10_fpc/aggregate.R
# Reads sim_output/arc_runs/R10_fpc/*.rds, writes results/arc/R10_fpc_*.csv.
# Non-destructive: touches nothing the locked pipeline owns.
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("codes", "codes")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))
RUN_DIR <- file.path(REPO_ROOT, "codes", "arc_runs", "R10_fpc")
source(file.path(RUN_DIR, "fpc_helpers.R"))

ID  <- "R10_fpc"
OUT <- file.path(DATA_ROOT, "arc_runs", ID)
arc_csv_dir <- file.path(RESULTS_DIR, "arc")
dir.create(arc_csv_dir, recursive = TRUE, showWarnings = FALSE)

# ---- PART A: FPC vs no-FPC on R1 (L1_param, L2_smooth) -----------------------
fa_files <- list.files(OUT, pattern = "_partA(_SMOKE)?\\.rds$", full.names = TRUE)
if (length(fa_files)) {
  oA  <- readRDS(fa_files[[1]])                  # PART A is a single task
  rows <- oA$per_rep; Psi <- oA$Psi
  summ <- do.call(rbind, lapply(split(rows, list(rows$rung, rows$method, rows$fpc),
                                      drop = TRUE), function(d)
    summarize_arm(d$b, d$se, d$df, Psi, "R1", d$rung[1], d$method[1], d$fpc[1])))
  rownames(summ) <- NULL
  summ <- summ[order(summ$rung, summ$method, summ$fpc), ]
  num <- sapply(summ, is.numeric); summ[num] <- round(summ[num], 4)
  csvA <- file.path(arc_csv_dir, sprintf("%s_partA_summary.csv", ID))
  write.csv(summ, csvA, row.names = FALSE)
  cat("PART A (fpc vs no-fpc, R1):  ->", csvA, "\n")
  print(summ, row.names = FALSE)
} else cat("PART A: no _partA RDS found in", OUT, "\n")

# ---- PART B: first-stage-fraction sweep (canonical no-fpc CF, L1_param) -------
fb_files <- list.files(OUT, pattern = "_partB_J[0-9]+(_SMOKE)?\\.rds$", full.names = TRUE)
if (length(fb_files)) {
  rowsB <- do.call(rbind, lapply(fb_files, function(f) {
    o <- readRDS(f)
    summ <- summarize_arm(o$per_rep$b, o$per_rep$se, o$per_rep$df, o$Psi,
                          "R1", o$rung, "Fully-Aware-CF", fpc = "no")
    summ$J <- o$J; summ$frac <- o$frac; summ$N_pop <- o$params$N_pop
    summ
  }))
  rowsB <- rowsB[order(rowsB$frac), ]            # decreasing fraction f -> 0
  num <- sapply(rowsB, is.numeric); rowsB[num] <- round(rowsB[num], 4)
  csvB <- file.path(arc_csv_dir, sprintf("%s_fraction_sweep.csv", ID))
  write.csv(rowsB, csvB, row.names = FALSE)
  cat("\nPART B (fraction sweep f->0, R1/L1_param CF):  ->", csvB, "\n")
  print(rowsB, row.names = FALSE)
} else cat("PART B: no _partB RDS found in", OUT, "\n")
