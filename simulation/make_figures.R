# =====================================================================
# make_figures.R  —  Figure 1 (ladder coverage curve) + Table 1 (compact
#   simulation summary) for the Biometric Methodology manuscript.
#
#   Reads the LOCKED full-run results from sim_full_combined.rds ($summary;
#   falls back to results/sim_full_summary.csv if the RDS is unavailable).
#
#   Outputs:
#     outputs/figures/fig1_ladder_coverage.png   (300 dpi, B&W)
#     outputs/figures/fig1_ladder_coverage.pdf   (vector, for LaTeX)
#     results/table1_main.csv                     (the table data)
#     outputs/tables/table1_main.tex              (booktabs fragment)
#
#   Usage:  Rscript simulation/make_figures.R   (run from the repo root;
#           config.R is located robustly, as in aggregate.R)
# =====================================================================

# ---- locate config.R robustly (mirrors aggregate.R) ------------------
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from the repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

# ---- load the LOCKED summary -----------------------------------------
rds_path <- file.path(DATA_RESULTS, "sim_full_combined.rds")
csv_path <- file.path(RESULTS_DIR, "sim_full_summary.csv")
if (file.exists(rds_path)) {
  summ <- readRDS(rds_path)$summary
  cat("read summary from", rds_path, "\n")
} else if (file.exists(csv_path)) {
  summ <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("RDS not found; read summary from", csv_path, "\n")
} else {
  stop("neither ", rds_path, " nor ", csv_path, " exists")
}

# ---- canonical orderings / labels ------------------------------------
RUNGS     <- c("L1_param", "L2_smooth", "L3_adaptive", "L4_aggressive")
RUNG_LAB  <- c("L1\nGLM", "L2\nsmooth", "L3\n+RF", "L4\ndeep RF")
SCEN      <- c("standard", "R1")
SCEN_LAB  <- c("(a) Design A (6 PSUs/stratum)", "(b) Design B (2 PSUs/stratum)")

# the four headline arms for Figure 1 (CV is left to a Web Figure; it only
# runs at L2/L3 and would clutter the headline -- per plan-phase3.md Step 1)
FIG_ARMS  <- c("Fully-Aware", "Fully-Aware-CF", "Partially-Aware", "Non-Aware")
FIG_LEG   <- c("Fully-Aware (single-fit)", "Fully-Aware-CF (cross-fitted)",
               "Partially-Aware", "Non-Aware")
FIG_LTY   <- c(2, 1, 3, 4)        # dashed / solid / dotted / dot-dash
FIG_PCH   <- c(1, 19, 2, 5)       # open circle / filled circle / triangle / diamond
FIG_LWD   <- c(2, 2.4, 2, 2)

# panel (c) data: R20 runs all five arms at a deployable multi-learner non-Donsker
# library (glm + earth + deep RF), where the internal-CV arm is defined. It shows
# that internal CV under-covers like single-fit while only cross-fitting holds, at
# a realistic ensemble rather than the lone-forest L4. Same B&W lty/pch coding as
# panels (a)/(b); coverage is read from the R20 summary CSV.
r20 <- tryCatch(read.csv(file.path(RESULTS_DIR, "R20_cv_vs_cf_nondonsker_summary.csv"),
                         stringsAsFactors = FALSE), error = function(e) NULL)
C_ARMS <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV")
C_LEG  <- c("Fully-Aware (single-fit)", "Fully-Aware-CF (cross-fitted)", "Fully-Aware-CV")
C_LTY  <- c(2, 1, 3); C_PCH <- c(1, 19, 6); C_LWD <- c(2, 2.4, 2)
covc <- function(scen, arm) {
  if (is.null(r20)) return(NA_real_)
  v <- r20$coverage[r20$scenario == scen & r20$method == arm]
  if (length(v) == 1) v else NA_real_
}

cov_at <- function(scen, arm) {
  # coverage vector along the ladder L1..L4 for one scenario/arm (NA if absent)
  sapply(RUNGS, function(r) {
    v <- summ$coverage[summ$scenario == scen & summ$rung == r & summ$method == arm]
    if (length(v) == 1) v else NA_real_
  })
}

# ---- the plot routine (one call -> one device) -----------------------
draw_fig1 <- function() {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.82))
  on.exit(layout(1), add = TRUE)
  par(mar = c(4.4, 4.2, 2.4, 0.8), mgp = c(2.5, 0.8, 0),
      cex.axis = 0.85, cex.lab = 0.92, las = 1)
  for (j in seq_along(SCEN)) {
    scen <- SCEN[j]
    plot(NA, xlim = c(0.9, 4.1), ylim = c(0, 1), xaxt = "n",
         xlab = "Super Learner library complexity",
         ylab = "Empirical 95% CI coverage", main = SCEN_LAB[j], cex.main = 1.0)
    axis(1, at = 1:4, labels = RUNG_LAB, padj = 0.5)
    abline(h = 0.95, lty = 1, col = "grey55", lwd = 1)
    abline(h = c(0, 1), col = "grey85", lwd = 0.6)
    text(4.05, 0.965, "0.95", cex = 0.7, col = "grey35", pos = 2)
    for (k in seq_along(FIG_ARMS)) {
      y <- cov_at(scen, FIG_ARMS[k]); x <- 1:4
      ok <- !is.na(y)
      lines(x[ok], y[ok], lty = FIG_LTY[k], lwd = FIG_LWD[k])
      points(x[ok], y[ok], pch = FIG_PCH[k], cex = 1.1, bg = "white", lwd = 1.5)
    }
    if (j == 1) {
      legend("bottomleft", legend = FIG_LEG, lty = FIG_LTY, pch = FIG_PCH,
             lwd = FIG_LWD, bty = "n", cex = 0.72, pt.cex = 1.0,
             pt.bg = "white", seg.len = 2.6, y.intersp = 1.1)
    }
  }
  # --- panel (c): deployable multi-learner ensemble (non-Donsker), from R20 ---
  if (!is.null(r20)) {
    plot(NA, xlim = c(0.8, 2.2), ylim = c(0, 1), xaxt = "n",
         xlab = "Design",
         ylab = "Empirical 95% CI coverage",
         main = "(c) Multi-learner ensemble\n(non-Donsker)", cex.main = 1.0)
    axis(1, at = 1:2, labels = c("A", "B"), padj = 0.5)
    abline(h = 0.95, lty = 1, col = "grey55", lwd = 1)
    abline(h = c(0, 1), col = "grey85", lwd = 0.6)
    text(2.18, 0.965, "0.95", cex = 0.7, col = "grey35", pos = 2)
    for (k in seq_along(C_ARMS)) {
      y <- c(covc("standard", C_ARMS[k]), covc("R1", C_ARMS[k])); x <- 1:2
      ok <- !is.na(y)
      lines(x[ok], y[ok], lty = C_LTY[k], lwd = C_LWD[k])
      points(x[ok], y[ok], pch = C_PCH[k], cex = 1.1, bg = "white", lwd = 1.5)
    }
    legend("bottomleft", legend = C_LEG, lty = C_LTY, pch = C_PCH,
           lwd = C_LWD, bty = "n", cex = 0.62, pt.cex = 1.0,
           pt.bg = "white", seg.len = 2.4, y.intersp = 1.1)
  }
}

img_dir <- file.path(REPO_ROOT, "outputs", "figures")
if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE)

png(file.path(img_dir, "fig1_ladder_coverage.png"),
    width = 10.5, height = 4.0, units = "in", res = 300)
draw_fig1(); dev.off()

pdf(file.path(img_dir, "fig1_ladder_coverage.pdf"),
    width = 10.5, height = 4.0)
draw_fig1(); dev.off()

cat("wrote Figure 1 ->", file.path(img_dir, "fig1_ladder_coverage.{png,pdf}"), "\n")

# =====================================================================
# Table 1 — compact summary: rungs {L1, L3, L4} x {standard, R1} x arms.
#   CV arm shown ONLY at L3 (the rung where it runs and where CV != CF is
#   sharpest), keeping it in the main paper (author decision) without
#   making the table ragged elsewhere.
#   Columns per scenario: bias, empirical SD, mean SE, coverage.
# =====================================================================
T_RUNGS <- c("L1_param", "L3_adaptive", "L4_aggressive")
T_RUNG_LAB <- c("L1 (GLM)", "L3 (+RF)", "L4 (deep RF)")
arms_for <- function(rung) {
  base <- c("Fully-Aware", "Fully-Aware-CF", "Partially-Aware", "Non-Aware")
  if (rung == "L3_adaptive")
    c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
  else base
}
ARM_LAB <- c("Fully-Aware"      = "Fully-Aware (single-fit)",
             "Fully-Aware-CF"   = "Fully-Aware-CF",
             "Fully-Aware-CV"   = "Fully-Aware-CV",
             "Partially-Aware"  = "Partially-Aware",
             "Non-Aware"        = "Non-Aware")

get1 <- function(scen, rung, arm, col)
  summ[[col]][summ$scenario == scen & summ$rung == rung & summ$method == arm]

# SE/SD ratio and coverage MCSE: prefer the locked columns, else derive them
# (se_ratio = mean SE / empirical SD; mcse = sqrt(p(1-p)/n_reps)) so the table
# is correct whether the source RDS/CSV carries the columns or not.
se_ratio_at <- function(scen, rung, arm) {
  v <- get1(scen, rung, arm, "se_ratio")
  if (length(v) == 1 && is.finite(v)) return(v)
  get1(scen, rung, arm, "mean_se") / get1(scen, rung, arm, "emp_sd")
}
mcse_at <- function(scen, rung, arm) {
  v <- get1(scen, rung, arm, "mcse_cov")
  if (length(v) == 1 && is.finite(v)) return(v)
  p <- get1(scen, rung, arm, "coverage"); n <- get1(scen, rung, arm, "n_reps")
  sqrt(p * (1 - p) / n)
}

tbl <- do.call(rbind, lapply(seq_along(T_RUNGS), function(i) {
  rung <- T_RUNGS[i]
  do.call(rbind, lapply(arms_for(rung), function(arm) {
    data.frame(
      rung = T_RUNG_LAB[i], arm = ARM_LAB[[arm]],
      bias_std = get1("standard", rung, arm, "bias"),
      sd_std   = get1("standard", rung, arm, "emp_sd"),
      se_std   = get1("standard", rung, arm, "mean_se"),
      ser_std  = se_ratio_at("standard", rung, arm),
      cov_std  = get1("standard", rung, arm, "coverage"),
      mcse_std = mcse_at("standard", rung, arm),
      bias_r1  = get1("R1", rung, arm, "bias"),
      sd_r1    = get1("R1", rung, arm, "emp_sd"),
      se_r1    = get1("R1", rung, arm, "mean_se"),
      ser_r1   = se_ratio_at("R1", rung, arm),
      cov_r1   = get1("R1", rung, arm, "coverage"),
      mcse_r1  = mcse_at("R1", rung, arm),
      stringsAsFactors = FALSE)
  }))
}))
rownames(tbl) <- NULL

# write the data CSV (3 dp for bias/sd/se, 3 dp for coverage)
tbl_csv <- tbl
num <- sapply(tbl_csv, is.numeric)
tbl_csv[num] <- lapply(tbl_csv[num], function(x) round(x, 3))
write.csv(tbl_csv, file.path(RESULTS_DIR, "table1_main.csv"), row.names = FALSE)
cat("wrote", file.path(RESULTS_DIR, "table1_main.csv"), "\n")
print(tbl_csv, row.names = FALSE)

# ---- booktabs LaTeX fragment (horizontal rules only; aligned decimals) --
f3 <- function(x) sprintf("%.3f", x)
f3s <- function(x) { s <- sprintf("%+.3f", x); s }   # signed bias
f2 <- function(x) sprintf("%.2f", x)                 # SE/SD ratio (2 dp)
covcell <- function(p, mc)                           # coverage with MCSE in parens
  sprintf("%.3f\\,{\\scriptsize(%.3f)}", p, mc)
lines <- c(
  "% Auto-generated by simulation/make_figures.R -- do not edit by hand.",
  "\\begin{tabular}{ll rrrrr rrrrr}",
  "\\toprule",
  "\\multicolumn{12}{@{}l}{\\footnotesize Kang--Schafer misspecified DGP (\\texttt{model\\_type=complex}); correctly-specified control in Web Appendix~D} \\\\",
  "\\addlinespace[2pt]",
  " & & \\multicolumn{5}{c}{Design A} & \\multicolumn{5}{c}{Design B} \\\\",
  "\\cmidrule(lr){3-7}\\cmidrule(lr){8-12}",
  "Library & Estimator & Bias & SD & SE & SE/SD & Cov & Bias & SD & SE & SE/SD & Cov \\\\",
  "\\midrule")
cur_rung <- ""
for (r in seq_len(nrow(tbl))) {
  row <- tbl[r, ]
  if (row$rung != cur_rung) {
    if (cur_rung != "") lines <- c(lines, "\\addlinespace")
    cur_rung <- row$rung
  }
  # put the rung label only on the first row of each block
  first_in_block <- (r == 1) || (tbl$rung[r] != tbl$rung[r - 1])
  rcell <- if (first_in_block) row$rung else ""
  lines <- c(lines, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    rcell, row$arm,
    f3s(row$bias_std), f3(row$sd_std), f3(row$se_std), f2(row$ser_std), covcell(row$cov_std, row$mcse_std),
    f3s(row$bias_r1),  f3(row$sd_r1),  f3(row$se_r1),  f2(row$ser_r1),  covcell(row$cov_r1, row$mcse_r1)))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
tex_path <- file.path(REPO_ROOT, "outputs", "tables", "table1_main.tex")
writeLines(lines, tex_path)
cat("wrote", tex_path, "\n")

cat("\nDONE. Figure 1 + Table 1 generated from locked results.\n")
