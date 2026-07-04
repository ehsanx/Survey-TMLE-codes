# =====================================================================
# render_fig_ladder_allarms.R  --  Web-Appendix companion to the B&W
# main-text Figure 1 (ladder coverage). This colour version adds the
# Fully-Aware-CV (internal cross-validation) arm that Figure 1 omits,
# present only at L2/L3 (a short segment), so the appendix shows that
# internal CV does NOT rescue coverage the way cross-fitting (CF) does.
#
# Reuses the talk-figure logic (Writing/present/make_fig1_cv.R) but writes
# to the APPENDIX images dir and reads results/sim_full_summary.csv via a
# repo-relative path (Rscript cwd MUST be the repo root).
#
# Reads  results/sim_full_summary.csv
# Writes outputs/figures/fig_ladder_allarms.png  (300 dpi)
#        outputs/figures/fig_ladder_allarms.pdf
# Palette = Okabe-Ito (colourblind-safe), shared with the CSEB talk so the
# slides and the paper agree.
# =====================================================================

summ <- read.csv("results/sim_full_summary.csv", stringsAsFactors = FALSE)

# ---- SHARED arm palette (Okabe-Ito, colourblind-safe) ----------------
ARM_COL <- c("Non-Aware"        = "#E69F00",   # orange
             "Partially-Aware"  = "#009E73",   # bluish green
             "Fully-Aware"      = "#D55E00",   # vermillion (single-fit)
             "Fully-Aware-CV"   = "#CC79A7",   # reddish purple (internal CV)
             "Fully-Aware-CF"   = "#0072B2")   # blue (cross-fitted = proposed)
# shapes: sq, tri, diamond, inverted-tri, circle
ARM_PCH <- c("Non-Aware" = 15, "Partially-Aware" = 17, "Fully-Aware" = 18,
             "Fully-Aware-CV" = 25, "Fully-Aware-CF" = 19)

RUNGS    <- c("L1_param", "L2_smooth", "L3_adaptive", "L4_aggressive")
RUNG_LAB <- c("L1\nGLM", "L2\nsmooth", "L3\n+RF", "L4\ndeep RF")
SCEN     <- c("standard", "R1")
SCEN_LAB <- c("(a) Design A (6 PSUs/stratum)", "(b) Design B (2 PSUs/stratum)")

# draw order (CF last = on top); legend order is pedagogical
ARMS    <- c("Non-Aware", "Partially-Aware", "Fully-Aware", "Fully-Aware-CV", "Fully-Aware-CF")
LEG_ORD <- c("Fully-Aware-CF", "Fully-Aware", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
LEG_LAB <- c("Fully-Aware-CF (cross-fitted)", "Fully-Aware (single-fit)",
             "Fully-Aware-CV (internal CV)", "Partially-Aware", "Non-Aware")
LWD <- c("Non-Aware" = 2.2, "Partially-Aware" = 2.2, "Fully-Aware" = 2.4,
         "Fully-Aware-CV" = 2.8, "Fully-Aware-CF" = 3.6)

cov_at <- function(scen, arm)
  sapply(RUNGS, function(r) {
    v <- summ$coverage[summ$scenario == scen & summ$rung == r & summ$method == arm]
    if (length(v) == 1) v else NA_real_
  })

draw <- function() {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1, 2), mar = c(4.8, 4.8, 2.8, 0.8), mgp = c(2.9, 0.85, 0),
      cex.axis = 1.05, cex.lab = 1.2, font.lab = 2, las = 1)
  for (j in seq_along(SCEN)) {
    scen <- SCEN[j]
    plot(NA, xlim = c(0.9, 4.1), ylim = c(0, 1), xaxt = "n", yaxt = "n",
         xlab = "Super Learner library complexity",
         ylab = "Empirical 95% CI coverage", main = SCEN_LAB[j], cex.main = 1.2)
    axis(1, at = 1:4, labels = RUNG_LAB, padj = 0.5)
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 0.95))   # 0.95 labelled ON the axis
    abline(h = 0.95, lty = 1, col = "grey60", lwd = 1.2)
    for (arm in ARMS) {
      y <- cov_at(scen, arm); x <- 1:4; ok <- !is.na(y)
      lines(x[ok], y[ok], lty = 1, lwd = LWD[[arm]], col = ARM_COL[[arm]])
      points(x[ok], y[ok], pch = ARM_PCH[[arm]], cex = 1.5, lwd = 1.8,
             col = ARM_COL[[arm]], bg = ARM_COL[[arm]])
    }
    if (j == 1)
      legend("bottomleft", legend = LEG_LAB, col = ARM_COL[LEG_ORD],
             pch = ARM_PCH[LEG_ORD], pt.bg = ARM_COL[LEG_ORD],
             lwd = LWD[LEG_ORD], lty = 1, bty = "n", cex = 0.92,
             pt.cex = 1.35, seg.len = 2.4, y.intersp = 1.18)
  }
}

img <- file.path("outputs", "figures")
if (!dir.exists(img)) dir.create(img, recursive = TRUE)

png(file.path(img, "fig_ladder_allarms.png"),
    width = 8, height = 4.2, units = "in", res = 300)
draw(); dev.off()

pdf(file.path(img, "fig_ladder_allarms.pdf"), width = 8, height = 4.2)
draw(); dev.off()

# ---- sanity-check report: CV at L2/L3 + CF at L4, both designs ----
cv <- summ[summ$method == "Fully-Aware-CV",
           c("scenario", "rung", "coverage")]
cf4 <- summ[summ$method == "Fully-Aware-CF" & summ$rung == "L4_aggressive",
            c("scenario", "rung", "coverage")]
cat("wrote outputs/figures/fig_ladder_allarms.{png,pdf}\n")
cat("--- Fully-Aware-CV coverage (L2/L3 only) ---\n")
print(cv[order(cv$scenario, cv$rung), ], row.names = FALSE)
cat("--- Fully-Aware-CF coverage at L4 ---\n")
print(cf4[order(cf4$scenario), ], row.names = FALSE)
