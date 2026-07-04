# =====================================================================
# nhanes/05_figures.R  —  clean forest plots of the five arms per example
#
# Standard forest layout: arm labels at left, point + 95% CI in the data region,
# the numeric estimate in a RIGHT-HAND GUTTER (no overlap with the lines), a dashed
# reference at 0, and the primary Fully-Aware-CF drawn as a filled marker.
# Reads nhanes/results/<id>/estimators_all_arms.rds; writes the report figure
# nhanes/figures/forest_all.{png,pdf} AND the manuscript Figure 2
# outputs/figures/fig2_nhanes_forest.{png,pdf}.
# =====================================================================

REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
res_dir <- file.path(REPO_ROOT, "nhanes", "results")
fig_dir <- file.path(REPO_ROOT, "nhanes", "figures"); dir.create(fig_dir, showWarnings = FALSE)
ms_dir  <- file.path(REPO_ROOT, "outputs", "figures"); dir.create(ms_dir, recursive = TRUE, showWarnings = FALSE)
ids <- c("E1", "E2", "E3", "E4")

ARMS <- c("Non-Aware", "Partially-Aware", "Fully-Aware", "Fully-Aware-CV", "Fully-Aware-CF")
PCH  <- c(1, 2, 0, 5, 19)        # CF = filled circle (primary)
# match render_table2_tex.R's rounding (round-then-format) so Figure 2 endpoints == Table 3
f3 <- function(x) sprintf("%.3f", round(x, 3))

# Figure 2 plots the SAME m=40 multiple-imputation estimates as main-text Table 2
# (single source of truth: results/R06_mi_summary.csv), so the figure and the
# table are identical by construction. The per-example RDS is read only for the arm
# labels/structure; its point estimates and CIs are overwritten from the MI summary.
.mi_fig <- read.csv(file.path(REPO_ROOT, "results", "R06_mi_summary.csv"), stringsAsFactors = FALSE)
load_res <- function(id) {
  f <- file.path(res_dir, id, "estimators_all_arms.rds"); if (!file.exists(f)) return(NULL)
  o <- readRDS(f); m <- .mi_fig[.mi_fig$example == id, ]
  for (i in seq_len(nrow(o$results))) {
    j <- which(m$method == o$results$method[i])
    if (length(j) == 1) { o$results$b[i] <- m$b[j]; o$results$lcl[i] <- m$lcl[j]; o$results$ucl[i] <- m$ucl[j] }
  }
  o
}

draw_forest <- function(o, main = NULL) {
  rr <- o$results[match(ARMS, o$results$method), ]; rr <- rr[!is.na(rr$method), ]
  k <- nrow(rr); yp <- rev(seq_len(k))
  dr  <- range(c(rr$lcl, rr$ucl, 0), na.rm = TRUE); w <- diff(dr); if (w == 0) w <- 1
  pad <- 0.10 * w
  txt_x <- dr[2] + 0.18 * w                         # estimate column starts here (right of all CIs)
  xlim  <- c(dr[1] - pad, dr[2] + 1.45 * w)         # reserve a right gutter for the text
  par(mar = c(4, 8.2, 2.6, 0.6))
  plot(NA, xlim = xlim, ylim = c(0.4, k + 0.6), yaxt = "n", xaxt = "n", bty = "n",
       xlab = "Risk difference (95% CI)", ylab = "", main = main %||% o$results$label[1], cex.main = 0.96)
  ax <- pretty(dr); ax <- ax[ax >= dr[1] - pad & ax <= dr[2] + pad]
  axis(1, at = ax, cex.axis = 0.85)
  abline(v = 0, lty = 2, col = "grey55")
  axis(2, at = yp, labels = rr$method, las = 1, cex.axis = 0.82, tick = FALSE, line = -0.5)
  for (i in seq_len(k)) {
    pchi <- PCH[match(rr$method[i], ARMS)]
    segments(rr$lcl[i], yp[i], rr$ucl[i], yp[i], lwd = 2)
    points(rr$b[i], yp[i], pch = pchi, cex = 1.25, lwd = 1.7, bg = "white")
    text(txt_x, yp[i], sprintf("%s (%s, %s)", f3(rr$b[i]), f3(rr$lcl[i]), f3(rr$ucl[i])),
         pos = 4, cex = 0.66, col = "grey20", offset = 0)
  }
}
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

res <- setNames(lapply(ids, load_res), ids)
have <- ids[!vapply(res, is.null, logical(1))]
if (!length(have)) stop("no results in ", res_dir, " -- run 03 first")

draw_panels <- function() {
  op <- par(mfrow = c(2, 2)); on.exit(par(op))
  for (id in have) {
    lab <- sub("\\s*\\(hardened\\)", "", res[[id]]$results$label[1])   # drop "(hardened)" from the panel title
    draw_forest(res[[id]], main = paste0(id, ":  ", lab))
  }
}

## report figure
png(file.path(fig_dir, "forest_all.png"), 3200, 2200, res = 300); draw_panels(); dev.off()
pdf(file.path(fig_dir, "forest_all.pdf"), 12, 8.5); draw_panels(); dev.off()
## manuscript Figure 2 (same figure, into the manuscript image dir)
png(file.path(ms_dir, "fig2_nhanes_forest.png"), 3200, 2200, res = 300); draw_panels(); dev.off()
pdf(file.path(ms_dir, "fig2_nhanes_forest.pdf"), 12, 8.5); draw_panels(); dev.off()

## per-example single-panel forests (report)
for (id in have) {
  dir.create(file.path(fig_dir, id), showWarnings = FALSE)
  png(file.path(fig_dir, id, "forest.png"), 2200, 1200, res = 300); draw_forest(res[[id]]); dev.off()
}
cat("Forest figures written:\n -", file.path(fig_dir, "forest_all.png"),
    "\n -", file.path(ms_dir, "fig2_nhanes_forest.png"), "\n")
