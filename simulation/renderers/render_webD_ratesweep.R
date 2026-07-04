# =====================================================================
# render_webD_ratesweep.R  —  Web Appendix D PSU-rate-sweep artifacts.
# Emits BOTH:
#   (1) outputs/tables/webD_ratesweep.tex   (tab:webD_ratesweep)
#   (2) outputs/figures/fig_ratesweep.{png,pdf}
# from results/R19_rate_sweep_summary.csv (200-rep ARC sweep, Design A).
# Numbers are read from the CSV; no hand transcription.
# =====================================================================
source("simulation/renderers/render_helpers_color.R")

src <- "results/R19_rate_sweep_summary.csv"
d <- read.csv(src, check.names = FALSE)
d <- d[d$scenario == "standard", ]

rung_lab <- c(L1_param = "L1 (GLM)", L2_smooth = "L2 (+GAM, MARS)",
              L3_adaptive = "L3 (+RF)", L4_aggressive = "L4 (deep RF)")
rung_ord <- names(rung_lab)
m_ord <- c(60, 120, 200, 300)

# ---------------------------------------------------------------------
# (1) Table: rows = rungs; cols = mean_prod_sqrtm at m=60/120/200/300
#     plus the trend ratio (trend_vs_m6) at m=300.
# ---------------------------------------------------------------------
emit_table <- function() {
  out <- c("\\begin{table}[H]\\centering\\small")
  out <- c(out, "\\caption{Realized product nuisance error $\\sqrt{m}\\,\\lVert\\widehat Q-Q_0\\rVert\\,\\lVert\\widehat g-g_0\\rVert$ across the PSU sweep $m_{\\mathrm{total}}=60\\to300$ (200 replicates, Design A), by learner rung. Columns give the mean realized product at each sampled-PSU count; the final column is the trend ratio relative to the baseline $m=60$ ($m_{\\mathrm{total}}=300$ value $\\div$ $m_{\\mathrm{total}}=60$ value, so $1$ means flat). The product is flat-to-declining at L2--L3 (the product rate is plausibly met, $o(m^{-1/2})$), but grows at L4 (rate fails) and at L1 (the misspecified parametric rung's approximation-error floor). This extends the single-$m$ diagnostic in Web Table~\\ref{tab:webD_rate}.}")
  out <- c(out, "\\label{tab:webD_ratesweep}")
  out <- c(out, "\\begin{tabular}{@{}l cccc c@{}}", "\\toprule")
  out <- c(out, hdr_row("Rung & $m{=}60$ & $m{=}120$ & $m{=}200$ & $m{=}300$ & Trend ($m{=}300/60$) \\\\"))
  out <- c(out, "\\midrule")
  out <- c(out, sprintf("\\multicolumn{6}{@{}l}{\\textit{%s}}\\\\",
                        "$\\sqrt{m}\\,\\lVert\\widehat Q-Q_0\\rVert\\,\\lVert\\widehat g-g_0\\rVert$ (mean over 200 reps)"))
  for (r in rung_ord) {
    vals <- sapply(m_ord, function(mm) d$mean_prod_sqrtm[d$rung == r & d$m_total == mm])
    trend <- d$trend_vs_m6[d$rung == r & d$m_total == 300]
    out <- c(out, sprintf("\\quad %s & %s & %s & %s & %s & %s \\\\",
                          rung_lab[[r]],
                          .f3(vals[1]), .f3(vals[2]), .f3(vals[3]), .f3(vals[4]),
                          .f2(trend)))
  }
  c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

tex <- paste(c(gen_banner("simulation/renderers/render_webD_ratesweep.R", src), "",
               paste(emit_table(), collapse = "\n")), collapse = "\n")
writeLines(tex, "outputs/tables/webD_ratesweep.tex")
cat(tex, "\n")
cat("\n--- wrote outputs/tables/webD_ratesweep.tex ---\n")

# ---------------------------------------------------------------------
# (2) Figure: one coloured line per rung; x = m_total, y = mean_prod_sqrtm.
#     Okabe-Ito palette; base graphics; appendix figure (colour is free).
# ---------------------------------------------------------------------
rung_col <- c(L1_param = "#E69F00", L2_smooth = "#009E73",
              L3_adaptive = "#56B4E9", L4_aggressive = "#D55E00")

draw_fig <- function() {
  yvals <- d$mean_prod_sqrtm
  ylim <- range(yvals)
  pad <- 0.04 * diff(ylim)
  ylim <- c(ylim[1] - pad, ylim[2] + pad)
  par(mar = c(4.4, 4.8, 1.0, 1.0), mgp = c(2.7, 0.8, 0), cex = 1.0)
  plot(NA, xlim = range(m_ord), ylim = ylim,
       xlab = "sampled PSUs (m_total)",
       ylab = expression(sqrt(m) * " " * "product nuisance error"),
       xaxt = "n")
  axis(1, at = m_ord, labels = m_ord)
  for (r in rung_ord) {
    sub <- d[d$rung == r, ]
    sub <- sub[order(sub$m_total), ]
    lines(sub$m_total, sub$mean_prod_sqrtm, col = rung_col[[r]], lwd = 2.4)
    points(sub$m_total, sub$mean_prod_sqrtm, col = rung_col[[r]],
           pch = 19, cex = 1.2)
  }
  legend("topleft", legend = rung_lab[rung_ord],
         col = rung_col[rung_ord], lwd = 2.4, pch = 19,
         bty = "n", cex = 0.9, seg.len = 1.6)
}

png_path <- "outputs/figures/fig_ratesweep.png"
pdf_path <- "outputs/figures/fig_ratesweep.pdf"
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

png(png_path, width = 6, height = 4, units = "in", res = 300)
draw_fig(); dev.off()
pdf(pdf_path, width = 6, height = 4)
draw_fig(); dev.off()
cat(sprintf("--- wrote %s and %s ---\n", png_path, pdf_path))
