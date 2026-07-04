# =====================================================================
# Nhanes/R/render_floorsens_tex.R
# R17 propensity-floor sensitivity for E4 (GDM history -> hypertension),
# primary Fully-Aware-CF arm -> outputs/tables/webF_floorsens.tex.
# Reads results/R17_floor_sensitivity.csv (one row per floor: 0.05, 0.025,
# 0.01). Shows that the primary point estimate is mildly floor-sensitive
# while the null conclusion is stable, and that the large exposed mass below
# the 0.05 floor motivates the overlap-restricted estimand reframe.
# =====================================================================
source("simulation/renderers/render_helpers_color.R")

src <- "results/R17_floor_sensitivity.csv"
d <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)

# order rows by descending floor (0.05 -> 0.025 -> 0.01)
d <- d[order(-d$floor), ]

# CI cell with proper LaTeX minus signs
ci <- function(l, u) sprintf("(%s, %s)",
                             gsub("-", "$-$", .f3(l)), gsub("-", "$-$", .f3(u)))

L <- c(
  gen_banner("Nhanes/renderers/render_floorsens_tex.R", src),
  "\\begin{table}[H]\\centering\\small",
  paste0("\\caption{Propensity-floor sensitivity for E4 (gestational-diabetes ",
         "history $\\to$ hypertension) under the primary Fully-Aware-CF estimator. ",
         "Columns: the propensity floor, the risk-difference point estimate, its ",
         "95\\% confidence interval, the weighted share of all units clipped at the ",
         "floor, and the weighted exposed mass below 0.05. The point estimate is ",
         "mildly floor-sensitive (0.016 at 0.05 to 0.022 at 0.01), but the null ",
         "conclusion (CI covers 0) is stable across all three floors. The large ",
         "weighted exposed mass below 0.05 ($\\approx 0.21$) motivates the ",
         "overlap-restricted (truncated) estimand reframe; see the by-example ",
         "overlap diagnostics in Web Table~\\ref{tab:webF_share}.}"),
  "\\label{tab:webF_floorsens}",
  "\\begin{tabular}{@{}l c c c c@{}}",
  "\\toprule",
  hdr_row("Floor & RD $\\widehat{b}$ & 95\\% CI & Share at floor & Exposed mass $<0.05$ \\\\"),
  "\\midrule")
for (i in seq_len(nrow(d))) {
  r <- d[i, ]
  L <- c(L, sprintf("%s & %s & %s & %s & %s \\\\",
                    .f3(r$floor),
                    gsub("-", "$-$", .bsign(r$b)),
                    ci(r$lcl, r$ucl),
                    .f3(r$share_clip_w),
                    .f3(r$expmass05_w)))
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")

writeLines(paste(L, collapse = "\n"), "outputs/tables/webF_floorsens.tex")
cat(paste(L, collapse = "\n"), "\n")
cat("\n--- wrote outputs/tables/webF_floorsens.tex ---\n")
