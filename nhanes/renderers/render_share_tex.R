# =====================================================================
# Nhanes/R/render_share_tex.R
# R17 by-example overlap shares at the 0.05 propensity floor ->
# outputs/tables/webF_share.tex. Reads results/R17_floor_share.csv
# (one row per example at floor 0.05). Shows the exposed-mass-below-floor
# overlap diagnostic that decides whether the overlap-restricted/truncated
# estimand caveat matters: it is material only for the rare exposure E4.
# =====================================================================
source("simulation/renderers/render_helpers_color.R")

src <- "results/R17_floor_share.csv"
d <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)

ord <- c("E2", "E3", "E4")
d <- d[match(ord, d$example), ]

# human-readable example labels (keyed by example)
exlab <- c(E2 = "E2: food insecurity $\\to$ depression",
           E3 = "E3: e-cigarette use $\\to$ hypertension",
           E4 = "E4: GDM history $\\to$ hypertension")

ci <- function(l, u) sprintf("(%s, %s)",
                             gsub("-", "$-$", .f3(l)), gsub("-", "$-$", .f3(u)))

L <- c(
  gen_banner("Nhanes/renderers/render_share_tex.R", src),
  "\\begin{table}[H]\\centering\\small",
  paste0("\\caption{Overlap diagnostics by example at the 0.05 propensity floor. ",
         "Columns: the example, the (weighted) exposure prevalence, the headline ",
         "overlap column---the weighted \\emph{exposed} mass with estimated ",
         "propensity below 0.05---the weighted share of \\emph{all} units clipped ",
         "at the floor, and the primary Fully-Aware-CF risk difference (95\\% CI). ",
         "The overlap-restricted (truncated) estimand caveat matters chiefly for the ",
         "rare exposure E4, where roughly 21\\% of the exposed weighted mass falls ",
         "below the floor; for E2 and E3 the exposed mass below 0.05 is minor ",
         "($\\approx 0.03$--$0.04$). The full E4 floor sweep is in Web ",
         "Table~\\ref{tab:webF_floorsens}.}"),
  "\\label{tab:webF_share}",
  "\\begin{tabular}{@{}l c c c c@{}}",
  "\\toprule",
  hdr_row("Example & Exposure prev. & Exposed mass $<0.05$ & Share $<0.05$ (all units) & RD (95\\% CI) \\\\"),
  "\\midrule")
for (i in seq_len(nrow(d))) {
  r <- d[i, ]
  L <- c(L, sprintf("%s & %s & %s & %s & %s %s \\\\",
                    exlab[[r$example]],
                    .f3(r$A_prev),
                    .f3(r$expmass05_w),
                    .f3(r$share_clip_w),
                    gsub("-", "$-$", .bsign(r$b)),
                    ci(r$lcl, r$ucl)))
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")

writeLines(paste(L, collapse = "\n"), "outputs/tables/webF_share.tex")
cat(paste(L, collapse = "\n"), "\n")
cat("\n--- wrote outputs/tables/webF_share.tex ---\n")
