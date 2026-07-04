# =====================================================================
# Nhanes/R/render_aipwe1_tex.R
# R15 external survey-weighted AIPW on E1 (short sleep -> obesity) ->
# outputs/tables/webF_aipwe1.tex. Reads results/R15_aipw_nhanes_E1.csv
# (AIPW-SF single-fit + AIPW-CF cross-fitted). The cross-fitted AIPW agrees
# with the primary TMLE-CF, and the jackknife and linearization SEs
# coincide on this well-overlapped example. The primary cross-fitted arm
# (AIPW-CF) is tinted.
# =====================================================================
source("simulation/renderers/render_helpers_color.R")

src <- "results/R15_aipw_nhanes_E1.csv"
d <- read.csv(src, stringsAsFactors = FALSE, check.names = FALSE)

ord <- c("AIPW-SF", "AIPW-CF")
d <- d[match(ord, d$method), ]

ci <- function(l, u) sprintf("(%s, %s)",
                             gsub("-", "$-$", .f3(l)), gsub("-", "$-$", .f3(u)))

L <- c(
  gen_banner("Nhanes/renderers/render_aipwe1_tex.R", src),
  "\\begin{table}[H]\\centering\\small",
  paste0("\\caption{External survey-weighted augmented inverse-probability-weighted ",
         "(AIPW) estimator on E1 (short sleep $\\to$ obesity), as a methods cross-check ",
         "against the primary TMLE. Columns: the AIPW arm, the risk-difference point ",
         "estimate, the jackknife standard error, the linearization standard error, and ",
         "the 95\\% confidence interval. The cross-fitted AIPW-CF estimate (0.036) agrees ",
         "with the primary Fully-Aware-CF TMLE (0.036; Web Table~\\ref{tab:webF_arms}), ",
         "and the jackknife and linearization SEs coincide to three decimals on this ",
         "well-overlapped example. AIPW-CF (cross-fitted) is the primary arm.}"),
  "\\label{tab:webF_aipwe1}",
  "\\begin{tabular}{@{}l c c c c@{}}",
  "\\toprule",
  hdr_row("Estimator & RD $\\widehat{b}$ & Jackknife SE & Linearization SE & 95\\% CI \\\\"),
  "\\midrule")
for (i in seq_len(nrow(d))) {
  r <- d[i, ]
  L <- c(L, sprintf("%s%s & %s & %s & %s & %s \\\\",
                    cf_row(r$method, "AIPW-CF"),
                    r$method,
                    gsub("-", "$-$", .bsign(r$b)),
                    .f3(r$se_jkn),
                    .f3(r$se_lin),
                    ci(r$lcl, r$ucl)))
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")

writeLines(paste(L, collapse = "\n"), "outputs/tables/webF_aipwe1.tex")
cat(paste(L, collapse = "\n"), "\n")
cat("\n--- wrote outputs/tables/webF_aipwe1.tex ---\n")
