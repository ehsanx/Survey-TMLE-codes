source("simulation/renderers/render_helpers_color.R")

# =====================================================================
# render_webD_deployable.R  -- Web Appendix D table for the deployable
# non-Donsker libraries L6 / L7, generated from R21:
#   results/R21_deployable_cvcf_summary.csv  (no hand transcription).
# Output: outputs/tables/webD_deployable.tex  (label tab:webD_deploy)
#   The 3x2 decision view {single-fit, cluster-CV, cross-fit} x {TMLE, AIPW},
#   coverage AND se/sd, at the deployable libraries L6 (5-learner, tuned
#   XGBoost) and L7 (L6 + HAL), both designs. The se/sd column is mandatory:
#   the honest claim is ONE-SIDED (single-fit/cluster-CV under-cover via
#   se/sd<1; only cross-fitting reaches se/sd~1).
# =====================================================================

r21 <- read.csv("results/R21_deployable_cvcf_summary.csv", stringsAsFactors = FALSE)
g <- function(rg, sc, m, col) {
  r <- r21[[col]][r21$rung == rg & r21$scenario == sc & r21$method == m]
  if (length(r) == 1) r else NA_real_
}

blocks <- list(
  list(rg = "L6_deployable", sc = "standard", lab = "L6 (GLM, MARS, glmnet, tuned XGBoost, deep RF) --- Design A"),
  list(rg = "L6_deployable", sc = "R1",       lab = "L6 (as above) --- Design B"),
  list(rg = "L7_hal",        sc = "standard", lab = "L7 ($=$ L6 $+$ HAL) --- Design A"),
  list(rg = "L7_hal",        sc = "R1",       lab = "L7 ($=$ L6 $+$ HAL) --- Design B"))

# protocol -> (TMLE method, AIPW method)
prot <- list(
  list(lab = "single-fit",  t = "Fully-Aware",    a = "AIPW-SF", cf = FALSE),
  list(lab = "cluster-CV",  t = "Fully-Aware-CV", a = "AIPW-CV", cf = FALSE),
  list(lab = "cross-fit (CF)", t = "Fully-Aware-CF", a = "AIPW-CF", cf = TRUE))

out <- c("\\begin{table}[H]\\centering\\small")
out <- c(out, paste0(
  "\\caption{Deployable non-Donsker libraries: the 3$\\times$2 decision view ",
  "\\{single-fit, cluster-CV, cross-fit\\} $\\times$ \\{TMLE, AIPW\\} at L6 (a five-learner ",
  "Super Learner with tuned gradient boosting and a deep random forest --- the library one would ",
  "actually deploy) and L7 ($=$ L6 $+$ the highly adaptive lasso), over 1{,}000 replicates in both designs. ",
  "Each cell gives Wald coverage and the mean-SE-to-empirical-SD ratio (se/sd). Only cross-fitting reaches ",
  "se/sd $\\approx 1$ and nominal coverage; single-fit and cluster-CV stay anti-conservative (se/sd $<1$) for ",
  "\\emph{both} estimators. The L4 over-coverage optic (cross-fit se/sd $\\approx 1.44$) is gone on a deployable ",
  "library. Coverage $< 0.90$ is shaded; the primary cross-fit row is tinted.}"))
out <- c(out, "\\label{tab:webD_deploy}")
out <- c(out, "\\begin{tabular}{@{}l cc cc@{}}", "\\toprule")
out <- c(out, hdr_row("Protocol & \\multicolumn{2}{c}{TMLE} & \\multicolumn{2}{c}{AIPW} \\\\"))
out <- c(out, hdr_row(" & coverage & se/sd & coverage & se/sd \\\\"))
out <- c(out, "\\midrule")

for (b in blocks) {
  out <- c(out, sprintf("\\multicolumn{5}{@{}l}{\\textit{%s}}\\\\", b$lab))
  for (p in prot) {
    tint <- if (p$cf) "\\rowcolor{cfTint} " else ""
    out <- c(out, sprintf("%s%s & %s & %s & %s & %s \\\\",
                          tint, p$lab,
                          cov_cell(g(b$rg, b$sc, p$t, "coverage")), .f2(g(b$rg, b$sc, p$t, "se_ratio")),
                          cov_cell(g(b$rg, b$sc, p$a, "coverage")), .f2(g(b$rg, b$sc, p$a, "se_ratio"))))
  }
  if (!identical(b, blocks[[length(blocks)]])) out <- c(out, "\\addlinespace")
}
out <- c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")

tex <- paste(c(gen_banner("simulation/renderers/render_webD_deployable.R",
               "results/R21_deployable_cvcf_summary.csv"), out), collapse = "\n")
writeLines(tex, "outputs/tables/webD_deployable.tex")
cat(tex)
cat("\n\n--- wrote outputs/tables/webD_deployable.tex ---\n")
