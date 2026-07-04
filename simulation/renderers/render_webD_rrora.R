source("simulation/renderers/render_helpers_color.R")

# =====================================================================
# render_webD_rrora.R  -- Web Appendix D table for the RR/OR + continuous
# corollary (Corollary cor:rrora), generated from the A21 demonstration:
#   results/A21_rr_or_continuous_{rr,rd}_summary.csv  (no hand transcription).
# Output: outputs/tables/webD_rrora.tex  (label tab:webD_rrora)
#   Coverage of the delta-method 95% Wald intervals for RR, OR, and RD, by
#   outcome (binary / continuous) x rung (L1, L2) x arm (single-fit / cross-fit).
# =====================================================================

rr <- read.csv("results/A21_rr_or_continuous_rr_summary.csv", check.names = FALSE)
rd <- read.csv("results/A21_rr_or_continuous_rd_summary.csv", check.names = FALSE)

cov_of <- function(df, oc, rg, m, es) {
  r <- df[df$outcome == oc & df$rung == rg & df$method == m & df$estimand == es, ]
  if (!nrow(r)) return(NA_real_); r$coverage[1]
}
ser_of <- function(df, oc, rg, m, es) {
  r <- df[df$outcome == oc & df$rung == rg & df$method == m & df$estimand == es, ]
  if (!nrow(r)) return(NA_real_); r$se_ratio[1]
}

blocks <- list(
  list(oc = "binary",     rg = "L1_param",  lab = "Binary outcome, L1 (GLM)"),
  list(oc = "binary",     rg = "L2_smooth", lab = "Binary outcome, L2 (+GAM, MARS)"),
  list(oc = "continuous", rg = "L1_param",  lab = "Continuous outcome, L1 (GLM)"),
  list(oc = "continuous", rg = "L2_smooth", lab = "Continuous outcome, L2 (+GAM, MARS)"))
meth_ord <- c("Fully-Aware", "Fully-Aware-CF")
meth_lab <- c("Fully-Aware (single-fit)", "Fully-Aware-CF"); names(meth_lab) <- meth_ord
primary <- "Fully-Aware-CF"

out <- c("\\begin{table}[H]\\centering\\small")
out <- c(out, paste0(
  "\\caption{Coverage of the delta-method 95\\% Wald intervals for the risk ratio (RR), odds ratio (OR), ",
  "and risk difference (RD) over 800 replicates on the standard design, for the ",
  "single-fit (Fully-Aware) and cross-fitted (Fully-Aware-CF) arms (Corollary~\\ref{cor:rrora}). With a flexible ",
  "ensemble (L2) the ratio estimands and the continuous-outcome RD attain nominal coverage, and ",
  "$\\mathrm{se}/\\mathrm{sd}\\approx1$ confirms the delta-method variance; the parametric rung (L1) inherits the ",
  "same Kang--Schafer misspecification bias as the risk difference. The odds ratio is specific to a binary ",
  "outcome (`---' in the continuous case, where RD reads as the mean difference and RR as the ratio of means). ",
  "Coverage $< 0.90$ is shaded.}"))
out <- c(out, "\\label{tab:webD_rrora}")
out <- c(out, "\\begin{tabular}{@{}l ccc c@{}}", "\\toprule")
out <- c(out, hdr_row("Estimator & RR & OR & RD & se/sd (RR) \\\\"))
out <- c(out, "\\midrule")

for (b in blocks) {
  out <- c(out, sprintf("\\multicolumn{5}{@{}l}{\\textit{%s}}\\\\", b$lab))
  for (m in meth_ord) {
    rrc <- cov_cell(cov_of(rr, b$oc, b$rg, m, "RR"))
    orc <- if (b$oc == "binary") cov_cell(cov_of(rr, b$oc, b$rg, m, "OR")) else "---"
    rdc <- cov_cell(cov_of(rd, b$oc, b$rg, m, "RD"))
    ser <- .f2(ser_of(rr, b$oc, b$rg, m, "RR"))
    out <- c(out, sprintf("%s%s & %s & %s & %s & %s \\\\",
                          cf_row(m, primary), meth_lab[[m]], rrc, orc, rdc, ser))
  }
  if (!identical(b, blocks[[length(blocks)]])) out <- c(out, "\\addlinespace")
}
out <- c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")

tex <- paste(c(gen_banner("simulation/renderers/render_webD_rrora.R",
               "results/A21_rr_or_continuous_{rr,rd}_summary.csv"), out), collapse = "\n")
writeLines(tex, "outputs/tables/webD_rrora.tex")
cat(tex)
cat("\n\n--- wrote outputs/tables/webD_rrora.tex ---\n")
