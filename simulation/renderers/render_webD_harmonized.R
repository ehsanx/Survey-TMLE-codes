# =====================================================================
# render_webD_harmonized.R  —  Web Appendix D: harmonized-floor robustness
# check. With a SINGLE harmonized propensity floor (0.05) + Q-truncation
# applied uniformly across ALL arms, does the single-fit collapse at L4
# survive? Generated from the locked ARC summary CSV; mirrors the shared
# colour contract (simulation/renderers/render_helpers_color.R).
# Output: outputs/tables/webD_harmonized.tex
# Focused table: rungs L3 and L4 only (where the gap matters), per design.
# =====================================================================
source("simulation/renderers/render_helpers_color.R")

src <- "results/R16_harmonized_sim_summary.csv"
d <- read.csv(src, check.names = FALSE)

rung_lab <- c(L3_adaptive = "L3 (+RF)", L4_aggressive = "L4 (deep RF)")
rung_ord <- names(rung_lab)
# Per-rung row sets: CV-h exists only at the multi-learner rung L3.
rung_meth <- list(
  L3_adaptive   = c("Fully-Aware-h", "Fully-Aware-CF-h", "Fully-Aware-CV-h"),
  L4_aggressive = c("Fully-Aware-h", "Fully-Aware-CF-h"))
meth_lab <- c("Fully-Aware-h"    = "Fully-Aware-h (single-fit)",
              "Fully-Aware-CF-h" = "Fully-Aware-CF-h",
              "Fully-Aware-CV-h" = "Fully-Aware-CV-h")
scen_lab <- c(standard = "Design A ($H=10$ strata, $\\approx 6$ PSUs/stratum, $n=1{,}518$)",
              R1 = "Design B ($H=50$ strata, 2 PSUs/stratum, NHANES-like, $n=2{,}034$)")
scen_ord <- c("standard", "R1")

emit_design <- function(scn) {
  out <- sprintf("\\multicolumn{5}{@{}l}{\\textit{%s}}\\\\", scen_lab[[scn]])
  for (r in rung_ord) {
    out <- c(out, sprintf("\\multicolumn{5}{@{}l}{\\quad\\textit{%s}}\\\\", rung_lab[[r]]))
    sub <- d[d$scenario == scn & d$rung == r, ]
    for (m in rung_meth[[r]]) {
      row <- sub[sub$method == m, ]
      if (nrow(row) == 0) next
      out <- c(out, sprintf("%s\\quad\\quad %s & %s & %s & %s & %s \\\\",
                            cf_row(m, "Fully-Aware-CF-h"), meth_lab[[m]],
                            .bsign(row$bias), .f3(row$emp_sd),
                            .f2(row$se_ratio), cov_cell(row$coverage)))
    }
    if (r != rung_ord[length(rung_ord)]) out <- c(out, "\\addlinespace")
  }
  out
}

cap <- paste0(
  "Harmonized-floor robustness check (1{,}000 replicates). With a \\emph{single} ",
  "harmonized propensity floor of 0.05 and identical $Q$-truncation applied ",
  "uniformly across \\emph{all} arms, the single-fit collapse at the aggressive ",
  "rung L4 \\emph{survives}: the single-fit Fully-Aware-h attains coverage only ",
  "0.388 (Design A) and 0.328 (Design B), whereas its cross-fitted counterpart ",
  "Fully-Aware-CF-h holds nominal coverage (0.986 Design A, 0.991 Design B). The ",
  "headline gap of Figure~1 is therefore a cross-fitting effect, not an artifact ",
  "of differing propensity floors. Only the rungs where the gap matters (L3, L4) ",
  "are shown; Fully-Aware-CV-h is defined only at the multi-learner rung L3. ",
  "Columns are bias, empirical SD, the SE-to-SD ratio, and Wald coverage of ",
  "nominal 0.95 intervals. The primary cross-fitted arm is shaded; coverage ",
  "below 0.90 is highlighted. At L1 the single-fit Fully-Aware-h has substantial ",
  "non-convergence ($\\approx 17\\%$ of replicates under Design A; those rungs are ",
  "omitted here). Monte Carlo standard error of coverage is about 0.01. Compare ",
  "with the unharmonized arms in Web Tables~\\ref{tab:webD_standard} (Design A) ",
  "and~\\ref{tab:webD_R1} (Design B).")

tab <- c(
  "\\begin{table}[H]\\centering\\small",
  sprintf("\\caption{%s}", cap),
  "\\label{tab:webD_harmonized}",
  "\\begin{tabular}{@{}l rr r r@{}}",
  "\\toprule",
  hdr_row("Estimator & Bias & SD & SE:SD & Coverage \\\\"),
  "\\midrule")
for (i in seq_along(scen_ord)) {
  tab <- c(tab, emit_design(scen_ord[i]))
  if (i != length(scen_ord)) tab <- c(tab, "\\addlinespace")
}
tab <- c(tab, "\\bottomrule", "\\end{tabular}", "\\end{table}")

lines <- c(gen_banner("simulation/renderers/render_webD_harmonized.R", src), "", tab)
writeLines(paste(lines, collapse = "\n"), "outputs/tables/webD_harmonized.tex")
cat(paste(lines, collapse = "\n"))
cat("\n\n--- wrote outputs/tables/webD_harmonized.tex ---\n")
