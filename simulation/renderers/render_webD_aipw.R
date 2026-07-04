# =====================================================================
# render_webD_aipw.R  —  Web Appendix D: external survey-weighted AIPW
# benchmark (single-fit vs cross-fitted), generated from the locked ARC
# summary CSV. Mirrors simulation/renderers/render_webD_tex.R layout + the shared colour
# contract (simulation/renderers/render_helpers_color.R). Output: outputs/tables/webD_aipw.tex
# Uses ONLY se_type == "jkn" rows (the linearized SE is numerically
# identical to the jackknife here by construction).
# =====================================================================
source("simulation/renderers/render_helpers_color.R")

src <- "results/R15_aipw_benchmark_summary.csv"
d <- read.csv(src, check.names = FALSE)
d <- d[d$se_type == "jkn", ]

rung_lab <- c(L1_param = "L1 (GLM)", L2_smooth = "L2 (+GAM, MARS)",
              L3_adaptive = "L3 (+RF)", L4_aggressive = "L4 (deep RF)")
rung_ord <- names(rung_lab)
meth_ord <- c("AIPW-SF", "AIPW-CF")
meth_lab <- c("AIPW-SF" = "AIPW (single-fit)", "AIPW-CF" = "AIPW-CF (cross-fitted)")
scen_lab <- c(standard = "Design A ($H=10$ strata, $\\approx 6$ PSUs/stratum, $n=1{,}518$)",
              R1 = "Design B ($H=50$ strata, 2 PSUs/stratum, NHANES-like, $n=2{,}034$)")
scen_ord <- c("standard", "R1")

emit_design <- function(scn) {
  out <- sprintf("\\multicolumn{6}{@{}l}{\\textit{%s}}\\\\", scen_lab[[scn]])
  for (r in rung_ord) {
    out <- c(out, sprintf("\\multicolumn{6}{@{}l}{\\quad\\textit{%s}}\\\\", rung_lab[[r]]))
    sub <- d[d$scenario == scn & d$rung == r, ]
    for (m in meth_ord) {
      row <- sub[sub$method == m, ]
      if (nrow(row) == 0) next
      out <- c(out, sprintf("%s\\quad\\quad %s & %s & %s & %s & %s & %s \\\\",
                            cf_row(m, "AIPW-CF"), meth_lab[[m]],
                            .bsign(row$bias), .f3(row$emp_sd), .f3(row$mean_se),
                            .f2(row$se_ratio), cov_cell(row$coverage)))
    }
    if (r != rung_ord[length(rung_ord)]) out <- c(out, "\\addlinespace")
  }
  out
}

cap <- paste0(
  "External survey-weighted augmented inverse-probability-weighted (AIPW) ",
  "benchmark over 1{,}000 replicates: bias, empirical standard deviation (SD), ",
  "mean estimated standard error (SE), the SE-to-SD ratio, and Wald coverage of ",
  "nominal 0.95 intervals, for the single-fit AIPW estimator (AIPW-SF) and its ",
  "cross-fitted counterpart (AIPW-CF), under both survey designs. The linearized ",
  "standard error is numerically identical to the survey jackknife here by ",
  "construction, so only the jackknife rows are shown. Cross-fitting is necessary ",
  "for AIPW as well as for TMLE: AIPW-SF collapses at the aggressive rung L4 ",
  "(coverage 0.621 for Design A, 0.630 for Design B, with SE:SD $\\approx 0.43$--$0.47$), ",
  "while AIPW-CF holds nominal coverage at every rung ($\\approx 0.96$ at L4). Compare ",
  "with the corresponding TMLE arms in Web Table~\\ref{tab:webD_standard} (Design A) ",
  "and Web Table~\\ref{tab:webD_R1} (Design B). The primary cross-fitted arm is ",
  "shaded; coverage below 0.90 is highlighted. Monte Carlo standard error of ",
  "coverage is about 0.01.")

tab <- c(
  "\\begin{table}[H]\\centering\\small",
  sprintf("\\caption{%s}", cap),
  "\\label{tab:webD_aipw}",
  "\\begin{tabular}{@{}l rrrr r@{}}",
  "\\toprule",
  hdr_row("Estimator & Bias & SD & SE & SE:SD & Coverage \\\\"),
  "\\midrule")
for (i in seq_along(scen_ord)) {
  tab <- c(tab, emit_design(scen_ord[i]))
  if (i != length(scen_ord)) tab <- c(tab, "\\addlinespace")
}
tab <- c(tab, "\\bottomrule", "\\end{tabular}", "\\end{table}")

lines <- c(gen_banner("simulation/renderers/render_webD_aipw.R", src), "", tab)
writeLines(paste(lines, collapse = "\n"), "outputs/tables/webD_aipw.tex")
cat(paste(lines, collapse = "\n"))
cat("\n\n--- wrote outputs/tables/webD_aipw.tex ---\n")
