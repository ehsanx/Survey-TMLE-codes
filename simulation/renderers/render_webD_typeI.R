source("simulation/renderers/render_helpers_color.R")

# =====================================================================
# render_webD_typeI.R  -- Web Appendix D type-I / power table (R13),
# generated from results/R13_null_typeI_summary.csv (no hand transcription).
# Output: outputs/tables/webD_typeI.tex  (label tab:webD_typeI)
#   - one rung-blocked table; for each estimator, the rejection rate of
#     H0: psi = 0 at the true null (te = 0, target 0.05) and the power
#     (te = 0.3, true psi ~ 0.036), for Design A then Design B (NHANES-like).
# =====================================================================

src <- "results/R13_null_typeI_summary.csv"
d <- read.csv(src, check.names = FALSE)

rung_lab <- c(L1_param = "L1 (GLM)", L2_smooth = "L2 (+GAM, MARS)",
              L3_adaptive = "L3 (+RF)", L4_aggressive = "L4 (deep RF)")
rung_ord <- names(rung_lab)
meth_ord <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV",
              "Partially-Aware", "Non-Aware")
meth_lab <- c("Fully-Aware (single-fit)", "Fully-Aware-CF", "Fully-Aware-CV",
              "Partially-Aware", "Non-Aware")
names(meth_lab) <- meth_ord
primary <- "Fully-Aware-CF"

# helper: fetch a single reject_rate for (te, scenario, rung, method)
get_rr <- function(te, scn, r, m) {
  row <- d[d$te == te & d$scenario == scn & d$rung == r & d$method == m, ]
  if (nrow(row) == 0) return(NA_real_)
  row$reject_rate[1]
}

out <- c("\\begin{table}[H]\\centering\\small")
out <- c(out, paste0(
  "\\caption{Type-I error (rejection of $H_0\\!: \\psi = 0$ at the true null, ",
  "target 0.05) and power (true $\\psi \\approx 0.036$) over 1{,}000 replicates, ",
  "by rung and design. The primary Fully-Aware-CF controls type-I where the ",
  "design-naive arms over-reject and the single-fit estimator collapses at L4; ",
  "cells with type-I $> 0.075$ are shaded. The small excess for the cross-fitted ",
  "arm at the Design B parametric rung reflects the known finite-sample bias under the misspecified DGP ",
  "(Web Table~\\ref{tab:webD_standard}), not a variance defect. Fully-Aware-CV is ",
  "defined only at the multi-learner rungs L2--L3.}"))
out <- c(out, "\\label{tab:webD_typeI}")
out <- c(out, "\\begin{tabular}{@{}l cc cc@{}}", "\\toprule")
out <- c(out, "\\multicolumn{1}{@{}l}{} & \\multicolumn{2}{c}{Design A} & \\multicolumn{2}{c}{Design B (NHANES-like)} \\\\")
out <- c(out, "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}")
out <- c(out, hdr_row("Estimator & type-I & power & type-I & power \\\\"))
out <- c(out, "\\midrule")

for (r in rung_ord) {
  out <- c(out, sprintf("\\multicolumn{5}{@{}l}{\\textit{%s}}\\\\", rung_lab[[r]]))
  for (m in meth_ord) {
    # skip absent rows (e.g. Fully-Aware-CV at L1/L4)
    if (is.na(get_rr(0, "standard", r, m)) && is.na(get_rr(0, "R1", r, m))) next
    a_t1 <- typeI_cell(get_rr(0,   "standard", r, m))
    a_pw <- .f3(      get_rr(0.3, "standard", r, m))
    b_t1 <- typeI_cell(get_rr(0,   "R1",       r, m))
    b_pw <- .f3(      get_rr(0.3, "R1",       r, m))
    out <- c(out, sprintf("%s\\quad %s & %s & %s & %s & %s \\\\",
                          cf_row(m, primary), meth_lab[[m]],
                          a_t1, a_pw, b_t1, b_pw))
  }
  if (r != rung_ord[length(rung_ord)]) out <- c(out, "\\addlinespace")
}
out <- c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")

tex <- paste(c(gen_banner("simulation/renderers/render_webD_typeI.R", src), out), collapse = "\n")
writeLines(tex, "outputs/tables/webD_typeI.tex")
cat(tex)
cat("\n\n--- wrote outputs/tables/webD_typeI.tex ---\n")
