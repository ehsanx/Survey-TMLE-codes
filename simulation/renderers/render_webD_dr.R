source("simulation/renderers/render_helpers_color.R")
# =====================================================================
# render_webD_dr.R  --  Web Appendix D: R14 double-robustness factorial
# (Q and g each correctly/incorrectly specified) + the informative-sampling
# unweighted-out-of-fold threat. Output: outputs/tables/webD_dr.tex.
# Source: results/R14_dr_factorial_summary.csv (locked ARC run).
# Two stacked panels in ONE table:
#   Panel A: CF-u (paper default), L1, non-informative -- 2x2 (Q,g) DR cells,
#            Bias | Coverage. Bias ~ 0 unless BOTH wrong (WW).
#   Panel B: informative sampling, L1 -- CF-u Bias | CF-u Coverage | CF-w Bias
#            (n_diverged). CF-u (unweighted, default) stays unbiased & well-
#            covered and converges 1000/1000; weighted CF-w diverges -> the
#            unweighted-OOF threat does NOT materialise; weighted fits unstable.
# =====================================================================
src <- "results/R14_dr_factorial_summary.csv"
d <- read.csv(src, check.names = FALSE)

# (Q,g) cells in display order, with their plain-English labels
cells <- list(
  c(q = "C", g = "C", lab = "Q correct, g correct"),
  c(q = "C", g = "W", lab = "Q correct, g wrong"),
  c(q = "W", g = "C", lab = "Q wrong, g correct"),
  c(q = "W", g = "W", lab = "Q wrong, g wrong")
)

pick <- function(qv, gv, samp, meth) {
  r <- d[d$q_spec == qv & d$g_spec == gv & d$sampling == samp &
         d$rung == "L1_param" & d$method == meth, ]
  stopifnot(nrow(r) == 1)
  r
}

out <- c(gen_banner("simulation/renderers/render_webD_dr.R", src))
out <- c(out, "\\begin{table}[H]\\centering\\small")
out <- c(out, paste0(
  "\\caption{Double robustness and the informative-sampling unweighted-out-of-fold ",
  "threat (R14, parametric rung L1, 1{,}000 target replicates per cell). ",
  "\\textbf{Panel A} is the $2\\times2$ double-robustness factorial under ",
  "non-informative sampling for the unweighted cross-fitted estimator CF-u (the ",
  "paper default): the outcome model $Q$ and the propensity $g$ are each correctly ",
  "or incorrectly specified, and bias stays near zero unless \\emph{both} are wrong. ",
  "\\textbf{Panel B} probes the threat that unweighted out-of-fold nuisance fits ",
  "could bias the estimator under informative sampling; it reports CF-u bias and ",
  "Wald coverage of nominal $0.95$ intervals alongside the weighted variant CF-w ",
  "(with its number of diverged replicates). Under the paper's unweighted-OOF ",
  "cross-fitting the estimator remains doubly robust and the informative-sampling ",
  "threat does not arise (CF-u converges $1{,}000/1{,}000$), whereas the weighted ",
  "nuisance fits are numerically unstable and diverge in $100$--$200$ replicates. ",
  "Companion full-factorial results across designs and rungs are in Web ",
  "Table~\\ref{tab:webD_standard}.}"))
out <- c(out, "\\label{tab:webD_dr}")
out <- c(out, "\\begin{tabular}{@{}l rr r@{}}", "\\toprule")
out <- c(out, hdr_row("Specification & Bias & Coverage & CF-w bias ($n_{\\mathrm{div}}$) \\\\"))
out <- c(out, "\\midrule")

# ---- Panel A: CF-u, L1, non-informative -- Bias | Coverage (cols 2-3) ----
out <- c(out, "\\multicolumn{4}{@{}l}{\\textit{Panel A: double robustness (CF-u, L1, non-informative sampling)}}\\\\")
for (cc in cells) {
  r <- pick(cc[["q"]], cc[["g"]], "noninfo", "CF-u")
  out <- c(out, sprintf("\\quad %s & %s & %s & \\\\",
                        cc[["lab"]], .bsign(r$bias), cov_cell(r$coverage)))
}
out <- c(out, "\\addlinespace")

# ---- Panel B: informative -- CF-u Bias | CF-u Coverage | CF-w Bias (n_div) ----
out <- c(out, "\\multicolumn{4}{@{}l}{\\textit{Panel B: informative sampling, L1 (does the unweighted-OOF threat arise?)}}\\\\")
for (cc in cells) {
  ru <- pick(cc[["q"]], cc[["g"]], "info", "CF-u")
  rw <- pick(cc[["q"]], cc[["g"]], "info", "CF-w")
  out <- c(out, sprintf("\\quad %s & %s & %s & %s (%d) \\\\",
                        cc[["lab"]], .bsign(ru$bias), cov_cell(ru$coverage),
                        .bsign(rw$bias), as.integer(rw$n_diverged)))
}

out <- c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")

tex <- paste(out, collapse = "\n")
writeLines(tex, "outputs/tables/webD_dr.tex")
cat(tex)
cat("\n\n--- wrote outputs/tables/webD_dr.tex ---\n")

# ---- echo key numbers for cross-check ----
cat("\nKEY NUMBERS:\n")
cat("Panel A (CF-u, noninfo, L1) bias  CC/CW/WC/WW:",
    .bsign(pick("C","C","noninfo","CF-u")$bias),
    .bsign(pick("C","W","noninfo","CF-u")$bias),
    .bsign(pick("W","C","noninfo","CF-u")$bias),
    .bsign(pick("W","W","noninfo","CF-u")$bias), "\n")
cat("Panel A (CF-u, noninfo, L1) cov   CC/CW/WC/WW:",
    .f3(pick("C","C","noninfo","CF-u")$coverage),
    .f3(pick("C","W","noninfo","CF-u")$coverage),
    .f3(pick("W","C","noninfo","CF-u")$coverage),
    .f3(pick("W","W","noninfo","CF-u")$coverage), "\n")
cat("Panel B (info, L1) CF-u bias      CC/CW/WC/WW:",
    .bsign(pick("C","C","info","CF-u")$bias),
    .bsign(pick("C","W","info","CF-u")$bias),
    .bsign(pick("W","C","info","CF-u")$bias),
    .bsign(pick("W","W","info","CF-u")$bias), "\n")
cat("Panel B (info, L1) CF-w bias      CC/CW/WC/WW:",
    .bsign(pick("C","C","info","CF-w")$bias),
    .bsign(pick("C","W","info","CF-w")$bias),
    .bsign(pick("W","C","info","CF-w")$bias),
    .bsign(pick("W","W","info","CF-w")$bias), "\n")
cat("Panel B (info, L1) CF-w n_div     CC/CW/WC/WW:",
    pick("C","C","info","CF-w")$n_diverged,
    pick("C","W","info","CF-w")$n_diverged,
    pick("W","C","info","CF-w")$n_diverged,
    pick("W","W","info","CF-w")$n_diverged, "\n")
cat("Panel B (info, L1) CF-u n_div     CC/CW/WC/WW:",
    pick("C","C","info","CF-u")$n_diverged,
    pick("C","W","info","CF-u")$n_diverged,
    pick("W","C","info","CF-u")$n_diverged,
    pick("W","W","info","CF-u")$n_diverged, "\n")
