# =====================================================================
# install_packages.R  —  verify (and only if needed, install) the R packages
# this repository uses. Run once from the repository root:
#   Rscript install_packages.R
# `statmod` is intentionally NOT required (Gauss-Hermite is computed in base R
# in codes/dgp.R). `nhanesA` is needed only for the NHANES download step.
# =====================================================================
pkgs <- c("SuperLearner", "tmle", "survey", "surveyCV",
          "earth", "gam", "glmnet", "ranger",          # estimator engine + simulation ladder
          "mice", "dplyr", "tidyr", "purrr", "nhanesA", # NHANES build / impute / download
          "rsimsum")                                    # optional: simulation summaries
ip   <- rownames(installed.packages())
need <- setdiff(pkgs, ip)
if (length(need)) {
  cat("installing missing:", paste(need, collapse = ", "), "\n")
  install.packages(need, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages already present (no install needed).\n")
}
for (p in pkgs) cat(sprintf("%-14s %s\n", p,
  if (p %in% rownames(installed.packages())) as.character(packageVersion(p)) else "MISSING"))
