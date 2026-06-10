# =====================================================================
# Nhanes/R/04_descriptives.R  —  survey-weighted Table 1 by exposure, per example
#
# Builds the design on the FULL MEC frame and subset()s to the domain (proper
# sub-population), then a survey-weighted descriptive table of covariates + outcome
# by exposure arm, with standardized mean differences (SMD). Saves a tidy
# data.frame per example + a combined markdown-ready table.
# =====================================================================

suppressMessages({ library(survey); library(tableone) })
REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
ana_dir <- file.path(REPO_ROOT, "Nhanes", "analytic")
res_dir <- file.path(REPO_ROOT, "Nhanes", "results"); dir.create(res_dir, showWarnings = FALSE)
options(survey.lonely.psu = "adjust")

ids <- c("E1", "E2", "E3", "E4")
table1 <- list()
for (id in ids) {
  obs  <- readRDS(file.path(ana_dir, paste0(id, "_imputed.rds")))
  covs <- attr(obs, "covs")
  obs$Af <- factor(obs$A, c(0, 1), c("Unexposed", "Exposed"))
  obs$Yf <- factor(obs$Y, c(0, 1), c("No", "Yes"))
  des  <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_POOLED,
                    nest = TRUE, data = obs)
  des_sub <- subset(des, inpop)
  vars <- c(covs, "Yf")
  t1 <- tryCatch(
    svyCreateTableOne(vars = vars, strata = "Af", data = des_sub, test = FALSE, addOverall = TRUE),
    error = function(e) { cat(id, "tableone error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(t1)) {
    m <- print(t1, smd = TRUE, printToggle = FALSE, varLabels = TRUE, noSpaces = TRUE)
    df <- data.frame(variable = rownames(m), m, row.names = NULL, check.names = FALSE)
    table1[[id]] <- df
    dir.create(file.path(res_dir, id), showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(tableone = t1, df = df, example = attr(obs, "example")),
            file.path(res_dir, id, "table1_descriptives.rds"))
    cat(sprintf("\n===== Table 1: %s (%s) =====\n", id, attr(obs, "example")))
    print(df, row.names = FALSE)
  }
}
cat("\nDONE. Per-example Table 1 saved under", res_dir, "\n")
