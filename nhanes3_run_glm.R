# =================================================================
# Script to Fit Survey-Weighted GLMs for OR, RR, and RD (All Awareness Levels)
# =================================================================

# --- 1. SETUP ---

# Load necessary libraries
# install.packages(c("dplyr", "survey"))
library(dplyr)
library(survey)

# Create a directory to store results if it doesn't exist
dir.create("base_results", showWarnings = FALSE)

# --- Load and Prepare Data ---
imputed_data_complete <- readRDS("E:/GitHub/survey-tmle/Real data analysis/server codes/nhanes_singly_imputed_data.rds")

# Define the outcome, exposure, and covariates
Y <- as.numeric(imputed_data_complete$HYPERTENSION) - 1
A <- as.numeric(imputed_data_complete$GDM_HISTORY) - 1
W_df <- imputed_data_complete %>%
  select(RIDAGEYR, RIDRETH3, DMDEDUC2, INDFMPIR, BMXBMI, RHQ172) %>%
  mutate(across(where(is.factor), as.numeric))

# Combine all necessary variables into a single data frame
analysis_data <- data.frame(
  Y, A, W_df,
  SDMVPSU = imputed_data_complete$SDMVPSU,
  SDMVSTRA = imputed_data_complete$SDMVSTRA,
  WTMEC2YR = imputed_data_complete$WTMEC2YR
)


# --- 2. MODEL FITTING AND EXTRACTION ---

# =========================================================
# a. Fully-Aware: Using weights, strata, and clusters
# =========================================================
cat("--- Running Fully-Aware Models ---\n")

# Create the full survey design object
svy_design_full <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  data = analysis_data,
  nest = TRUE
)

# Fit models using the full survey design
fit_or_full <- svyglm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, design = svy_design_full, family = quasibinomial(link = "logit"))
fit_rr_full <- svyglm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, design = svy_design_full, family = quasipoisson(link = "log"))
fit_rd_full <- svyglm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, design = svy_design_full, family = gaussian(link = "identity"))

# Extract and combine results
results_or_full <- tibble(Parameter = "OR", Estimate = exp(coef(fit_or_full)["A"]), `95% CI Lower` = exp(confint(fit_or_full)["A", 1]), `95% CI Upper` = exp(confint(fit_or_full)["A", 2]), `P-Value` = summary(fit_or_full)$coefficients["A", 4], `Std. Error` = summary(fit_or_full)$coefficients["A", 2])
results_rr_full <- tibble(Parameter = "RR", Estimate = exp(coef(fit_rr_full)["A"]), `95% CI Lower` = exp(confint(fit_rr_full)["A", 1]), `95% CI Upper` = exp(confint(fit_rr_full)["A", 2]), `P-Value` = summary(fit_rr_full)$coefficients["A", 4], `Std. Error` = summary(fit_rr_full)$coefficients["A", 2])
results_rd_full <- tibble(Parameter = "ATE", Estimate = coef(fit_rd_full)["A"], `95% CI Lower` = confint(fit_rd_full)["A", 1], `95% CI Upper` = confint(fit_rd_full)["A", 2], `P-Value` = summary(fit_rd_full)$coefficients["A", 4], `Std. Error` = summary(fit_rd_full)$coefficients["A", 2])

final_results_full <- bind_rows(results_rd_full, results_rr_full, results_or_full) %>% mutate(Analysis = "Survey GLM Fully-Aware", .before = 1)
print(final_results_full)
saveRDS(final_results_full, file = "base_results/survey_glm_fully_aware.rds")


# =========================================================
# b. Partially-Aware: Using weights only
# =========================================================
cat("\n--- Running Partially-Aware Models ---\n")

# Create a survey design object with weights only (no strata or clusters)
svy_design_partial <- svydesign(
  id = ~1,
  strata = NULL,
  weights = ~WTMEC2YR,
  data = analysis_data
)

# Fit models using the partial survey design
fit_or_partial <- svyglm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, design = svy_design_partial, family = quasibinomial(link = "logit"))
fit_rr_partial <- svyglm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, design = svy_design_partial, family = quasipoisson(link = "log"))
fit_rd_partial <- svyglm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, design = svy_design_partial, family = gaussian(link = "identity"))

# Extract and combine results
results_or_partial <- tibble(Parameter = "OR", Estimate = exp(coef(fit_or_partial)["A"]), `95% CI Lower` = exp(confint(fit_or_partial)["A", 1]), `95% CI Upper` = exp(confint(fit_or_partial)["A", 2]), `P-Value` = summary(fit_or_partial)$coefficients["A", 4], `Std. Error` = summary(fit_or_partial)$coefficients["A", 2])
results_rr_partial <- tibble(Parameter = "RR", Estimate = exp(coef(fit_rr_partial)["A"]), `95% CI Lower` = exp(confint(fit_rr_partial)["A", 1]), `95% CI Upper` = exp(confint(fit_rr_partial)["A", 2]), `P-Value` = summary(fit_rr_partial)$coefficients["A", 4], `Std. Error` = summary(fit_rr_partial)$coefficients["A", 2])
results_rd_partial <- tibble(Parameter = "ATE", Estimate = coef(fit_rd_partial)["A"], `95% CI Lower` = confint(fit_rd_partial)["A", 1], `95% CI Upper` = confint(fit_rd_partial)["A", 2], `P-Value` = summary(fit_rd_partial)$coefficients["A", 4], `Std. Error` = summary(fit_rd_partial)$coefficients["A", 2])

final_results_partial <- bind_rows(results_rd_partial, results_rr_partial, results_or_partial) %>% mutate(Analysis = "Survey GLM Partially-Aware", .before = 1)
print(final_results_partial)
saveRDS(final_results_partial, file = "base_results/survey_glm_partially_aware.rds")


# =========================================================
# c. Non-Aware: Using no survey features (standard GLM)
# =========================================================
cat("\n--- Running Non-Aware Models ---\n")

# Fit standard GLMs with no survey adjustments
fit_or_non <- glm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, data = analysis_data, family = binomial(link = "logit"))
fit_rr_non <- glm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, data = analysis_data, family = poisson(link = "log"))
fit_rd_non <- glm(Y ~ A + RIDAGEYR + RIDRETH3 + DMDEDUC2 + INDFMPIR + BMXBMI + RHQ172, data = analysis_data, family = gaussian(link = "identity"))

# Extract and combine results
results_or_non <- tibble(Parameter = "OR", Estimate = exp(coef(fit_or_non)["A"]), `95% CI Lower` = exp(confint.default(fit_or_non)["A", 1]), `95% CI Upper` = exp(confint.default(fit_or_non)["A", 2]), `P-Value` = summary(fit_or_non)$coefficients["A", 4], `Std. Error` = summary(fit_or_non)$coefficients["A", 2])
results_rr_non <- tibble(Parameter = "RR", Estimate = exp(coef(fit_rr_non)["A"]), `95% CI Lower` = exp(confint.default(fit_rr_non)["A", 1]), `95% CI Upper` = exp(confint.default(fit_rr_non)["A", 2]), `P-Value` = summary(fit_rr_non)$coefficients["A", 4], `Std. Error` = summary(fit_rr_non)$coefficients["A", 2])
results_rd_non <- tibble(Parameter = "ATE", Estimate = coef(fit_rd_non)["A"], `95% CI Lower` = confint.default(fit_rd_non)["A", 1], `95% CI Upper` = confint.default(fit_rd_non)["A", 2], `P-Value` = summary(fit_rd_non)$coefficients["A", 4], `Std. Error` = summary(fit_rd_non)$coefficients["A", 2])

final_results_non <- bind_rows(results_rd_non, results_rr_non, results_or_non) %>% mutate(Analysis = "Survey GLM Non-Aware", .before = 1)
print(final_results_non)
saveRDS(final_results_non, file = "base_results/survey_glm_non_aware.rds")

cat("\nAll GLM analyses complete.\n")
