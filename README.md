# NHANES Analysis: GDM and Hypertension Risk

This repository contains the R scripts for an analysis pipeline investigating the association between a history of gestational diabetes mellitus (GDM) and the risk of hypertension later in life, using data from the National Health and Nutrition Examination Survey (NHANES).

The analysis is structured into three main scripts that should be run sequentially.

## Project Overview

The goal of this project is to estimate the causal effect of GDM on hypertension while accounting for the complex survey design of NHANES. The pipeline performs the following major steps:

1.  **Data Acquisition and Preparation**: Downloads, merges, and cleans NHANES data from 2007-2018.
2.  **Imputation**: Handles missing data in the analytic sample using single imputation.
3.  **Modeling**:
    * Fits various Targeted Maximum Likelihood Estimation (TMLE) models to estimate causal effects.
    * Fits several Generalized Linear Models (GLMs) as a comparison.
4.  **Results**: Saves all model outputs as `.rds` files for further analysis or reporting.

---

## Scripts

### 1. `1_save_analytic_data.R`

This is the foundational script that prepares the data for all subsequent analyses.

**Key Functions:**

* **Download Data**: Downloads six cycles of NHANES data (2007-2018) using the `nhanesA` package. It fetches data from the Demographics (DEMO), Diabetes (DIQ), Reproductive Health (RHQ), Body Measures (BMX), and Blood Pressure (BPX) questionnaires.
* **Merge and Clean**: Merges the datasets and recodes variables to create a consistent, analysis-ready dataset. This includes:
    * Creating the binary outcome variable `HYPERTENSION`.
    * Creating the binary exposure variable `GDM_HISTORY`.
    * Harmonizing and simplifying categorical variables like race/ethnicity (`RIDRETH3`) and education (`DMDEDUC2`).
* **Define Eligibility**: Creates an eligibility flag (`is_eligible`) to include only non-pregnant females aged 20 and over without pre-existing diabetes.
* **Impute Missing Data**:
    * Filters the data to an `analytic_sample` of eligible participants with non-missing exposure and outcome.
    * Uses the `mice` package to perform a **single imputation** (`m=1`) for missing values in covariates (e.g., `BMXBMI`, `INDFMPIR`). Predictive mean matching (`pmm`) and polytomous regression (`polyreg`) are used as imputation methods.
* **Save Output**: Saves the final, clean, and imputed dataset as `nhanes_singly_imputed_data.rds`. This file is the input for the two modeling scripts.

### 2. `2_run_tmle.R`

This script uses the imputed data to perform a series of Targeted Maximum Likelihood Estimation (TMLE) analyses to estimate the causal effect of GDM on hypertension. TMLE is a doubly-robust method that uses machine learning to reduce bias.

**Key Functions:**

* **Load Data**: Loads the `nhanes_singly_imputed_data.rds` file.
* **Setup SuperLearner**: Defines several libraries of machine learning algorithms (e.g., `SL.glm`, `SL.gam`, `SL.earth`, `SL.glmnet`) to be used by TMLE for modeling the outcome and exposure mechanisms.
* **Define Analysis Tasks**: Creates a comprehensive list of different TMLE models to run, varying by:
    * **Nuisance Model Specfication**: Using a simple GLM or a SuperLearner with different numbers of algorithms.
    * **Survey Feature Awareness**:
        * **Non-Aware**: Ignores all survey design features.
        * **Partially-Aware**: Accounts for survey weights but ignores clustering and stratification.
        * **Fully-Aware**: Accounts for survey weights and uses the influence curve to calculate variance estimates that respect the complex survey design (clustering and stratification).
* **Parallel Execution**: Runs all the defined TMLE analyses in parallel to speed up computation time.
* **Save Results**: Saves the output of each TMLE analysis into a separate `.rds` file in the `results/` directory (e.g., `glm_weighted.rds`, `SL3_fully_aware.rds`).

### 3. `3_run_glm.R`

This script serves as a comparison to the TMLE analyses. It fits more traditional Generalized Linear Models (GLMs) to the imputed data to estimate the association between GDM and hypertension.

**Key Functions:**

* **Load Data**: Loads the `nhanes_singly_imputed_data.rds` file.
* **Fit Models**: Fits three types of GLMs to estimate different parameters:
    * **Odds Ratio (OR)**: Using a quasibinomial model with a logit link.
    * **Risk Ratio (RR)**: Using a quasipoisson model with a log link.
    * **Risk Difference / Average Treatment Effect (ATE)**: Using a Gaussian model with an identity link.
* **Vary Survey Awareness**: For each parameter (OR, RR, ATE), it fits the model in three ways:
    * **Non-Aware**: A standard `glm()` that ignores all survey features.
    * **Partially-Aware**: A `svyglm()` that only accounts for survey weights.
    * **Fully-Aware**: A `svyglm()` that accounts for weights, strata, and primary sampling units (PSUs).
* **Save Results**: Extracts the coefficients, confidence intervals, and p-values for the exposure (`A`) from each model and saves them to separate `.rds` files in the `base_results/` directory (e.g., `survey_glm_fully_aware.rds`).

---

## How to Run

1.  Place all three scripts in the same root directory.
2.  Run `1_save_analytic_data.R` first to generate the necessary `nhanes_singly_imputed_data.rds` file.
3.  Run `2_run_tmle.R` to perform the TMLE analyses. Results will be saved in the `results/` folder.
4.  Run `3_run_glm.R` to perform the GLM analyses. Results will be saved in the `base_results/` folder.
