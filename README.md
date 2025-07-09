# A Practical Guide to Valid Causal Inference with TMLE in Complex Surveys

This repository contains the R code to reproduce the analyses and simulations presented in the manuscript: *"A Practical Guide to Valid Causal Inference with TMLE in Complex Surveys: A Survey Design-Aware Variance Estimator and Application"*.

## Project Overview

Complex population-based surveys are a cornerstone of public health research, yet applying modern causal inference methods such as Targeted Maximum Likelihood Estimation (TMLE) to these datasets is challenging. The absence of a formal framework for incorporating survey design features—namely weighting, stratification, and clustering—can lead to invalid statistical inference.

This project introduces and validates a practical and theoretically-grounded data analysis strategy to address this gap. The work is divided into two main parts:

1. **Real-Data Analysis**: We apply the proposed "Fully-Aware" TMLE methodology to a real-world research question using data from the National Health and Nutrition Examination Survey (NHANES).

2. **Simulation Study**: We conduct an extensive simulation study that mimics the complexity of NHANES to quantify the performance of our proposed method against naïve approaches that ignore survey design features.

The goal is to provide a robust, replicable framework for applied researchers to conduct valid causal inference using TMLE with complex survey data.

---

## Repository Structure

The repository is organized into two main components: the real-data analysis and the simulation study.

* `/`: Contains the primary R scripts.

* `/base_results/`: Stores the output from the GLM-based conditional models in the real-data analysis (shown in the app, if selected).

* `/results/`: Stores the output from the TMLE-based models in the real-data analysis.

* `/data/`: Stores the output from the simulation studies.

---

## Part 1: Real-Data Analysis (NHANES)

This part of the project investigates the association between a history of gestational diabetes mellitus (GDM) and the subsequent risk of developing hypertension in adult women in the United States, using six cycles of NHANES data (2007-2018).

### Real-Data Scripts

* `1_save_analytic_data.R`: The foundational script. It downloads, merges, and cleans the raw NHANES data. It then creates the necessary analysis variables, defines the eligible study population, performs a single imputation for missing covariate data using `mice`, and saves the final analytic dataset (`nhanes_singly_imputed_data.rds`). **This script must be run first.**

* `2_run_tmle.R`: Loads the imputed data and runs a series of TMLE analyses. It uses different `SuperLearner` libraries to estimate the causal effect and compares three approaches: "Non-Aware" (ignores all survey features), "Partially-Aware" (uses weights only), and "Fully-Aware" (uses weights, strata, and clusters for variance estimation). Results are saved to the `/results` directory.

* `3_run_glm.R`: Serves as a comparison to the TMLE analyses. It fits more traditional survey-weighted Generalized Linear Models (GLMs) to estimate odds ratios, risk ratios, and risk differences, also comparing the "Non-Aware", "Partially-Aware", and "Fully-Aware" approaches. Results are saved to the `/base_results` directory.

### How to Run the Real-Data Analysis

1. Place `1_save_analytic_data.R`, `2_run_tmle.R`, and `3_run_glm.R` in the same root directory.

2. Run `1_save_analytic_data.R` to generate the `nhanes_singly_imputed_data.rds` file.

3. Run `2_run_tmle.R` to perform the TMLE analyses.

4. Run `3_run_glm.R` to perform the comparative GLM analyses.

---

## Part 2: Simulation Study

This part of the project validates the proposed methodology through a Monte Carlo simulation study. The simulation generates data that mimics the key features of NHANES, including stratification, clustering, and informative sampling.

### Simulation Scenarios

The simulations test estimator performance under two scenarios:

1. **Simple Scenario**: The relationships between confounders, treatment, and outcome are linear. Parametric models are correctly specified.

2. **Complex Scenario**: The relationships are highly non-linear, designed to challenge the estimators and test the flexibility of machine learning approaches under model misspecification.

### Simulation Scripts

The simulation is broken into multiple scripts, each testing a different set of models. They all follow a similar structure: a data generation function (`create.data`), a main simulation function (`run_one_simulation`), and a parallel backend to run 1,000 repetitions for both the simple and complex scenarios.

* `sim1glm.R`: Runs a TMLE analysis where the nuisance models (Q and g) are estimated using a standard `glm`.

* `sim1svyglm.R`: Runs a parametric g-computation analysis using `svyglm` as a benchmark.

* `sim3sl.R`, `sim4sl.R`, `sim5sl.R`, `sim8sl.R`: These scripts run TMLE using `SuperLearner` with increasingly complex libraries of machine learning algorithms. The number in the filename (3, 4, 5, 8) corresponds to the number of learners in the library, as detailed in the manuscript. Each script also compares the "Non-Aware", "Partially-Aware", "Fully-Aware", and a special "Fully-Aware-CV" (with survey-aware cross-validation) approach.

### How to Run the Simulations

Each simulation script is self-contained and can be run independently. For example, to run the simulation for the Super Learner with 3 algorithms, simply execute `sim3sl.R`. The script will create the necessary directories and save the aggregated results as an `.rds` file in the `/data` directory (e.g., `simulation_results_3.rds`).

---

## Dependencies

This project requires the following R packages. You can install them using `install.packages()`.

* `pacman`
* `tidyverse`
* `nhanesA`
* `here`
* `survey`
* `gtsummary`
* `naniar`
* `DataExplorer`
* `mice`
* `tmle`
* `SuperLearner`
* `earth`
* `gam`
* `glmnet`
* `parallel`
* `rsimsum`
* `surveyCV`
