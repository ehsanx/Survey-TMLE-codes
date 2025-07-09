<div align="center">
  <h1 align="center">A Practical Guide to Valid Causal Inference with TMLE in Complex Surveys</h1>
  <p align="center">
    <em>Code to reproduce the analyses and simulations for the manuscript: "A Practical Guide to Valid Causal Inference with TMLE in Complex Surveys: A Survey Design-Aware Variance Estimator and Application".</em>
  </p>
</div>

---

### 📖 **Project Overview**

Complex population-based surveys are a cornerstone of public health research, yet applying modern causal inference methods such as **Targeted Maximum Likelihood Estimation (TMLE)** to these datasets is challenging. The absence of a formal framework for incorporating survey design features—namely weighting, stratification, and clustering—can lead to invalid statistical inference.

This project introduces and validates a practical and theoretically-grounded data analysis strategy to address this gap. The work is divided into two main parts:

1.  **Real-Data Analysis**: We apply the proposed "Fully-Aware" TMLE methodology to a real-world research question using data from the National Health and Nutrition Examination Survey (NHANES).
2.  **Simulation Study**: We conduct an extensive simulation study that mimics the complexity of NHANES to quantify the performance of our proposed method against naïve approaches that ignore survey design features.

The goal is to provide a robust, replicable framework for applied researchers to conduct valid causal inference using TMLE with complex survey data.

---

### ✨ **Key Features**

| Feature                 | Description                                                                                                                              |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| 📊 **Real-Data Analysis** | Investigates the link between gestational diabetes and hypertension using NHANES data (2007-2018).                                       |
| ⚙️ **Simulation Engine** | Validates the methodology with a Monte Carlo simulation mimicking complex survey designs (stratification, clustering, informative sampling). |
| 🤖 **Machine Learning** | Implements TMLE with `SuperLearner` to allow for flexible, data-adaptive modeling of nuisance parameters.                                  |
| ⚖️ **Method Comparison** | Compares "Fully-Aware", "Partially-Aware", and "Non-Aware" approaches for both TMLE and traditional GLMs.                                  |
|  Reproducible         | All scripts are provided to ensure full reproducibility of the manuscript's findings.                                                    |

---

### 📂 **Repository Structure**

The repository is organized into the two main project components, with dedicated folders for scripts and results.

```
.
├── 📜 nhanes1_save_analytic_data.R
├── 📜 nhanes2_run_tmle.R
├── 📜 nhanes3_run_glm.R
├── 📜 sim1glm.R
├── 📜 sim1svyglm.R
├── 📜 sim3sl.R
├── 📜 sim4sl.R
├── 📜 sim5sl.R
├── 📜 sim8sl.R
├── 📁 base_results/
│   └── (Outputs from conditional GLM models)
├── 📁 results/
│   └── (Outputs from marginal TMLE models)
└── 📁 data/
    └── (Outputs from all simulation studies)
```

---

### 🚀 **How to Run the Analyses**

#### **Part 1: Real-Data Analysis (NHANES)**

This analysis investigates the association between gestational diabetes (GDM) and hypertension using six cycles of NHANES data.

**Scripts:**

* `nhanes1_save_analytic_data.R`: **(Run First)** Prepares the analytic dataset.
* `nhanes2_run_tmle.R`: Runs the primary TMLE analyses.
* `nhanes3_run_glm.R`: Runs comparative GLM analyses.

**Execution Order:**

1.  Ensure all three `nhanes*.R` scripts are in the same directory.
2.  Execute the scripts in numerical order:

    ```bash
    # Step 1: Create the imputed analytic dataset
    Rscript nhanes1_save_analytic_data.R

    # Step 2: Run the TMLE analyses
    Rscript nhanes2_run_tmle.R

    # Step 3: Run the comparative GLM analyses
    Rscript nhanes3_run_glm.R
    ```

#### **Part 2: Simulation Study**

This study validates the methodology using simulated data mimicking NHANES.

**Scenarios:**

* **Simple**: Nuisance models are linear and correctly specified.
* **Complex**: Nuisance models are non-linear, testing robustness to misspecification.

**Scripts:**

The simulation is broken into multiple scripts, each testing a different modeling strategy (`glm` (conditional), `svyglm` (g-computation), or `SuperLearner` with 3, 4, 5, or 8 algorithms).

**Execution:**

Each simulation script is self-contained and can be run independently. For example, to run the simulation for the Super Learner with 3 algorithms:

```bash
Rscript sim3sl.R
```

The script will create the necessary directories and save the aggregated results as an `.rds` file in the `/data` directory (e.g., `simulation_results_3.rds`).

---

### 📦 **Dependencies**

This project requires the following R packages. You can install them all at once using `pacman`.

```R
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, nhanesA, here, survey, gtsummary, naniar, DataExplorer,
  mice, tmle, SuperLearner, earth, gam, glmnet, parallel, rsimsum, surveyCV
)
```

---

### 🎓 **Citation**

If you use the code or methods from this repository in your work, please cite the original manuscript.

> Karim, ME. (2025). A Practical Guide to Valid Causal Inference with TMLE in Complex Surveys: A Survey Design-Aware Variance Estimator and Application. *submitted*.


---

### 📬 **Contact**

For questions, collaborations, or to report any issues, please feel free to contact the author:

- M Ehsan Karim
- 🌐 *Website*: [ehsank.com](https://ehsank.com/)
- 📧 *Email*: [ehsan.karim@ubc.ca](mailto:ehsan.karim@ubc.ca)