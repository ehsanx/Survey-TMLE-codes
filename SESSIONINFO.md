# Reproducibility — software environment

The analyses were produced in two environments, and both reproduce the number audit
(`Rscript verify_numbers.R` → `460 checks, 460 pass, 0 FAIL`):

- **Cluster (as run):** the simulation study and the NHANES enhancement/sensitivity
  analyses were run on a compute cluster under **R 4.4.0**. The package versions
  below are taken from the run manifests recorded at execution time; the exact
  per-run versions are stored alongside each run's outputs.
- **Local:** figures and tables are regenerated, and the number audit is run, under
  **R 4.5.1**.

Minor package-version differences between the two environments do not change the
quoted results — every number is re-derived from the committed summary CSVs by
`verify_numbers.R`.

## Cluster environment (as run — produced the paper's numbers)

`R version 4.4.0 (2024-04-24)`

| Package | Version |
|---|---|
| SuperLearner | 2.0.40 |
| tmle | 2.1.1 |
| survey | 4.5 |
| surveyCV | 0.2.0 |
| earth | 5.3.5 |
| gam | 1.22.7 |
| glmnet | 5.0 |
| ranger | 0.18.0 |

(The NHANES build/imputation/download packages — `mice`, `dplyr`, `tidyr`, `purrr`,
`nhanesA` — are recorded in the NHANES run manifests; local versions are listed below.)

## Local environment (regenerates exhibits + runs the audit)

`R version 4.5.1 (2025-06-13 ucrt)`

| Package | Version |
|---|---|
| SuperLearner | 2.0.29 |
| tmle | 2.1.1 |
| survey | 4.4.8 |
| surveyCV | 0.2.0 |
| earth | 5.3.4 |
| gam | 1.22.6 |
| glmnet | 4.1.10 |
| ranger | 0.18.0 |
| mice | 3.18.0 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| purrr | 1.1.0 |
| nhanesA | 1.4 |
| rsimsum | 0.13.0 |

Install the current versions with `Rscript install_packages.R`. To regenerate this
file, run `sessionInfo()` (or `utils::packageVersion()` over the package list) in your
own environment.
