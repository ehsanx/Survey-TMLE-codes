# =================================================================
# Main Script to Run TMLE Analyses
# =================================================================

# --- 1. SETUP ---

# Load necessary libraries
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(tidyverse, here, survey, tmle, SuperLearner, earth, gam, glmnet, parallel)
require(dplyr)
require(here)
require(survey)
require(tmle)
require(SuperLearner)
require(earth)
require(gam)
require(glmnet)
require(parallel)

# Create a directory to store results if it doesn't exist
if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}

# Load the singly imputed dataset
imputed_data_complete <- readRDS(here("nhanes_singly_imputed_data.rds"))

# Prepare common data components for TMLE
Y <- as.numeric(imputed_data_complete$HYPERTENSION) - 1
A <- as.numeric(imputed_data_complete$GDM_HISTORY) - 1
W <- imputed_data_complete %>%
  select(RIDAGEYR, RIDRETH3, DMDEDUC2, INDFMPIR, BMXBMI, RHQ172) %>%
  mutate(across(where(is.factor), as.numeric))
svy_weights <- imputed_data_complete$WTMEC2YR


# --- 2. DEFINE THE WORKER FUNCTION ---

# This function will run a single TMLE analysis based on a configuration list
run_analysis <- function(config) {
  
  # Load libraries required by each parallel worker
  library(tmle)
  library(SuperLearner)
  library(here)
  library(survey)
  library(earth)
  library(gam)
  library(glmnet)
  library(MASS)
  
  # Print a message to track progress
  cat("Starting analysis:", config$name, "\n")
  
  # Run the specified TMLE analysis
  tmle_fit <- tmle(
    Y = Y,
    A = A,
    W = W,
    family = "binomial", # family="binomial" ensures RR and OR are calculated
    obsWeights = config$weights,
    Q.SL.library = config$library,
    g.SL.library = config$library
  )
  
  # --- Post-processing for survey variance ---
  # If the method requires complex survey variance, calculate it
  if (config$se_method == "complex") {
    
    # Extract all relevant influence curves (for ATE, log(RR), and log(OR))
    ic_ate <- tmle_fit$estimates$IC$IC.ATE
    ic_log_rr <- tmle_fit$estimates$IC$IC.logRR
    ic_log_or <- tmle_fit$estimates$IC$IC.logOR
    
    # Add all ICs to the data frame to be used in svydesign
    data_for_svy <- imputed_data_complete
    data_for_svy$ic_ate <- ic_ate
    data_for_svy$ic_log_rr <- ic_log_rr
    data_for_svy$ic_log_or <- ic_log_or
    
    # Create survey design object once
    svy_design_tmle <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~svy_weights, data = data_for_svy, nest = TRUE)
    
    # Calculate complex SE for each parameter
    se_ate_complex <- as.numeric(SE(svymean(~ic_ate, svy_design_tmle)))
    se_log_rr_complex <- as.numeric(SE(svymean(~ic_log_rr, svy_design_tmle)))
    se_log_or_complex <- as.numeric(SE(svymean(~ic_log_or, svy_design_tmle)))
    
    # Store the complex SEs in a named list
    se_complex_all <- list(
      ATE = se_ate_complex,
      logRR = se_log_rr_complex,
      logOR = se_log_or_complex
    )
    
    # Store the final result object
    result_to_save <- list(
      name = config$name,
      estimates = tmle_fit$estimates,
      se_complex = se_complex_all # Save the list of SEs
    )

  } else {
    # For other methods, just save the standard TMLE output
    result_to_save <- list(
      name = config$name,
      estimates = tmle_fit$estimates,
      se_complex = NULL # Mark that complex SE was not calculated
    )
  }
  
  # Save the result object to its own RDS file
  output_filename <- here("results", paste0(config$name, ".rds"))
  saveRDS(result_to_save, file = output_filename)
  
  cat("Finished and saved:", config$name, "\n")
  
  return(TRUE) # Return TRUE on success
}


# --- 3. DEFINE ANALYSIS CONFIGURATIONS ---

# Define Super Learner libraries
SL.library3 <- c("SL.glm", "SL.gam", "SL.earth")
SL.library4 <- c("SL.glm", "SL.gam", "SL.earth", "SL.glm.interaction")
SL.library5 <- c("SL.glm", "SL.gam", "SL.earth", "SL.glm.interaction", "SL.stepAIC")

SL.library8 <- c("SL.glm", "SL.gam", "SL.earth", "SL.glm.interaction", "SL.stepAIC", "SL.glmnet") 

libraries <- list(SL3 = SL.library3, SL4 = SL.library4, SL5 = SL.library5, SL8 = SL.library8)

# Create a list of all analysis configurations to run
analysis_tasks <- list(
  # GLM-based analyses
  list(name = "glm_unweighted",          weights = NULL,          library = "SL.glm", se_method = "tmle_default"),
  list(name = "glm_weighted",            weights = svy_weights,   library = "SL.glm", se_method = "tmle_default"),
  list(name = "glm_weighted_psu_strata", weights = svy_weights,   library = "SL.glm", se_method = "complex")
)

# Add Super Learner tasks
for (lib_name in names(libraries)) {
  
  current_lib <- libraries[[lib_name]]
  
  analysis_tasks <- c(analysis_tasks, list(
    # Non-Aware (no weights)
    list(name = paste0(lib_name, "_non_aware"), weights = NULL, library = current_lib, se_method = "tmle_default"),
    # Fully-Aware (survey weights + complex SE)
    list(name = paste0(lib_name, "_fully_aware"), weights = svy_weights, library = current_lib, se_method = "complex"),
    # Partially-Aware (survey weights, simple SE)
    list(name = paste0(lib_name, "_partially_aware"), weights = svy_weights, library = current_lib, se_method = "tmle_default")
  ))
}


# --- 4. RUN ANALYSES IN PARALLEL ---

# Detect number of cores
no_cores <- detectCores() - 1
if (no_cores < 1) no_cores <- 1

# Create cluster
cl <- makeCluster(no_cores)

# Export necessary objects and data to the cluster workers
clusterExport(cl, varlist = c("imputed_data_complete", "Y", "A", "W", "svy_weights", "run_analysis"))

# Run the analyses in parallel
cat("\n--- Starting all TMLE analyses in parallel ---\n")
parallel_results <- parLapply(cl, analysis_tasks, function(task) {
  # Wrap the function call in a try block to catch any errors
  tryCatch({
    run_analysis(task)
  }, error = function(e) {
    # Return the error message if something goes wrong
    return(paste("Error in", task$name, ":", e$message))
  })
})

# Stop the cluster
stopCluster(cl)

# --- 5. FINAL REPORT ---
cat("\n--- Parallel processing complete ---\n")
print("Results for each analysis have been saved to individual .rds files in the 'results' directory.")

# Check for any errors during the parallel run
errors <- parallel_results[sapply(parallel_results, is.character)]
if (length(errors) > 0) {
  cat("\nSome analyses failed with errors:\n")
  for (err in errors) {
    cat("-", err, "\n")
  }
} else {
  cat("\nAll analyses completed successfully.\n")
}