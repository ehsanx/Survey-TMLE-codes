# --------------------------------------------------------------------------- #
# SCRIPT START: Generate Singly Imputed NHANES Analytic Dataset
# --------------------------------------------------------------------------- #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PART 1: LOAD PACKAGES and PREPARE RAW DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# -- Step 1.1: Install and load all necessary packages --
# Install pacman if you don't have it, a package manager for R
if (!require("pacman")) install.packages("pacman")

# Use pacman to load all necessary packages from both scripts
pacman::p_load(tidyverse, nhanesA, here, survey, gtsummary, naniar, DataExplorer, mice)


# -- Step 1.2: Download, Merge, and Select Data (2007-2018) --
# Define the NHANES data components
components <- c("DEMO", "DIQ", "RHQ", "BMX", "BPX")
cycles <- c("E", "F", "G", "H", "I", "J")
cycle_years <- c("2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018")

# Define the specific variables we need
vars_to_keep <- c(
  "SEQN", "RIAGENDR", "RIDAGEYR", 
  "RIDRETH1", "RIDRETH3", # Both race/ethnicity variables
  "DMDEDUC2", "INDFMPIR",
  "DIQ010", "RHQ131", "RHQ162", "RHQ172", "BMXBMI", "BPXSY1", "BPXDI1",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA"
)

cycle_map <- setNames(cycle_years, cycles)
all_cycles_data <- list()

for (cycle_letter in cycles) {
  cycle_suffix <- paste0("_", cycle_letter)
  tables_for_cycle <- map(components, ~ nhanes(paste0(.x, cycle_suffix)))
  merged_cycle_data <- tables_for_cycle %>% reduce(full_join, by = "SEQN")
  selected_data <- merged_cycle_data %>%
    mutate(CYCLE = cycle_map[cycle_letter]) %>%
    select(any_of(c("SEQN", "CYCLE", vars_to_keep)))
  all_cycles_data[[cycle_letter]] <- selected_data
}

# Combine all survey cycles into one raw dataframe
nhanes_raw <- bind_rows(all_cycles_data)


# -- Step 1.3: Create Analysis Variables and Eligibility Indicator --
nhanes_final <- nhanes_raw %>%
  mutate(
    # Combine RIDRETH1 and RIDRETH3 into a single, consistent variable
    RACE_ETHNICITY = coalesce(RIDRETH1, RIDRETH3)
  ) %>%
  mutate(
    # Create an eligibility indicator based on study criteria
    is_eligible = (
      RIAGENDR == "Female" &
        RIDAGEYR >= 20 &
        !is.na(RHQ162) &  # Proxy for ever-pregnant
        DIQ010 != "Yes"      # Exclude those with pre-existing diabetes
    ),
    
    # Outcome: Hypertension
    HYPERTENSION = case_when(
      BPXSY1 >= 140 | BPXDI1 >= 90 ~ "Yes",
      BPXSY1 < 140 & BPXDI1 < 90 ~ "No",
      TRUE ~ NA_character_
    ),
    
    # Exposure: GDM History
    GDM_HISTORY = case_when(
      RHQ162 == "Yes" ~ "Yes",
      RHQ162 == "No" ~ "No",
      TRUE ~ NA_character_
    )
  ) %>%
  # Recode categorical variables into simpler groups
  mutate(
    RACE_ETHNICITY = case_when(
      RACE_ETHNICITY %in% c("Mexican American", "Other Hispanic") ~ "Hispanic",
      RACE_ETHNICITY == "Non-Hispanic White" ~ "Non-Hispanic White",
      RACE_ETHNICITY == "Non-Hispanic Black" ~ "Non-Hispanic Black",
      RACE_ETHNICITY %in% c("Non-Hispanic Asian", "Other Race - Including Multi-Racial") ~ "Other",
      TRUE ~ NA_character_
    ),
    DMDEDUC2 = case_when(
      DMDEDUC2 %in% c("Less Than 9th Grade", "9-11th Grade (Includes 12th grade with no diploma)",
                      "Less than 9th grade", "9-11th grade (Includes 12th grade with no diploma)") ~ "Less than High School",
      DMDEDUC2 %in% c("High School Grad/GED or Equivalent", "High school graduate/GED or equivalent") ~ "High School Graduate / GED",
      DMDEDUC2 %in% c("Some College or AA degree", "Some college or AA degree", "Some College") ~ "Some College",
      DMDEDUC2 %in% c("College Graduate or above", "College graduate or above") ~ "College Graduate or above",
      TRUE ~ NA_character_
    ),
    RHQ172 = case_when(
      RHQ172 == "Yes" ~ "Yes",
      RHQ172 == "No" ~ "No",
      TRUE ~ NA_character_
    )
  ) %>%
  # Convert cleaned character vectors into factors
  mutate(
    HYPERTENSION = factor(HYPERTENSION, levels = c("No", "Yes")),
    GDM_HISTORY = factor(GDM_HISTORY, levels = c("No", "Yes")),
    RIDRETH3 = factor(RACE_ETHNICITY, levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other")),
    DMDEDUC2 = factor(DMDEDUC2, levels = c("Less than High School", "High School Graduate / GED", "Some College", "College Graduate or above")),
    RHQ172 = factor(RHQ172, levels = c("No", "Yes"))
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PART 2: PERFORM IMPUTATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# -- Step 2.1: Prepare Analytic Dataset for Imputation --
# Create the analytic sample: eligible participants with non-missing exposure and outcome.
analytic_sample <- subset(nhanes_final,
                          is_eligible == TRUE &
                            !is.na(GDM_HISTORY) &
                            !is.na(HYPERTENSION))

# Select only the variables to be used in the imputation and subsequent analysis
imputation_vars <- c(
  "SEQN", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", # Survey vars
  "GDM_HISTORY", "HYPERTENSION",           # Outcome/Exposure
  "RIDAGEYR", "RIDRETH3", "DMDEDUC2",       # Covariates
  "INDFMPIR", "BMXBMI", "RHQ172"
)

dat_to_impute <- analytic_sample %>%
  select(any_of(imputation_vars))


# -- Step 2.2: Perform Single Imputation using `mice` --
# Pre-process factors to remove any unused levels which can cause mice errors
dat_to_impute_processed <- dat_to_impute %>%
  mutate(across(where(is.factor), droplevels))

# Initialize mice to get the predictor matrix and method vector
ini <- mice(dat_to_impute_processed, maxit = 0, print = FALSE)
pred <- ini$predictorMatrix
meth <- ini$method

# Do not use ID and survey weight variables as predictors
pred[, c("SEQN", "WTMEC2YR", "SDMVPSU")] <- 0

# Set imputation methods for variables that have missing data
meth["BMXBMI"] <- "pmm"
meth["INDFMPIR"] <- "pmm"
meth["DMDEDUC2"] <- "polyreg"
meth["RHQ172"] <- "polyreg" 

# Run the imputation
imputation <- mice(
  data = dat_to_impute_processed,
  seed = 123,
  predictorMatrix = pred,
  method = meth,
  m = 1,      # m=1 for a single imputation
  maxit = 5,  # 5 iterations is usually sufficient
  print = TRUE
)

# Create the final, complete dataset from the single imputation
imputed_data_complete <- mice::complete(imputation, action = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PART 3: SAVE THE FINAL IMPUTED DATASET
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# -- Step 3.1: Save the Singly Imputed Dataset --
# Create a directory if it doesn't exist to prevent errors
if (!dir.exists(here("Real data analysis"))) {
  dir.create(here("Real data analysis"), recursive = TRUE)
}

# Define the output file path
output_file <- here("Real data analysis/nhanes_singly_imputed_data.rds")

# Save the final data frame to an RDS file
saveRDS(imputed_data_complete, file = output_file)

# --- Final Success Message ---
print("---------------------------------------------------------------------")
print(paste("SUCCESS: The final singly imputed dataset has been saved to:", output_file))
print(paste("Final dataset contains", nrow(imputed_data_complete), "rows and", ncol(imputed_data_complete), "columns."))
print("---------------------------------------------------------------------")


# --------------------------------------------------------------------------- #
# SCRIPT END
# --------------------------------------------------------------------------- #