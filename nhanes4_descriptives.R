# Load data
imputed_data_complete <- readRDS("E:/GitHub/survey-tmle/Real data analysis/nhanes_singly_imputed_data.rds")
# Ensure necessary libraries are loaded
library(gtsummary)
library(survey)
library(smd)
library(broom)

imputed_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = imputed_data_complete
)

# --- Table 1a with p-values and SMDs ---
table1a_smd <- tbl_svysummary(
  imputed_design,
  by = HYPERTENSION,
  include = -c(SEQN, SDMVPSU, SDMVSTRA, WTMEC2YR),
  label = list(
    RIDAGEYR ~ "Age (years)", 
    RIDRETH3 ~ "Race/Ethnicity",
    DMDEDUC2 ~ "Education", 
    INDFMPIR ~ "Poverty-Income Ratio",
    BMXBMI   ~ "Body Mass Index (kg/m^2)", 
    GDM_HISTORY ~ "GDM History",
    RHQ172 ~ "History of baby > 9 lbs"
  )
) %>%
  add_overall() %>%
  add_p() %>%
  add_difference(test = everything() ~ "smd") %>% 
  modify_header(
    list(
      p.value ~ "**p-value**",
      estimate ~ "**SMD**"
    )
  ) %>%
  gtsummary::modify_caption("Table 1a: Survey-Weighted Characteristics by Hypertension (Imputed Data with SMDs)")

# Print the new table
table1a_smd

# --- Table 1a using svyCreateTableOne ---
# Ensure necessary libraries are loaded
library(tableone)
library(survey)

# Your survey design object is the same
imputed_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = imputed_data_complete
)

# --- Create Table 1 using svyCreateTableOne ---

# 1. Define the list of variables you want in the table
my_vars <- c(
  "RIDAGEYR", "RIDRETH3", "DMDEDUC2", "INDFMPIR", 
  "BMXBMI", "GDM_HISTORY", "RHQ172"
)

# 2. Identify which of those variables are categorical
categorical_vars <- c(
  "RIDRETH3", "DMDEDUC2", "GDM_HISTORY", "RHQ172"
)

# 3. Generate the table object
table1_obj <- svyCreateTableOne(
  vars = my_vars,
  strata = "HYPERTENSION",
  data = imputed_design,
  factorVars = categorical_vars,
  addOverall = FALSE
)

# 4. Print the table and specify smd = TRUE to display the SMD column
print(table1_obj, smd = TRUE)


#---------------------------------------------------------------------
# --- Table 1b with p-values and SMDs ---


# Ensure necessary libraries are loaded
library(gtsummary)
library(survey)

# Assumes 'imputed_design' object is already created

# --- Table stratified by GDM_HISTORY ---
table_by_gdm <- tbl_svysummary(
  imputed_design,
  by = GDM_HISTORY, # Stratify by GDM_HISTORY
  include = -c(SEQN, SDMVPSU, SDMVSTRA, WTMEC2YR),
  label = list(
    RIDAGEYR ~ "Age (years)", 
    RIDRETH3 ~ "Race/Ethnicity",
    DMDEDUC2 ~ "Education", 
    INDFMPIR ~ "Poverty-Income Ratio",
    BMXBMI   ~ "Body Mass Index (kg/m^2)", 
    HYPERTENSION ~ "Hypertension", 
    RHQ172 ~ "History of baby > 9 lbs"
  )
) %>%
  add_overall() %>%
  add_p() %>%
  add_difference(test = everything() ~ "smd") %>% 
  modify_header(
    list(
      p.value ~ "**p-value**",
      estimate ~ "**SMD**"
    )
  ) %>%
  gtsummary::modify_caption("Table 1b: Survey-Weighted Characteristics by GDM History (Imputed Data with SMDs)")

# Print the new table
table_by_gdm

# --- Create Table 1 using svyCreateTableOne ---

# Ensure necessary libraries are loaded
library(tableone)
library(survey)

# --- Create Table 1 using svyCreateTableOne stratified by GDM_HISTORY ---

# 1. Define the list of variables for the table
my_vars_gdm <- c(
  "RIDAGEYR", "RIDRETH3", "DMDEDUC2", "INDFMPIR", 
  "BMXBMI", "HYPERTENSION", "RHQ172"
)

# 2. Identify which of those variables are categorical
categorical_vars_gdm <- c(
  "RIDRETH3", "DMDEDUC2", "HYPERTENSION", "RHQ172"
)

# 3. Generate the table object
table1_obj_gdm <- svyCreateTableOne(
  vars = my_vars_gdm,
  strata = "GDM_HISTORY", # Change strata to GDM_HISTORY
  data = imputed_design,
  factorVars = categorical_vars_gdm
)

# 4. Print the table with SMDs
print(table1_obj_gdm, smd = TRUE)