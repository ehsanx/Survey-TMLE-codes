# =====================================================================
# Nhanes/R/examples.R  —  example registry for the survey-TMLE NHANES application
#
# Single source of truth for the four candidate illustrations (see Nhanes/nhanes.md).
# Used by 00_download_cache.R (component union to fetch) and 01_build_analytic.R
# (per-example exposure/outcome/covariate derivations). All four are MEC-based, so
# all use WTMEC2YR / (number of pooled cycles).  NOTE: exact NHANES variable names/
# availability are confirmed against the cached raw tables in 01 (some items, e.g.
# the sleep-hours variable and the e-cigarette items, change name across cycles).
# =====================================================================

`%+%` <- function(a, b) paste0(a, b)   # string concat helper (used in `verify` notes below)

## six 2-year cycles 2007-2008 .. 2017-2018 -> table suffixes _E .. _J
CYCLES <- c(E = "2007-2008", F = "2009-2010", G = "2011-2012",
            H = "2013-2014", I = "2015-2016", J = "2017-2018")
CYCLE_SUFFIXES <- names(CYCLES)
N_CYCLES <- length(CYCLES)                      # = 6  -> pooled weight = WTMEC2YR / 6

## per-component table name. CDC renamed the BP exam table to BPXO in 2017-2018
## (oscillometric; BPXOSY1/BPXODI1). Add other per-cycle exceptions here as found.
tbl_name <- function(comp, suf) {
  if (comp == "BPX" && suf == "J") return("BPXO_J")   # 2017-18 oscillometric BP
  paste0(comp, "_", suf)
}

## Design variables carried for every example (from DEMO each cycle).
DESIGN_VARS <- c("SEQN", "SDMVSTRA", "SDMVPSU", "WTMEC2YR", "WTINT2YR", "SDDSRVYR")

## Shared core confounders (present/derivable in every example).
## Race: use RIDRETH1 (4-level, present ALL cycles); RIDRETH3 (adds Non-Hispanic
## Asian) only exists from 2011-12 (_G) on, so it is NOT cross-cycle comparable
## for the 6-cycle examples -> derive race from RIDRETH1 in 01.
CORE_CONFOUNDERS <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "RIDRETH3",
                      "DMDEDUC2", "INDFMPIR", "SDDSRVYR")

## ---------------------------------------------------------------------
## The four examples. `components` = NHANES components to merge for that
## example (DEMO is always included). `verify` flags items whose exact
## variable name/cycle-availability must be checked in 01 against raw tables.
## ---------------------------------------------------------------------
EXAMPLES <- list(

  E1 = list(
    id = "E1", label = "Short sleep -> obesity",
    cycles = CYCLE_SUFFIXES,                  # all 6 -> pooled weight = WTMEC2YR / 6
    components = c("DEMO", "BMX", "SLQ", "PAQ", "ALQ", "SMQ", "DPQ", "DR1TOT", "HIQ"),
    subpop = "adults aged >= 20 with non-missing sleep + BMI",
    exposure = list(desc = "short sleep (<7 h) vs healthy (7-9 h); long sleepers handled in sensitivity",
                    vars = c("SLD010H", "SLD012")),     # name changes ~2015-16
    outcome  = list(desc = "obesity, BMXBMI >= 30", vars = "BMXBMI"),
    covariates = c(CORE_CONFOUNDERS, "BMXBMI_NA",      # BMI is the OUTCOME here, not a covariate
                   "PAQ_MET", "DR1TKCAL", "ALQ", "SMQ020", "SMQ040", "DPQ_PHQ9", "DMDMARTL", "HIQ011"),
    verify = "Sleep-hours variable: SLD010H (2007-2014) vs SLD012 (2015-2018). Depression (DPQ) here is a "
           %+% "CONFOUNDER of sleep->obesity, not a mediator."
  ),

  E2 = list(
    id = "E2", label = "Food insecurity -> depression",
    cycles = CYCLE_SUFFIXES,                  # all 6 -> pooled weight = WTMEC2YR / 6
    components = c("DEMO", "FSQ", "DPQ", "BMX", "MCQ", "HIQ", "SMQ"),
    subpop = "adults aged >= 20 with non-missing food-security + PHQ-9",
    exposure = list(desc = "food insecure (HFSSM >= 3 affirmatives) vs food secure",
                    vars = c("FSDAD", "FSDHH", "FSD032A")),   # adult/household scales
    outcome  = list(desc = "depression, PHQ-9 (DPQ010-090) summed >= 10", vars = paste0("DPQ", sprintf("%03d", seq(10, 90, 10)))),
    covariates = c(CORE_CONFOUNDERS, "DMDMARTL", "DMDHHSIZ", "HIQ011", "BMXBMI", "MCQ_chron"),
    verify = "Food-security scale variable name (FSDAD adult vs FSDHH household) and cut-point; confirm "
           %+% "DPQ present + scoring (>=10)."
  ),

  E3 = list(
    id = "E3", label = "E-cigarette use -> hypertension",
    cycles = c("I", "J"),                     # e-cig items SMQ900/905 only 2015-18 -> pooled weight = WTMEC2YR / 2
    components = c("DEMO", "SMQ", "SMQRTU", "BPX", "BPQ", "BMX", "ALQ", "PAQ", "DIQ", "MCQ"),
    subpop = "adults aged >= 20 with non-missing e-cig + BP",
    exposure = list(desc = "ever-vaped vs never (primary); current (past-30-day) vs never (sensitivity)",
                    vars = c("SMQ900", "SMQ905", "ECQ")),     # e-cig items; verify per cycle
    outcome  = list(desc = "hypertension: mean SBP>=130 or DBP>=80 OR BPQ020 (Dx) OR BPQ040A (meds)",
                    vars = c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4", "BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4",
                             "BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXODI1", "BPXODI2", "BPXODI3",
                             "BPQ020", "BPQ040A")),
    covariates = c(CORE_CONFOUNDERS, "SMQ020", "SMQ040", "SMD650", "BMXBMI", "ALQ", "PAQ_MET", "DIQ010"),
    verify = "E-cig items appear mainly 2015-2018 (SMQRTU/ECQ); confirm var names per cycle. Combustible "
           %+% "smoking (SMQ020/040/SMD650) is the key confounder. DIQ010 is a PRE-exposure confounder here."
  ),

  E4 = list(
    id = "E4", label = "GDM history -> hypertension (hardened)",
    cycles = CYCLE_SUFFIXES,                  # all 6 -> pooled weight = WTMEC2YR / 6
    components = c("DEMO", "RHQ", "DIQ", "BPX", "BPQ", "BMX", "SMQ", "MCQ"),
    subpop = "ever-pregnant (RHQ131==Yes), non-diabetic (DIQ010!=Yes) women aged >= 20",
    exposure = list(desc = "history of gestational diabetes (RHQ162==Yes) vs not", vars = "RHQ162"),
    outcome  = list(desc = "hypertension: mean SBP>=130 or DBP>=80 OR BPQ020 OR BPQ040A (hardened)",
                    vars = c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4", "BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4",
                             "BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXODI1", "BPXODI2", "BPXODI3",
                             "BPQ020", "BPQ040A")),
    covariates = c(CORE_CONFOUNDERS, "RHQ160", "RHQ171", "RHD180", "RHQ031", "MCQ300c", "SMQ020", "SMQ040", "BMXBMI"),
    # BMXBMI is a PROXY for unmeasured pre-pregnancy adiposity (partial-mediator caveat; sensitivity).
    exclude = c("DIQ010 (subsequent diabetes = MEDIATOR; use only as cohort exclusion)",
                "post-pregnancy weight gain", "antihypertensive meds (-> outcome def)",
                "prevalent CVD (MCQ160* = collider)", "fasting glucose/HbA1c/lipids (mediators)",
                "same-pregnancy preeclampsia (collider)"),
    verify = "RHQ162 (told had diabetes during pregnancy) availability per cycle; reproductive vars "
           %+% "RHQ160/RHQ171/RHD180/RHQ031. Exclusion: drop prevalent diabetes (DIQ010==Yes)."
  )
)

## convenience: the UNION of components to download/cache (DEMO always first)
COMPONENTS_ALL <- unique(c("DEMO", unlist(lapply(EXAMPLES, `[[`, "components"))))

invisible(NULL)
