# =====================================================================
# Nhanes/R/01_build_analytic.R  —  per-example analytic datasets from cached raw
#
# For each example E1..E4 (Nhanes/R/examples.R) this:
#   (1) reads the cached raw tables for that example's cycles and merges a curated
#       column set within cycle (by SEQN), then binds cycles;
#   (2) harmonizes cross-cycle variables (BP BPXSY*->BPXOSY*, sleep SLD010H->SLD012,
#       race RIDRETH1);
#   (3) decodes integer codes (translated=FALSE) with sentinel (7/9/77/99) -> NA;
#   (4) derives the binary exposure A, binary outcome Y, and a JUSTIFIED covariate set;
#   (5) builds the pooled MEC weight = WTMEC2YR / (number of cycles used) and the
#       sub-population indicator `inpop` (eligibility & observed A,Y) WITHOUT deleting
#       rows -- the FULL MEC sample is retained so the design keeps all PSUs/strata
#       (subset the DESIGN in 03, per EpiMethods surveydata8).
# Saves Nhanes/analytic/<id>_analytic.rds and prints a per-example summary.
#
# Covariate / eligibility rationale is documented inline per example and in
# Nhanes/nhanes.md (DAG: condition only on PRE-exposure common causes; never on
# descendants of exposure/outcome; BMI in E4 is a pre-pregnancy-adiposity PROXY).
# =====================================================================

suppressMessages({ library(dplyr); library(tidyr); library(purrr) })

REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
source(file.path(REPO_ROOT, "Nhanes", "R", "examples.R"))
raw_dir <- file.path(REPO_ROOT, "Nhanes", "raw")
out_dir <- file.path(REPO_ROOT, "Nhanes", "analytic"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rd <- function(tbl) { f <- file.path(raw_dir, paste0(tbl, ".rds")); if (file.exists(f)) readRDS(f) else NULL }
na_if_in <- function(x, codes) { x[x %in% codes] <- NA; x }     # map sentinel codes to NA
pick <- function(df, cols) df[, intersect(cols, names(df)), drop = FALSE]   # safe column select

## curated raw columns to keep per component (union across examples; missing ones ignored)
COLS <- list(
  DEMO = c("SEQN","SDMVSTRA","SDMVPSU","WTMEC2YR","SDDSRVYR","RIAGENDR","RIDAGEYR",
           "RIDRETH1","DMDEDUC2","INDFMPIR","DMDMARTL","DMDHHSIZ"),
  BMX  = c("SEQN","BMXBMI"),
  BPX  = c("SEQN","BPXSY1","BPXSY2","BPXSY3","BPXSY4","BPXDI1","BPXDI2","BPXDI3","BPXDI4",
           "BPXOSY1","BPXOSY2","BPXOSY3","BPXODI1","BPXODI2","BPXODI3"),
  BPQ  = c("SEQN","BPQ020","BPQ040A"),
  RHQ  = c("SEQN","RHQ131","RHQ162","RHQ160","RHQ171","RHD180","RHQ031"),
  DIQ  = c("SEQN","DIQ010"),
  SLQ  = c("SEQN","SLD010H","SLD012"),
  FSQ  = c("SEQN","FSDAD"),
  DPQ  = c("SEQN", paste0("DPQ", sprintf("%03d", seq(10, 90, 10)))),
  SMQ  = c("SEQN","SMQ020","SMQ040","SMD650","SMQ900","SMQ905"),
  PAQ  = c("SEQN","PAQ650","PAQ665"),
  MCQ  = c("SEQN","MCQ300C","MCQ300c","MCQ300A","MCQ300a"),
  HIQ  = c("SEQN","HIQ011")
)

## ---- merge one example's cached raw into a full MEC frame (curated cols) -------
merge_example <- function(ex) {
  comps <- unique(c("DEMO", ex$components))
  per_cycle <- lapply(ex$cycles, function(suf) {
    tabs <- lapply(comps, function(comp) {
      x <- rd(tbl_name(comp, suf)); if (is.null(x)) return(NULL)
      pick(x, COLS[[comp]] %||% names(x))
    })
    tabs <- Filter(Negate(is.null), tabs)
    df <- Reduce(function(a, b) full_join(a, b, by = "SEQN"), tabs)
    df$CYCLE_SUF <- suf
    df
  })
  bind_rows(per_cycle)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- shared harmonizers -------------------------------------------------------
harmonize_bp <- function(df) {
  sy <- intersect(c("BPXSY1","BPXSY2","BPXSY3","BPXSY4","BPXOSY1","BPXOSY2","BPXOSY3"), names(df))
  di <- intersect(c("BPXDI1","BPXDI2","BPXDI3","BPXDI4","BPXODI1","BPXODI2","BPXODI3"), names(df))
  zna <- function(M) { M[M == 0] <- NA; M }                       # a 0 BP reading = missing
  df$SBP <- if (length(sy)) rowMeans(zna(as.matrix(df[, sy, drop = FALSE])), na.rm = TRUE) else NA_real_
  df$DBP <- if (length(di)) rowMeans(zna(as.matrix(df[, di, drop = FALSE])), na.rm = TRUE) else NA_real_
  df$SBP[is.nan(df$SBP)] <- NA; df$DBP[is.nan(df$DBP)] <- NA
  df
}
race4 <- function(RIDRETH1)                                        # 4-level, all cycles
  factor(c("MexAm","OtherHisp","NHWhite","NHBlack","Other")[RIDRETH1],
         levels = c("NHWhite","NHBlack","MexAm","OtherHisp","Other"))
educ3 <- function(DMDEDUC2) {                                      # collapse to 3 levels
  e <- na_if_in(DMDEDUC2, c(7, 9))
  factor(ifelse(e <= 2, "<HS", ifelse(e == 3, "HS", "SomeColl+")),
         levels = c("<HS","HS","SomeColl+"))
}
smk3 <- function(SMQ020, SMQ040) {                                 # never/former/current
  ever <- na_if_in(SMQ020, c(7, 9))
  cur  <- na_if_in(SMQ040, c(7, 9))
  out <- rep(NA_character_, length(ever))
  out[ever == 2] <- "Never"
  out[ever == 1 & cur == 3] <- "Former"
  out[ever == 1 & cur %in% c(1, 2)] <- "Current"
  factor(out, levels = c("Never","Former","Current"))
}
married2 <- function(DMDMARTL) {
  m <- na_if_in(DMDMARTL, c(77, 99))
  factor(ifelse(m %in% c(1, 6), "Married/partner", "Not"), levels = c("Not","Married/partner"))
}
htn_outcome <- function(df) {                                     # hardened HTN
  bpq020 <- na_if_in(df$BPQ020, c(7, 9)); bpq040 <- na_if_in(df$BPQ040A, c(7, 9))
  hi_sbp <- df$SBP >= 130; hi_dbp <- df$DBP >= 80
  # %in% coerces NA -> FALSE per component so a missing meds answer does not blank Y
  y <- as.integer((hi_sbp %in% TRUE) | (hi_dbp %in% TRUE) | (bpq020 %in% 1) | (bpq040 %in% 1))
  # define only where there is at least a BP reading OR a self-report answer
  has_info <- (!is.na(df$SBP) | !is.na(df$DBP) | !is.na(bpq020))
  y[!has_info] <- NA
  y
}
phq9 <- function(df) {                                            # PHQ-9 >= 10, complete 9 items
  items <- paste0("DPQ", sprintf("%03d", seq(10, 90, 10)))
  M <- sapply(items, function(c) na_if_in(df[[c]], c(7, 9)))
  ok <- rowSums(is.na(M)) == 0
  s  <- rowSums(M); y <- as.integer(s >= 10); y[!ok] <- NA; y
}

## ============================ per-example derivations =========================
## Each returns the FULL MEC frame with: A, Y, covariate columns, design vars,
## CYCLE_SUF, and `inpop` (eligibility & !is.na(A) & !is.na(Y)). No rows dropped.

derive_E1 <- function(df) {                       # short sleep -> obesity (adults 20+)
  slp <- ifelse(is.na(df$SLD012), na_if_in(df$SLD010H, c(77, 99)), df$SLD012)  # hours
  df$A <- ifelse(slp < 7, 1L, ifelse(slp >= 7 & slp <= 9, 0L, NA_integer_))    # long(>9)->NA (primary)
  df$Y <- as.integer(df$BMXBMI >= 30)
  df$age <- df$RIDAGEYR; df$sex <- factor(df$RIAGENDR, 1:2, c("Male","Female"))
  df$race <- race4(df$RIDRETH1); df$educ <- educ3(df$DMDEDUC2)
  df$pir <- df$INDFMPIR; df$cycle <- factor(df$SDDSRVYR)
  df$smk <- smk3(df$SMQ020, df$SMQ040); df$married <- married2(df$DMDMARTL)
  df$pa <- factor(ifelse(na_if_in(df$PAQ650, c(7,9)) == 1 | na_if_in(df$PAQ665, c(7,9)) == 1, "Yes","No"))
  df$phq <- phq9(df)                              # depression = CONFOUNDER of sleep->obesity here
  df$inpop <- with(df, RIDAGEYR >= 20 & !is.na(A) & !is.na(Y))
  attr(df, "covs") <- c("age","sex","race","educ","pir","cycle","smk","married","pa","phq")
  df
}

derive_E2 <- function(df) {                       # food insecurity -> depression (adults 20+)
  fs <- na_if_in(df$FSDAD, c(7, 9))               # 1 full,2 marginal,3 low,4 very low
  df$A <- ifelse(fs >= 3, 1L, ifelse(fs %in% c(1,2), 0L, NA_integer_))
  df$Y <- phq9(df)
  df$age <- df$RIDAGEYR; df$sex <- factor(df$RIAGENDR, 1:2, c("Male","Female"))
  df$race <- race4(df$RIDRETH1); df$educ <- educ3(df$DMDEDUC2)
  df$pir <- df$INDFMPIR; df$cycle <- factor(df$SDDSRVYR)
  df$married <- married2(df$DMDMARTL); df$hhsize <- df$DMDHHSIZ
  df$insured <- factor(ifelse(na_if_in(df$HIQ011, c(7,9)) == 1, "Yes","No"))
  df$bmi <- df$BMXBMI; df$smk <- smk3(df$SMQ020, df$SMQ040)
  df$inpop <- with(df, RIDAGEYR >= 20 & !is.na(A) & !is.na(Y))
  attr(df, "covs") <- c("age","sex","race","educ","pir","cycle","married","hhsize","insured","bmi","smk")
  df
}

derive_E3 <- function(df) {                       # ever-vaped -> hypertension (adults 20+, 2015-18)
  ev <- na_if_in(df$SMQ900, c(7, 9))              # 1 ever e-cig, 2 never
  df$A <- ifelse(ev == 1, 1L, ifelse(ev == 2, 0L, NA_integer_))
  df$Y <- htn_outcome(df)
  df$age <- df$RIDAGEYR; df$sex <- factor(df$RIAGENDR, 1:2, c("Male","Female"))
  df$race <- race4(df$RIDRETH1); df$educ <- educ3(df$DMDEDUC2)
  df$pir <- df$INDFMPIR; df$cycle <- factor(df$SDDSRVYR)
  df$smk <- smk3(df$SMQ020, df$SMQ040)            # combustible smoking = key confounder
  df$cigday <- ifelse(df$smk == "Current", na_if_in(df$SMD650, c(777, 999)), 0)  # 0/day if not a current smoker (avoids structural NA)
  df$bmi <- df$BMXBMI
  df$pa <- factor(ifelse(na_if_in(df$PAQ650, c(7,9)) == 1 | na_if_in(df$PAQ665, c(7,9)) == 1, "Yes","No"))
  df$diab <- factor(ifelse(na_if_in(df$DIQ010, c(7,9)) == 1, "Yes","No"))   # pre-exposure confounder here
  df$inpop <- with(df, RIDAGEYR >= 20 & !is.na(A) & !is.na(Y))
  attr(df, "covs") <- c("age","sex","race","educ","pir","cycle","smk","cigday","bmi","pa","diab")
  df
}

derive_E4 <- function(df) {                       # GDM -> hypertension (ever-pregnant non-diabetic women 20+)
  gdm <- na_if_in(df$RHQ162, c(7, 9))             # 1 yes GDM, 2 no
  df$A <- ifelse(gdm == 1, 1L, ifelse(gdm == 2, 0L, NA_integer_))
  df$Y <- htn_outcome(df)
  everpreg <- na_if_in(df$RHQ131, c(7, 9)); diab <- na_if_in(df$DIQ010, c(7, 9))
  df$age <- df$RIDAGEYR; df$race <- race4(df$RIDRETH1); df$educ <- educ3(df$DMDEDUC2)
  df$pir <- df$INDFMPIR; df$cycle <- factor(df$SDDSRVYR)
  df$npreg <- na_if_in(df$RHQ160, c(77, 99)); df$parity <- na_if_in(df$RHQ171, c(77, 99))
  df$agefirst <- na_if_in(df$RHD180, c(7777, 9999))
  df$menop <- factor(ifelse(na_if_in(df$RHQ031, c(7,9)) == 2, "Postmenopausal","Premenopausal"))  # 2=no longer menstruating
  df$famdm <- factor(ifelse(na_if_in(df$MCQ300C %||% df$MCQ300c, c(7,9)) == 1, "Yes","No"))
  df$smk <- smk3(df$SMQ020, df$SMQ040); df$bmi <- df$BMXBMI   # BMI = pre-pregnancy adiposity PROXY (caveat)
  # eligibility: ever-pregnant, NON-diabetic (exclude prevalent diabetes), female 20+
  df$inpop <- with(df, RIAGENDR == 2 & RIDAGEYR >= 20 & everpreg == 1 & diab != 1 & !is.na(A) & !is.na(Y))
  df$inpop[is.na(df$inpop)] <- FALSE
  attr(df, "covs") <- c("age","race","educ","pir","cycle","npreg","parity","agefirst","menop","famdm","smk","bmi")
  df
}

DERIVERS <- list(E1 = derive_E1, E2 = derive_E2, E3 = derive_E3, E4 = derive_E4)

## ============================ build loop =====================================
summ <- list()
for (id in names(EXAMPLES)) {
  ex <- EXAMPLES[[id]]
  df <- merge_example(ex)
  df <- harmonize_bp(df)
  df <- DERIVERS[[id]](df)
  # restrict to the MEC-examined design set (WTMEC2YR>0); keep ALL such rows (domain via inpop)
  df <- df[!is.na(df$WTMEC2YR) & df$WTMEC2YR > 0, ]
  df$WTMEC_POOLED <- df$WTMEC2YR / length(ex$cycles)
  covs <- attr(df, "covs")
  keep <- c("SEQN","SDMVSTRA","SDMVPSU","WTMEC2YR","WTMEC_POOLED","SDDSRVYR","CYCLE_SUF",
            "A","Y","inpop", covs)
  out <- df[, intersect(keep, names(df))]
  attr(out, "covs") <- covs; attr(out, "example") <- ex$label
  saveRDS(out, file.path(out_dir, paste0(id, "_analytic.rds")))
  d <- out[out$inpop, ]
  summ[[id]] <- data.frame(
    example = id, label = ex$label, cycles = length(ex$cycles),
    n_MEC = nrow(out), n_domain = nrow(d),
    A_prev = round(mean(d$A), 3), Y_prev = round(mean(d$Y), 3),
    cov_anymiss = round(mean(rowSums(is.na(d[, covs, drop = FALSE])) > 0), 3))
}
cat("\n================ ANALYTIC DATASET SUMMARY ================\n")
print(do.call(rbind, summ), row.names = FALSE)
cat("\nSaved per-example analytic frames to", out_dir, "\n")
