# =====================================================================
# nhanes/sensitivity/bench/run_bench_s36.R
#   Web Table S36 -- compute-landscape benchmark: per-estimator / per-estimand
#   wall + CPU time, SINGLE imputation. Pairs with run_bench.R (S35) and is
#   meant to run on the SAME node so the two tables are hardware-comparable.
#
#   Rows: the 5 design-awareness arms (per-arm timing from the engine's
#   diag$timing), survey-weighted AIPW (aipw_helpers), risk-ratio/odds-ratio
#   (delta method on the primary CF arm), the continuous-outcome estimator
#   (Gaussian family, E1 BMI), and a deployable Super Learner (ranger-augmented,
#   ARC-safe). Primary arms + AIPW + RR/OR + deployable on E3 (smallest n);
#   continuous on E1.
#
#   Env: SIM_CODE, NH_ANA, NH_RAW, BENCH_SEED, BENCH_LIB, BENCH_DEPLOY_LIB,
#        BENCH_RESULTS.
# =====================================================================
t0 <- proc.time()[3]
CODE <- Sys.getenv("SIM_CODE", "R")
.this <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(f) else NA }, error = function(e) NA)
RUN_DIR <- if (!is.na(.this)) dirname(.this) else file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "sensitivity", "bench")
source(file.path(CODE, "estimators.R"))                                 # run_estimators (instrumented: diag$timing)
source(file.path(dirname(RUN_DIR), "R06_mi", "mi_helpers.R"))           # impute_m, encode_domain
source(file.path(CODE, "estimands_rr_or.R"))                            # delta_rr_or
AIPW_H <- file.path("simulation", "enhancements", "R15_aipw_benchmark", "aipw_helpers.R")
have_aipw <- file.exists(AIPW_H); if (have_aipw) source(AIPW_H)         # aipw_arms
suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(earth); library(glmnet); library(mice)
                   for (p in c("ranger","xgboost")) suppressWarnings(try(library(p, character.only = TRUE), silent = TRUE)) })
R21_H <- file.path("simulation", "enhancements", "R21_deployable_cvcf", "r21_helpers.R")   # optional: SL.xgboost.tuned / SL.hal9001 wrappers
if (file.exists(R21_H)) suppressMessages(tryCatch(source(R21_H), error = function(e) message("  ! r21_helpers: ", conditionMessage(e))))

SEED <- { x <- Sys.getenv("BENCH_SEED"); if (nzchar(x)) as.integer(x) else 20260607L }
ANA  <- Sys.getenv("NH_ANA", file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "analytic"))
RAW  <- Sys.getenv("NH_RAW", file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "raw"))
RES  <- Sys.getenv("BENCH_RESULTS", file.path(Sys.getenv("REPO_ROOT","."), "results", "arc"))
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE, showWarnings = FALSE)
LIB    <- { x <- Sys.getenv("BENCH_LIB"); if (nzchar(x)) trimws(strsplit(x, ",")[[1]]) else c("SL.glm","SL.earth","SL.glmnet") }
DEPLOY <- { x <- Sys.getenv("BENCH_DEPLOY_LIB"); if (nzchar(x)) trimws(strsplit(x, ",")[[1]]) else c("SL.glm","SL.earth","SL.glmnet","SL.ranger","SL.xgboost") }  # full deployable; ARC has ranger/xgboost/hal9001. Smoke: set BENCH_DEPLOY_LIB=SL.glm,SL.mean

timed <- function(expr) { st <- system.time(v <- tryCatch(expr, error = function(e) { message("  ! ", conditionMessage(e)); NULL }))
                          list(v = v, wall = as.numeric(st[["elapsed"]]), cpu = as.numeric(st[["user.self"]] + st[["sys.self"]])) }
NAr  <- function() list(v = NULL, wall = NA_real_, cpu = NA_real_)

# ---- E3 single imputation: 5 arms with per-arm timing -----------------------
cat("[s36] E3 single imputation: 5 arms (instrumented per-arm timing)\n")
ob0  <- readRDS(file.path(ANA, "E3_analytic.rds"))
comp <- impute_m(ob0, M = 1L, seed = SEED)[[1]]; enc <- encode_domain(comp); ob <- enc$obs; Wc <- enc$wnames
set.seed(SEED)
A   <- timed(run_estimators(ob, learners = LIB, V_cf = 5L, inner_cv_folds = 5L, W_cols = Wc, nest = TRUE, inpop = ob$inpop))
tim <- if (!is.null(A$v)) A$v$diagnostics$timing else numeric(0)
g   <- function(k) if (k %in% names(tim)) as.numeric(tim[k]) else NA_real_

# representative design-based SE cost (one svymean on the CF EIF, nest=TRUE)
t_se <- NA_real_
if (!is.null(A$v)) {
  arm <- A$v$diagnostics$arms[["Fully-Aware-CF"]]; cf_eif <- arm$eif1 - arm$eif0
  fe  <- numeric(nrow(ob)); fe[ob$inpop] <- cf_eif
  dd  <- data.frame(strata = ob$strata, cluster = ob$cluster, weight = ob$weight, e = fe, dom = ob$inpop)
  t_se <- timed({ des <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dd, nest = TRUE); svymean(~e, subset(des, dom)) })$wall
}

# ---- survey-weighted AIPW (independent competitor; SF + CF in one call) ------
cat("[s36] survey-weighted AIPW\n")
AW <- if (have_aipw) timed(aipw_arms(ob, learners = LIB, V_cf = 5L, W_cols = Wc, nest = TRUE, inpop = ob$inpop)) else NAr()
aipw_tim <- if (!is.null(AW$v) && !is.null(AW$v$diagnostics$timing)) AW$v$diagnostics$timing else c(sf = NA_real_, cf = NA_real_)

# ---- risk ratio / odds ratio (delta method on the primary CF arm) ---------
cat("[s36] RR/OR delta method\n")
RO <- NAr()
if (!is.null(A$v)) {
  cfa <- A$v$diagnostics$arms[["Fully-Aware-CF"]]
  RO  <- timed(delta_rr_or(cfa$psi1, cfa$psi0, cfa$eif1, cfa$eif0, ob$strata, ob$cluster, ob$weight,
                           clustered = TRUE, nest = TRUE, inpop = ob$inpop))
}

# ---- continuous outcome (Gaussian; E1 BMI) ----------------------------------
cat("[s36] continuous E1-BMI (Gaussian)\n")
CT <- timed({
  obc  <- readRDS(file.path(ANA, "E1_imputed.rds")); covs <- attr(obc, "covs")
  bmx  <- do.call(rbind, lapply(list.files(RAW, "^BMX_.*\\.rds$", full.names = TRUE),
                                function(f) { d <- readRDS(f); d[, c("SEQN", "BMXBMI")] }))
  obc$BMI <- bmx$BMXBMI[match(obc$SEQN, bmx$SEQN)]; obc$Y <- obc$BMI
  obc$inpop <- obc$inpop & !is.na(obc$BMI); sub <- which(obc$inpop)
  cd <- droplevels(obc[sub, covs, drop = FALSE]); mm <- model.matrix(~ ., data = cd)[, -1, drop = FALSE]
  wn <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wn
  for (j in seq_along(wn)) { obc[[wn[j]]] <- NA_real_; obc[[wn[j]]][sub] <- mm[, j] }
  obc$strata <- obc$SDMVSTRA; obc$cluster <- obc$SDMVPSU; obc$weight <- obc$WTMEC_POOLED
  set.seed(SEED); run_estimators(obc, learners = LIB, V_cf = 5L, inner_cv_folds = 5L, family = "gaussian",
                                 W_cols = wn, nest = TRUE, inpop = obc$inpop)
})
timC <- if (!is.null(CT$v)) CT$v$diagnostics$timing else numeric(0)
gC   <- function(k) if (k %in% names(timC)) as.numeric(timC[k]) else NA_real_

# ---- deployable Super Learner (ranger-augmented; ARC-safe), E3 --------------
cat("[s36] deployable SL {", paste(DEPLOY, collapse = ","), "}\n")
DP   <- timed({ set.seed(SEED); run_estimators(ob, learners = DEPLOY, V_cf = 5L, inner_cv_folds = 5L, W_cols = Wc, nest = TRUE, inpop = ob$inpop) })
timD <- if (!is.null(DP$v)) DP$v$diagnostics$timing else numeric(0)
gD   <- function(k) if (k %in% names(timD)) as.numeric(timD[k]) else NA_real_

# ---- assemble S36 -----------------------------------------------------------
mkrow <- function(estimator, V, refits, pw, dse) {
  pw <- as.numeric(pw); dse <- as.numeric(dse)
  data.frame(estimator = estimator, V = V, refits = refits,
             point_wall = round(pw, 2), point_cpu = round(pw, 2), cores = 1L,
             design_se = round(dse, 3), total = round(ifelse(is.na(pw), NA, pw + ifelse(is.na(dse), 0, dse)), 2),
             stringsAsFactors = FALSE)
}
S36 <- rbind(
  mkrow("Non-Aware (unweighted single-fit)",            1, 1, g("single_fit_unweighted"), 0.00),
  mkrow("Partially-Aware (weighted point)",             1, 1, g("single_fit_weighted"),   t_se),
  mkrow("Fully-Aware (weighted single-fit)",            1, 1, g("single_fit_weighted"),   t_se),
  mkrow("Fully-Aware-CV (internal cross-validation)",   5, 5, g("cv"),                    t_se),
  mkrow("Fully-Aware-CF (primary, cross-fit)",        5, 5, g("cf"),                    t_se),
  mkrow("AIPW, single-fit (survey-weighted)",           1, 1, if ('sf' %in% names(aipw_tim)) as.numeric(aipw_tim['sf']) else NA, t_se),
  mkrow("AIPW, cross-fit (survey-weighted)",            5, 5, if ('cf' %in% names(aipw_tim)) as.numeric(aipw_tim['cf']) else NA, t_se),
  mkrow("Risk ratio / odds ratio (delta)",              5, 0, RO$wall,                    NA),
  mkrow("Continuous outcome (Gaussian, E1)",            5, 5, gC("cf"),                   t_se),
  mkrow("Deployable Super Learner (CF)",                5, 5, gD("cf"),                   t_se)
)
write.csv(S36, file.path(RES, "bench_s36_compute.csv"), row.names = FALSE)
saveRDS(list(S36 = S36, tim_e3 = tim, tim_cont = timC, tim_deploy = timD, t_se = t_se,
             aipw = AW$wall, rror = RO$wall, deploy_lib = DEPLOY,
             arm_results = if (!is.null(A$v)) A$v$results else NULL),
        file.path(RES, "bench_s36.rds"))

cat("\n==== S36 COMPUTE LANDSCAPE (single imputation; E3, continuous on E1) ====\n")
print(S36, row.names = FALSE)
cat(sprintf("design-SE (one svymean) = %.3fs | 5-arm suite wall = %.1fs | AIPW = %.1fs | deployable = %.1fs\n",
            t_se, A$wall, AW$wall, DP$wall))
cat(sprintf("[s36] done in %.1f min.\n", (proc.time()[3] - t0) / 60))
