# =====================================================================
# nhanes_e1.R  —  R15_aipw_benchmark: survey-weighted AIPW on NHANES E1
# (spec item A7: "...and on one NHANES example (E1, best overlap)";
#  Writing/comments/phase4-arc-sim-specs.md)
#
# SINGLE task (no array). Mirrors the locked NHANES drivers exactly:
#   * data loading + one-hot covariate encoding byte-identical to
#     Nhanes/R/03_run_estimators.R (E1_imputed.rds, attr "covs"/"example",
#     model.matrix dummies, strata=SDMVSTRA, cluster=SDMVPSU,
#     weight=WTMEC_POOLED)
#   * domain handling: nuisances fit on inpop; variances on the FULL design
#     (nest=TRUE, SDMVPSU repeats across strata) -- the .se_des subset-the-
#     DESIGN logic for se_lin, and for se_jkn the replicate design is built
#     on the FULL design then subset() to the domain (fallback: domain-zeroed
#     full-sample ratio; realized mode recorded in the output as jkn_mode).
#   * B-split convention: seeds 20260607 + 1:B (NHANES_B, default 20) -- the
#     SAME split-seed base as Nhanes/R/nhanes_arc.R / R05, so split-MC
#     variability is comparable across runs.
#   * primary library c("SL.glm","SL.earth","SL.glmnet")
#     (Nhanes/R/03_run_estimators.R); override via R15_LIB_OVERRIDE
#     (comma-separated; the SMOKE wiring check uses "SL.glm").
#
# Output: results/arc/R15_aipw_nhanes_E1.csv  (per arm: b, b_split_sd,
# se_jkn, se_lin, df, lcl/ucl [JKn t-df], lcl_lin/ucl_lin, B, library,
# jkn_mode) + a full RDS and a manifest under $R15_E1_OUT.
#
# Env vars (local fallbacks):
#   SLURM_CPUS_PER_TASK  cores                          (default detectCores-1)
#   NHANES_B             cross-fit-split repeats        (default 20; SMOKE -> 2)
#   SIM_CODE             engine dir (codes/)            (default "codes")
#   NH_ANA               analytic dir (imputed rds)     (default REPO_ROOT/Nhanes/analytic)
#   R15_DIR              this run's folder              (default REPO_ROOT/codes/arc_runs/R15_aipw_benchmark)
#   R15_E1_OUT           output dir                     (default REPO_ROOT/Nhanes/nhanes_output/arc_runs/R15_e1)
#   R15_LIB_OVERRIDE     comma-sep SL library override  (default primary)
#   SMOKE                "1" -> wiring check: B=2, SL.glm, 2 cores, SMOKE_ outputs
# =====================================================================

CODE      <- Sys.getenv("SIM_CODE", "R")
REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
HERE      <- Sys.getenv("R15_DIR",
                        file.path(REPO_ROOT, "simulation", "enhancements", "R15_aipw_benchmark"))
source(file.path(CODE, "estimators.R"))               # .sl, make_cf_folds, .se_des (read-only)
source(file.path(HERE, "aipw_helpers.R"))             # aipw_arms(), .aipw_jkn()
source(file.path(CODE, "arc_runs", "_checkpoint.R"))  # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(SuperLearner) })

geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")
cores <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
B     <- geti("NHANES_B", 20L)
ANA   <- Sys.getenv("NH_ANA", file.path(REPO_ROOT, "Nhanes", "analytic"))
OUT   <- Sys.getenv("R15_E1_OUT",
                    file.path(REPO_ROOT, "Nhanes", "nhanes_output", "arc_runs", "R15_e1"))
MAN   <- file.path(OUT, "manifest")
ARC_RESULTS <- file.path(REPO_ROOT, "results", "arc")
for (d in c(OUT, MAN, ARC_RESULTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- library: primary, overridable -------------------------------
LIB <- c("SL.glm", "SL.earth", "SL.glmnet")           # primary (03_run_estimators.R)
ov  <- Sys.getenv("R15_LIB_OVERRIDE", "")
if (SMOKE && !nzchar(ov)) ov <- "SL.glm"              # cheap wiring-check default
if (nzchar(ov)) LIB <- trimws(strsplit(ov, ",")[[1]])
if (SMOKE) { B <- min(B, 2L); cores <- max(1L, min(cores, 2L)) }
G_FLOOR <- 0.05; V_CF <- 5L                           # engine OOF floor + folds

# ---- checkpoint (single task) ----------------------------------------------
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr15_nhanes_E1.rds", tag))
if (!SMOKE) arc_skip_if_done(fn, 1L)

# ---- load + encode E1 (byte-identical to Nhanes/R/03_run_estimators.R) ------
obs  <- readRDS(file.path(ANA, "E1_imputed.rds"))
covs <- attr(obs, "covs"); label <- attr(obs, "example"); sub <- which(obs$inpop)
covdf <- droplevels(obs[sub, covs, drop = FALSE])
mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
stopifnot(nrow(mm) == length(sub))
wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
for (j in seq_along(wnames)) { obs[[wnames[j]]] <- NA_real_; obs[[wnames[j]]][sub] <- mm[, j] }
obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED

cat(sprintf("[R15-E1] %s | n_domain=%d | A_prev=%.3f Y_prev=%.3f | B=%d | cores=%d | lib={%s} | floor=%.2f\n",
            label, length(sub), mean(obs$A[sub]), mean(obs$Y[sub]),
            B, cores, paste(LIB, collapse = ","), G_FLOOR))

# ---- B repeated cross-fit-split runs (mirror nhanes_arc.R) -------------------
one_split <- function(seed) {
  set.seed(seed)
  r <- aipw_arms(obs, learners = LIB, g_floor = G_FLOOR, V_cf = V_CF,
                 W_cols = wnames, nest = TRUE, inpop = obs$inpop)
  list(results = r$results, drow = r$diagnostics$drow)
}
seeds <- 20260607L + seq_len(B)                       # SAME seed base as nhanes_arc.R
ncl <- max(1L, min(cores, B))
if (ncl > 1L) {
  cl <- makeCluster(ncl); on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    HERE <- Sys.getenv("R15_DIR",
                       file.path(Sys.getenv("REPO_ROOT", "."),
                                 "simulation", "enhancements", "R15_aipw_benchmark"))
    source(file.path(CODE, "estimators.R"))
    source(file.path(HERE, "aipw_helpers.R"))
    suppressMessages({ library(survey); library(SuperLearner); library(earth); library(glmnet) })
    TRUE
  }))
  clusterExport(cl, c("obs", "LIB", "wnames", "G_FLOOR", "V_CF"), envir = environment())
  runs <- parLapply(cl, seeds, one_split)
} else {
  runs <- lapply(seeds, one_split)
}

# ---- summarize across the B splits: per-arm mean estimate + split-SD ---------
allres <- do.call(rbind, Map(function(r, b) cbind(split = b, r$results), runs, seq_along(runs)))
arms <- unique(allres$method)
summ <- do.call(rbind, lapply(arms, function(m) {
  s <- allres[allres$method == m, ]
  data.frame(example = "E1", label = label, method = m, B = B,
             library = paste(LIB, collapse = "+"),
             b = mean(s$b), b_split_sd = sd(s$b),
             se_jkn = mean(s$se_jkn), se_lin = mean(s$se_lin), df = s$df[1],
             stringsAsFactors = FALSE)
}))
crit <- qt(0.975, pmax(1, summ$df))
summ$lcl <- summ$b - crit * summ$se_jkn; summ$ucl <- summ$b + crit * summ$se_jkn   # JKn primary
summ$lcl_lin <- summ$b - crit * summ$se_lin; summ$ucl_lin <- summ$b + crit * summ$se_lin
summ$jkn_mode    <- runs[[1]]$drow$jkn_mode      # AIPW-SF arm's realized JKn mode
summ$jkn_mode_cf <- runs[[1]]$drow$jkn_mode_cf   # AIPW-CF arm's (divergence visible)

# ---- save: full RDS + CSV (SMOKE_ prefixed in smoke; locked CSV untouched) ---
drows <- do.call(rbind, Map(function(r, b) cbind(split = b, r$drow), runs, seq_along(runs)))
out <- list(run = "R15_aipw_benchmark_nhanes_E1", example = "E1", label = label,
            learners = LIB, B = B, seeds = seeds, g_floor = G_FLOOR, V_cf = V_CF,
            n_domain = length(sub), summary = summ, per_split = allres,
            per_split_diag = drows)
saveRDS(out, fn)
csv_fn <- file.path(ARC_RESULTS, paste0(tag, "R15_aipw_nhanes_E1.csv"))
csv <- summ; num <- sapply(csv, is.numeric); csv[num] <- lapply(csv[num], round, 6)
write.csv(csv, csv_fn, row.names = FALSE)
cat(sprintf("[R15-E1] saved %s\n[R15-E1] saved %s\n", fn, csv_fn))

# ---- manifest (mirror R05's NHANES manifest; sim POP seeds are N/A here) -----
manifest <- list(
  run = "R15_aipw_benchmark_nhanes_E1", example = "E1", label = label,
  learners = LIB, B = B, g_floor = G_FLOOR, V_cf = V_CF,
  seeds = list(split_seed_base = 20260607L, split_seeds = seeds),
  n_domain = length(sub), covariates = covs, dummy_cols = wnames,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]])
saveRDS(manifest, file.path(MAN, sprintf("%smanifest_nhanes_E1.rds", tag)))

# ---- print: AIPW summary + locked five-arm context (if available) ------------
cat("\n[R15-E1] AIPW arms (mean over", B, "splits):\n")
print(format(summ[, c("method","b","b_split_sd","se_jkn","se_lin","df","lcl","ucl","jkn_mode","jkn_mode_cf")],
             digits = 4), row.names = FALSE)
lk_fn <- file.path(REPO_ROOT, "Nhanes", "results", "E1", "estimators_all_arms.rds")
if (file.exists(lk_fn)) {
  lk <- readRDS(lk_fn)
  cat("\n[R15-E1] locked five-arm E1 results (single split, primary pipeline) for context:\n")
  print(format(lk$results[, c("method","b","se","lcl","ucl","df")], digits = 4),
        row.names = FALSE)
}

# ---- SMOKE gate ---------------------------------------------------------------
if (SMOKE) {
  ok_arms   <- setequal(summ$method, c("AIPW-SF", "AIPW-CF"))
  ok_finite <- all(is.finite(allres$b), is.finite(allres$se_jkn), is.finite(allres$se_lin))
  rat       <- summ$se_jkn / summ$se_lin
  ok_ratio  <- all(is.finite(rat)) && all(rat >= 0.5 & rat <= 2)
  ok_mode   <- all(drows$jkn_mode    %in% c("subset_repdes", "domain_zeroed_full")) &&
               all(drows$jkn_mode_cf %in% c("subset_repdes", "domain_zeroed_full"))
  if (ok_arms && ok_finite && ok_ratio && ok_mode) {
    cat(sprintf("[SMOKE-GATE] PASS: both arms finite; se_jkn/se_lin in [0.5,2] (%.3f, %.3f); jkn_mode=%s jkn_mode_cf=%s.\n",
                rat[1], rat[2], drows$jkn_mode[1], drows$jkn_mode_cf[1]))
  } else {
    cat(sprintf("[SMOKE-GATE] STOP: ok_arms=%s ok_finite=%s ok_ratio=%s ok_mode=%s -- inspect before submitting.\n",
                ok_arms, ok_finite, ok_ratio, ok_mode))
  }
}
