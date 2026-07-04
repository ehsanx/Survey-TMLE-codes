# =====================================================================
# nhanes/arc_runs/R06_mi/run.R  —  R06_mi enhancement run driver
#
# GOAL: replace the single imputation used for the headline NHANES Table 2
# (nhanes/02_impute.R uses mice m = 1) with MULTIPLE IMPUTATION m = 40 for the
# PRIMARY estimator on ALL FOUR examples, and combine the per-imputation
# estimates/SEs by Rubin's rules to produce MI-based intervals. m = 40 follows
# White, Royston & Wood (2011): m >= the percentage of incomplete cases (E4 = 29%).
#
# PRIMARY: the 3-learner near-Donsker library
#   LIB = c('SL.glm','SL.earth','SL.glmnet')
# i.e. the SAME `LIB` as nhanes/03_run_estimators.R (the L2-equivalent rung),
# NOT the full L1-L4 ladder. We deliberately do NOT MI the ladder (that would be
# 10x the expensive deep rungs).
#
# For each example:
#   1. impute m = 40 completed datasets (reuse 02_impute.R's mice setup, m = 40)
#   2. run run_estimators(... LIB, inpop = domain) on EACH completed dataset
#   3. pool point estimates + SEs ACROSS imputations BY ARM via Rubin's rules
#      (between + within variance, t reference with Barnard-Rubin df)
#
# The reported headline arm is Fully-Aware-CF (the primary estimator),
# but we pool ALL arms returned by run_estimators so the full Table-2 column is
# available with MI intervals.
#
# NON-DESTRUCTIVE: SOURCES the canonical engine (R/estimators.R) read-only;
# all new logic is in this run's mi_helpers.R. Outputs go to this run's own
# folders and NEVER clobber the locked single-imputation results.
#
# Env vars (local-test fallbacks):
#   SMOKE                '1' -> tiny run (M=2, examples E2,E3 only)   (default 0)
#   SLURM_ARRAY_TASK_ID  1-based task = example index 1..4             (default: all)
#   SLURM_CPUS_PER_TASK  cores (used to parallelize the M imputations) (default detectCores-1)
#   MI_M                 number of imputations                         (default 40)
#   MI_SEED              base mice seed                                (default 20260607)
#   SIM_CODE             engine dir (R/)                           (default R)
#   NH_ANA               analytic dir (*_analytic.rds)                 (default nhanes/analytic)
#   R06_OUT              per-example output dir   (default nhanes/nhanes_output/arc_runs/R06)
#   R06_RESULTS          summary CSV dir          (default results/arc)
#   R06_MANIFEST         manifest dir             (default nhanes/nhanes_output/arc_runs/R06/manifest)
# =====================================================================

t_start <- proc.time()[3]
RUN_ID <- "R06_mi"

# ---- locate this run's folder (for sourcing mi_helpers.R) --------------------
.this_file <- tryCatch({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else NA_character_
}, error = function(e) NA_character_)
RUN_DIR <- if (!is.na(.this_file)) {
             dirname(.this_file)
           } else file.path(Sys.getenv("REPO_ROOT", "."), "nhanes", "arc_runs", "R06_mi")

# ---- canonical engine (read-only) --------------------------------------------
CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "estimators.R"))      # exposes run_estimators + building blocks
source(file.path(RUN_DIR, "mi_helpers.R"))   # impute_m, rubin_pool, encode_domain
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner)
                   library(earth); library(glmnet); library(mice) })

# ---- config -----------------------------------------------------------------
SMOKE <- identical(Sys.getenv("SMOKE"), "1")
geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
M     <- geti("MI_M",   if (SMOKE) 2L else 40L)   # primary headline run used m = 40
SEED  <- geti("MI_SEED", 20260607L)
cores <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
ANA   <- Sys.getenv("NH_ANA", file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "analytic"))
OUT   <- Sys.getenv("R06_OUT",      file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "nhanes_output", "arc_runs", "R06"))
RES   <- Sys.getenv("R06_RESULTS",  file.path(Sys.getenv("REPO_ROOT","."), "results", "arc"))
MAN   <- Sys.getenv("R06_MANIFEST", file.path(OUT, "manifest"))
for (d in c(OUT, RES, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LIB <- c("SL.glm", "SL.earth", "SL.glmnet")   # primary (== 03_run_estimators.R)

# ---- which examples this invocation handles ----------------------------------
ALL_EX <- c("E1", "E2", "E3", "E4")
task   <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (SMOKE) {
  # decision-rule smoke: the two BORDERLINE examples first (E2 significant, E3 null)
  examples <- c("E2", "E3")
} else if (nzchar(task)) {
  ti <- as.integer(task); stopifnot(ti >= 1L, ti <= length(ALL_EX))
  examples <- ALL_EX[ti]
} else {
  examples <- ALL_EX                          # local: do all four serially
}

cat(sprintf("[%s] SMOKE=%s  M=%d  cores=%d  examples={%s}  LIB={%s}\n",
            RUN_ID, SMOKE, M, cores, paste(examples, collapse=","), paste(LIB, collapse=",")))

# =====================================================================
# per-example MI pipeline
# =====================================================================
run_example_mi <- function(ex_id) {
  # per-example checkpoint: reuse a finished example and CONTINUE the loop (do NOT
  # quit -- in serial/local mode all four examples are looped and the combined CSV
  # is written below; quitting on the first cached hit would skip that combine).
  fn <- file.path(OUT, sprintf("mi_%s.rds", ex_id))
  if (!SMOKE && file.exists(fn) && isTRUE(file.info(fn)$size > 200L)) {
    cached <- tryCatch(readRDS(fn), error = function(e) NULL)
    if (!is.null(cached)) {
      cat(sprintf("[%s] SKIP checkpoint: %s already complete -> reusing\n", ex_id, basename(fn)))
      return(cached)
    }
  }

  obs0 <- readRDS(file.path(ANA, paste0(ex_id, "_analytic.rds")))
  label <- attr(obs0, "example")
  cat(sprintf("\n=== %s (%s): imputing m=%d ...\n", ex_id, label, M))

  # (1) m completed datasets (mice m=M; same setup as 02_impute.R) ------------
  comp <- impute_m(obs0, M = M, seed = SEED)
  imp_covs <- attr(comp, "imputed_covs")

  # (2) run the primary on each completed dataset --------------------
  # Each imputation gets its OWN deterministic seed so the CF/SL splits are
  # reproducible; different across imputations (split variability is averaged
  # over, like the rest of the MI noise).
  one_imp <- function(k) {
    enc <- encode_domain(comp[[k]])
    set.seed(SEED + 1000L * k)
    r <- run_estimators(enc$obs, learners = LIB, V_cf = 5L, inner_cv_folds = 5L,
                        W_cols = enc$wnames, nest = TRUE, inpop = enc$obs$inpop)
    r$results$imp <- k
    r$results
  }
  ncl <- max(1L, min(cores, M))
  if (ncl > 1L) {
    cl <- makeCluster(ncl)
    on.exit(stopCluster(cl), add = TRUE)
    invisible(clusterEvalQ(cl, {
      CODE <- Sys.getenv("SIM_CODE", "R")
      source(file.path(CODE, "estimators.R"))
      suppressMessages({ library(survey); library(tmle); library(SuperLearner)
                         library(earth); library(glmnet) })
      TRUE
    }))
    clusterExport(cl, c("comp", "LIB", "SEED", "encode_domain"), envir = environment())
    per_imp <- parLapply(cl, seq_len(M), one_imp)
    stopCluster(cl); on.exit()
  } else {
    per_imp <- lapply(seq_len(M), one_imp)
  }
  per_imp_df <- do.call(rbind, per_imp)        # method,b,se,df,imp  (M rows per arm)

  # (3) Rubin pool per arm -----------------------------------------------------
  arms <- unique(per_imp_df$method)
  pooled <- do.call(rbind, lapply(arms, function(m) {
    s  <- per_imp_df[per_imp_df$method == m, ]
    rp <- rubin_pool(s$b, s$se, df_complete = s$df[1])  # all M share the design df
    data.frame(example = ex_id, label = label, method = m,
               b = rp$b, se = rp$se, df = rp$df,
               b_imp_sd = sd(s$b),                       # raw spread of point ests
               within_se = sqrt(rp$Ubar), between_sd = sqrt(rp$B),
               riv = rp$riv, fmi = rp$fmi, M = M,
               stringsAsFactors = FALSE)
  }))
  crit <- qt(0.975, pmax(1, pooled$df))
  pooled$lcl <- pooled$b - crit * pooled$se
  pooled$ucl <- pooled$b + crit * pooled$se
  pooled$signif <- (pooled$lcl > 0) | (pooled$ucl < 0)   # CI excludes 0

  # save per-example object (per-imputation rows + pooled + meta) --------------
  out <- list(example = ex_id, label = label, M = M, lib = LIB,
              imputed_covs = imp_covs, per_imp = per_imp_df, pooled = pooled,
              A_prev = mean(obs0$A[obs0$inpop]), Y_prev = mean(obs0$Y[obs0$inpop]),
              n_domain = sum(obs0$inpop))
  saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
  cat(sprintf("[%s] saved %s\n", ex_id, fn))

  # console echo: the primary + a quick all-arm view ----------------
  prim <- pooled[pooled$method == "Fully-Aware-CF", ]
  cat(sprintf("  PRIMARY (Fully-Aware-CF): b=%.4f  se=%.4f  CI=[%.4f, %.4f]  df=%.1f  fmi=%.3f  signif=%s\n",
              prim$b, prim$se, prim$lcl, prim$ucl, prim$df, prim$fmi, prim$signif))
  print(pooled[, c("method","b","se","lcl","ucl","df","fmi","signif")], row.names = FALSE)
  out
}

results <- lapply(examples, run_example_mi)

# =====================================================================
# combined summary CSV (matches the locked Table-2 column conventions)
# =====================================================================
comb <- do.call(rbind, lapply(results, `[[`, "pooled"))
# round numeric cols for the CSV (same convention as aggregate_nhanes.R)
num <- sapply(comb, is.numeric); comb_print <- comb; comb_print[num] <- round(comb_print[num], 4)
suffix  <- if (SMOKE) "_SMOKE" else ""
csv_all <- file.path(RES, sprintf("R06_mi_summary%s.csv", suffix))
write.csv(comb_print, csv_all, row.names = FALSE)

# primary-only headline CSV (the actual Table-2 MI column) ----------
prim_tbl <- comb_print[comb_print$method == "Fully-Aware-CF", ]
csv_prim <- file.path(RES, sprintf("R06_mi_primary%s.csv", suffix))
write.csv(prim_tbl, csv_prim, row.names = FALSE)

cat(sprintf("\nDONE. %d example(s) in %.1f min.\n -> %s\n -> %s\n",
            length(examples), (proc.time()[3] - t_start)/60, csv_all, csv_prim))
cat("\n==== PRIMARY (Fully-Aware-CF) MI Table-2 column ====\n")
print(prim_tbl[, c("example","label","b","se","lcl","ucl","df","fmi","signif")], row.names = FALSE)

# =====================================================================
# reproducibility manifest (mirrors nhanes_arc.R / run_sim.R)
# =====================================================================
manifest <- list(
  run_id = RUN_ID, examples = examples, M = M, seed = SEED, lib = LIB,
  smoke = SMOKE,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","glmnet","mice"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  elapsed_min = round((proc.time()[3] - t_start)/60, 2))
saveRDS(manifest, file.path(MAN, sprintf("manifest_%s%s.rds",
        paste(examples, collapse = "-"), suffix)))
cat("\nManifest saved to", MAN, "\n")
