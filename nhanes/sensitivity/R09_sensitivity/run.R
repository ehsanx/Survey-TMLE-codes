# =====================================================================
# nhanes/arc_runs/R09_sensitivity/run.R  —  R09_sensitivity driver
#
# GOAL (pre-empt an identification-grounds major revision): provide LOCAL re-runs
# (laptop-runnable, no SLURM needed) of the PRIMARY estimator on
# MODIFIED adjustment sets, and report the ATE/CI delta vs the main analysis:
#
#   (a) E4  with BMI OMITTED from the adjustment set
#       (the pre-specified mediator-vs-confounder sensitivity: BMI is a
#        pre-pregnancy adiposity PROXY -> if BMI sits on the GDM -> HTN path it is
#        a mediator and should be dropped; we show the estimate either way).
#   (b) E2  treating BMI as a MEDIATOR (drop bmi) vs confounder
#       (food insecurity -> depression; BMI may mediate via adiposity).
#   (c) E1  dropping PHQ-9 (alternate adjustment)
#       (short sleep -> obesity; depression as confounder vs questionable adj.).
#   (d) an empirical SDMVSTRA cross-cycle NON-OVERLAP CHECK across all four
#       examples (verify strata do not span cycles in the pooled design -> the
#       pooled stratification / nest=TRUE is valid). A quick tabulation, NOT an
#       estimator run.
#
# Each sensitivity = THE PRIMARY run on the modified adjustment set.
# PRIMARY library == nhanes/03_run_estimators.R:
#   LIB = c('SL.glm','SL.earth','SL.glmnet')      (the L2-equivalent near-Donsker)
# Headline arm = Fully-Aware-CF; we keep all five arms for completeness.
#
# NON-DESTRUCTIVE: SOURCES the canonical engine (R/estimators.R) read-only and
# the already-IMPUTED analytic frames (nhanes/analytic/<id>_imputed.rds). All new
# logic is in this run's sensitivity_helpers.R. Outputs go to this run's OWN
# folders and NEVER clobber the locked headline results.
#
# This is a covariate-set swap loop + one design check: LIGHT. Designed to finish
# on a laptop in a few minutes; an optional submit.slurm is included for
# convenience only.
#
# Env vars (all have local-test fallbacks):
#   SMOKE         '1' -> tiny run: only the (a) E4-no-BMI cell + the design check,
#                        cheap GLM-only library, finishes in 1-3 min   (default 0)
#   SIM_CODE      engine dir (R/)                       (default R)
#   NH_ANA        analytic dir (*_imputed.rds)              (default nhanes/analytic)
#   R09_OUT       per-cell output dir   (default nhanes/nhanes_output/arc_runs/R09)
#   R09_RESULTS   summary CSV dir       (default results/arc)
#   R09_MANIFEST  manifest dir          (default <R09_OUT>/manifest)
#   R09_SEED      run (CF/SL split) seed                    (default 20260607)
# =====================================================================

t_start <- proc.time()[3]
RUN_ID <- "R09_sensitivity"

# ---- locate this run's folder (to source sensitivity_helpers.R) --------------
.this_file <- tryCatch({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else NA_character_
}, error = function(e) NA_character_)
RUN_DIR <- if (!is.na(.this_file)) {
             dirname(.this_file)
           } else {
             file.path(Sys.getenv("REPO_ROOT", "."), "nhanes", "arc_runs", "R09_sensitivity")
           }

# ---- canonical engine (read-only) + this run's helpers -----------------------
CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "estimators.R"))               # run_estimators + building blocks
source(file.path(RUN_DIR, "sensitivity_helpers.R"))   # set_covs/encode_domain/run_primary/ci_delta/strata_cycle_check
source(file.path(CODE, "arc_runs", "_checkpoint.R"))  # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(survey); library(tmle); library(SuperLearner)
                   library(earth); library(glmnet) })

# ---- config ------------------------------------------------------------------
SMOKE <- identical(Sys.getenv("SMOKE"), "1")
geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SEED  <- geti("R09_SEED", 20260607L)
ANA   <- Sys.getenv("NH_ANA",      file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "analytic"))
OUT   <- Sys.getenv("R09_OUT",     file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "nhanes_output", "arc_runs", "R09"))
RES   <- Sys.getenv("R09_RESULTS", file.path(Sys.getenv("REPO_ROOT","."), "results", "arc"))
MAN   <- Sys.getenv("R09_MANIFEST", file.path(OUT, "manifest"))
for (d in c(OUT, RES, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# primary (== 03_run_estimators.R); SMOKE uses GLM-only for speed
LIB <- if (SMOKE) "SL.glm" else c("SL.glm", "SL.earth", "SL.glmnet")

# ---- the sensitivity grid ----------------------------------------------------
# Each entry: example, a short cell id, the covariate(s) to DROP, and a human
# label for the adjustment-set change. The MAIN analysis for each example is the
# full adjustment set (drop = none) -- recomputed here under the SAME LIB/seed so
# the delta is apples-to-apples (NOT differenced against the locked Table-2 CSV,
# whose seed/library could differ).
SENS_GRID <- list(
  list(example = "E4", cell = "E4_noBMI", drop = "bmi",
       desc = "(a) E4 BMI omitted: mediator-vs-confounder (pre-pregnancy adiposity proxy)"),
  list(example = "E2", cell = "E2_noBMI", drop = "bmi",
       desc = "(b) E2 BMI as mediator (dropped) vs confounder"),
  list(example = "E1", cell = "E1_noPHQ", drop = "phq",
       desc = "(c) E1 PHQ-9 dropped: alternate adjustment"),
  # E3 over-adjustment sensitivity (added per the e3-sensitivity adversarial review,
  # 2026-06-22): E3 is the one example whose conclusion flips, so its identification
  # deserves a robustness check. Drop the two adjustments whose confounder-vs-mediator
  # role is most contestable -- diabetes (a pre-exposure confounder in E3 but a mediator
  # in E4, the paper's own caveat) and BMI (a potential mediator e-cig -> adiposity ->
  # BP, and the influential adjuster) -- separately and combined, so a moved estimate is
  # attributable. Expected: the design-honest null survives every removal (covers 0).
  list(example = "E3", cell = "E3_noDIAB", drop = "diab",
       desc = "(e) E3 diabetes dropped: contested confounder-vs-mediator (confounder in E3, mediator in E4)"),
  list(example = "E3", cell = "E3_noBMI", drop = "bmi",
       desc = "(f) E3 BMI dropped: potential mediator (e-cig -> adiposity -> hypertension); the influential adjuster"),
  list(example = "E3", cell = "E3_noBOTH", drop = c("bmi", "diab"),
       desc = "(g) E3 BMI and diabetes both dropped: combined over-adjustment stress")
)
if (SMOKE) SENS_GRID <- SENS_GRID[1]   # smoke: just (a) E4-no-BMI

cat(sprintf("[%s] SMOKE=%s  SEED=%d  LIB={%s}  cells={%s}\n",
            RUN_ID, SMOKE, SEED, paste(LIB, collapse = ","),
            paste(sapply(SENS_GRID, `[[`, "cell"), collapse = ",")))

# =====================================================================
# (a)-(c): per-cell sensitivity = primary on the modified covs,
#          differenced against the SAME example's full-adjustment fit.
# =====================================================================
main_cache <- list()   # cache the full-adjustment fit per example (reuse across cells)

run_cell <- function(g) {
  ex <- g$example
  # ---- per-cell checkpoint: resume mid-loop if this cell already finished ----
  # NOT a SLURM array -- this is ONE job over many cells. So instead of quitting
  # the whole job (arc_skip_if_done), test the cell's own output rds and, if a
  # valid one exists, reuse it and skip the heavy run_primary fits for this cell.
  fn <- file.path(OUT, sprintf("sens_%s.rds", g$cell))
  if (!SMOKE && file.exists(fn) && isTRUE(file.info(fn)$size > 200L)) {
    cached <- tryCatch(readRDS(fn), error = function(e) NULL)
    if (!is.null(cached)) {
      cat(sprintf("[%s] SKIP checkpoint: %s already complete -> reusing\n",
                  RUN_ID, basename(fn)))
      return(cached)
    }
    cat(sprintf("[%s] checkpoint %s corrupt/truncated -> recomputing\n",
                RUN_ID, basename(fn)))
  }
  obs <- readRDS(file.path(ANA, paste0(ex, "_imputed.rds")))
  label <- attr(obs, "example")
  full_covs <- attr(obs, "covs")
  cat(sprintf("\n=== %s | %s ===\n     example: %s\n     full covs:    {%s}\n     dropping:     {%s}\n",
              g$cell, g$desc, label, paste(full_covs, collapse = ","), paste(g$drop, collapse = ",")))

  # MAIN: full adjustment set (cached per example) ----------------------------
  if (is.null(main_cache[[ex]])) {
    main <- run_primary(obs, LIB = LIB, seed = SEED)        # full covs
    main_cache[[ex]] <<- main
  } else main <- main_cache[[ex]]

  # SENSITIVITY: modified adjustment set --------------------------------------
  obs_s <- set_covs(obs, drop = g$drop)
  sens  <- run_primary(obs_s, LIB = LIB, seed = SEED)

  # tag and stack the two all-arm tables --------------------------------------
  main$adjust <- "main";        main$example <- ex; main$label <- label; main$cell <- g$cell
  sens$adjust <- "sensitivity"; sens$example <- ex; sens$label <- label; sens$cell <- g$cell
  arms_tbl <- rbind(main, sens)

  # ATE/CI delta on the primary (and on Fully-Aware single-fit too) --
  delta_cf <- ci_delta(main, sens, arm = "Fully-Aware-CF")
  delta_fa <- ci_delta(main, sens, arm = "Fully-Aware")
  delta <- rbind(delta_cf, delta_fa)
  delta$example <- ex; delta$label <- label; delta$cell <- g$cell; delta$dropped <- paste(g$drop, collapse = ",")

  # save the per-cell object --------------------------------------------------
  out <- list(cell = g$cell, example = ex, label = label, desc = g$desc,
              full_covs = full_covs, dropped = g$drop, lib = LIB, seed = SEED,
              arms = arms_tbl, delta = delta,
              g_fa_main = attr(main, "g_fa"), g_fa_sens = attr(sens, "g_fa"),
              n_domain = attr(sens, "n_domain"))
  saveRDS(out, fn)   # `fn` computed early for the per-cell checkpoint guard

  # console echo --------------------------------------------------------------
  cat("  -- all arms (main vs sensitivity) --\n")
  print(arms_tbl[, c("adjust","method","b","se","lcl","ucl","df","signif")], row.names = FALSE)
  cat("  -- ATE/CI delta (sensitivity - main) --\n")
  print(delta[, c("arm","b_main","b_sens","d_b","lcl_sens","ucl_sens",
                  "signif_main","signif_sens","conclusion_flip")], row.names = FALSE)
  out
}

cells <- lapply(SENS_GRID, run_cell)

# =====================================================================
# (d): empirical SDMVSTRA cross-cycle non-overlap check (all four examples).
#      Quick tabulation; runs even in SMOKE (it is the design-validity check
#      identification reviewers will ask for).
# =====================================================================
cat("\n=== (d) SDMVSTRA cross-cycle NON-OVERLAP check (full MEC design) ===\n")
CHK_EX <- c("E1", "E2", "E3", "E4")
design_check <- do.call(rbind, lapply(CHK_EX, function(ex) {
  obs <- readRDS(file.path(ANA, paste0(ex, "_imputed.rds")))
  c_full <- strata_cycle_check(obs, domain_only = FALSE)
  c_dom  <- strata_cycle_check(obs, domain_only = TRUE)
  data.frame(
    example = ex, label = attr(obs, "example"),
    n_strata = c_full$n_strata, n_cycles = c_full$n_cycles,
    max_cycles_spanned = c_full$max_cycles_spanned,
    n_strata_multicycle = c_full$n_strata_multicycle,
    overlap_clean = c_full$overlap_clean,
    strata_min = c_full$strata_range[1], strata_max = c_full$strata_range[2],
    dom_max_cycles_spanned = c_dom$max_cycles_spanned,
    dom_overlap_clean = c_dom$overlap_clean,
    stringsAsFactors = FALSE)
}))
print(design_check, row.names = FALSE)
all_clean <- all(design_check$overlap_clean)
cat(sprintf("\n  SDMVSTRA cross-cycle non-overlap: %s (every stratum lies in exactly ONE cycle: %s)\n",
            if (all_clean) "CLEAN" else "VIOLATED",
            if (all_clean) "pooled stratification / nest=TRUE valid"
            else "WARNING: some stratum spans >1 cycle -- re-examine pooled design"))

# =====================================================================
# combined summary CSVs (match the locked Table-2 column conventions where
# applicable: example,label,method,b,se,df,lcl,ucl plus the sensitivity tags)
# =====================================================================
suffix <- if (SMOKE) "_SMOKE" else ""

# (1) all arms, main vs sensitivity, every cell
arms_all <- do.call(rbind, lapply(cells, function(o)
  o$arms[, c("example","label","cell","adjust","method","b","se","df","lcl","ucl","signif")]))
num <- sapply(arms_all, is.numeric); arms_all[num] <- round(arms_all[num], 4)
csv_arms <- file.path(RES, sprintf("R09_sensitivity_arms%s.csv", suffix))
write.csv(arms_all, csv_arms, row.names = FALSE)

# (2) the headline delta table (the actual "report the ATE/CI delta" deliverable)
delta_all <- do.call(rbind, lapply(cells, `[[`, "delta"))
dnum <- sapply(delta_all, is.numeric); delta_all[dnum] <- round(delta_all[dnum], 4)
csv_delta <- file.path(RES, sprintf("R09_sensitivity_delta%s.csv", suffix))
write.csv(delta_all[, c("example","label","cell","dropped","arm",
                        "b_main","lcl_main","ucl_main","signif_main",
                        "b_sens","lcl_sens","ucl_sens","signif_sens",
                        "d_b","d_lcl","d_ucl","rel_d_b","conclusion_flip")],
          csv_delta, row.names = FALSE)

# (3) the design-validity check
csv_chk <- file.path(RES, sprintf("R09_design_check%s.csv", suffix))
write.csv(design_check, csv_chk, row.names = FALSE)

cat(sprintf("\nDONE. %d sensitivity cell(s) + design check in %.1f min.\n -> %s\n -> %s\n -> %s\n",
            length(cells), (proc.time()[3] - t_start) / 60, csv_arms, csv_delta, csv_chk))
cat("\n==== ATE/CI DELTA (primary, Fully-Aware-CF) ====\n")
print(delta_all[delta_all$arm == "Fully-Aware-CF",
                c("cell","dropped","b_main","b_sens","d_b","signif_main","signif_sens","conclusion_flip")],
      row.names = FALSE)

# =====================================================================
# reproducibility manifest (mirrors nhanes_arc.R / run_sim.R / R06_mi)
# =====================================================================
manifest <- list(
  run_id = RUN_ID, smoke = SMOKE, seed = SEED, lib = LIB,
  cells = lapply(SENS_GRID, function(g) list(cell = g$cell, example = g$example,
                                             dropped = g$drop, desc = g$desc)),
  design_check = design_check, design_check_clean = all_clean,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","glmnet","ranger","mice"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  elapsed_min = round((proc.time()[3] - t_start) / 60, 2))
saveRDS(manifest, file.path(MAN, sprintf("manifest%s.rds", suffix)))
cat("\nManifest saved to", MAN, "\n")
