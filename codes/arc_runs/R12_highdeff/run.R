# =====================================================================
# R12_highdeff / run.R  —  High-DEFF "Design C" stress run (SLURM array)
#
# GOAL: realise an EIF clustering design effect (deff_clust on the
# Fully-Aware EIF) of ~2.5-4 -- far above the headline ~1.25-1.40 -- and
# show that Fully-Aware / Fully-Aware-CF stay ~nominal while
# Partially-Aware (no clustering in the SE) degrades MORE than at low
# DEFF, sharpening the "clustering matters" point.
#
# NON-DESTRUCTIVE: this driver SOURCES the canonical engine
# (config/dgp/estimators/diagnostics/learners) read-only and gets all
# Design C behaviour from design_c.R (this folder), which only PASSES
# OVERRIDDEN ARGUMENTS to the canonical make_population()/draw_sample().
# Nothing in codes/*.R is modified. Outputs go to a private folder and a
# private summary CSV; the locked results are never touched.
#
# Grid: ONE scenario (Design C) x {L1_param, L3_adaptive, L4_aggressive}
#       x rep-chunks. (L2 is dropped: the spec asks PA vs FA vs CF at the
#       parametric, adaptive, and aggressive rungs.)
#
# Env vars (with local fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index                 (default 1)
#   SLURM_CPUS_PER_TASK  cores for this task                (default detectCores-1)
#   SIM_N_REPS           total reps per rung                (default 1000)
#   SIM_CHUNK            reps per array task                (default 100)
#   SIM_CODE             canonical code dir                 (default codes)
#   SIM_OUT              per-task RDS dir   (default sim_output/arc_runs/R12_highdeff)
#   SMOKE                "1" -> tiny single-cell run (100 reps, L4, 1 chunk)
#   DC_SIGMA2_Y/_A/_C, DC_BASE_N0, DC_BASE_M, DC_ALPHA_STRAT  Design C knobs
#     (read in design_c.R; export to RETUNE without editing code)
# =====================================================================

RUN_ID <- "R12_highdeff"

# ---- locate + source the canonical engine (read-only) ----------------
CODE <- Sys.getenv("SIM_CODE", "codes")
if (!file.exists(file.path(CODE, "config.R"))) {
  cand <- c("codes", "codes")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))      # POP_SEED, SAMPLE_SEED_BASE, paths, RNGkind
source(file.path(CODE, "dgp.R"))         # make_population, draw_sample, audit_population
source(file.path(CODE, "estimators.R"))  # run_estimators + building blocks
source(file.path(CODE, "diagnostics.R")) # deff_clust
source(file.path(CODE, "learners.R"))    # SL_LADDER, custom wrappers
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# ---- Design C wrapper (this folder; reuses canonical dgp.R via args) --
# resolve THIS script's directory so design_c.R is found regardless of cwd
.this_dir <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else file.path(REPO_ROOT, "codes/arc_runs", RUN_ID)
}, error = function(e) file.path(REPO_ROOT, "codes/arc_runs", RUN_ID))
source(file.path(.this_dir, "design_c.R"))

# ---- env helpers -----------------------------------------------------
geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- Sys.getenv("SMOKE") == "1"

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
OUT    <- Sys.getenv("SIM_OUT",
                     file.path(DATA_ROOT, "arc_runs", RUN_ID))   # private per-task RDS dir
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# private manifest + summary dirs (never clobber locked outputs)
RUN_MANIFEST <- file.path(MANIFEST_DIR, "arc_runs", RUN_ID)
RUN_SUMMARY  <- file.path(RESULTS_DIR, "arc")
for (d in c(RUN_MANIFEST, RUN_SUMMARY)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- the 3 rungs (parametric / adaptive / aggressive) ----------------
RUNGS_R12 <- c("L1_param", "L3_adaptive", "L4_aggressive")
stopifnot(all(RUNGS_R12 %in% names(SL_LADDER)))

# ---- SMOKE mode: 100 reps, L4 only, 1 chunk (the spec's gate cell) ----
if (SMOKE) {
  N_REPS <- geti("SIM_N_REPS", 100L)   # smoke default 100 (env can override down)
  CHUNK  <- N_REPS
  RUNGS_R12 <- "L4_aggressive"         # heaviest rung; DEFF + PA-degradation gate
  cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
  cat(sprintf("[SMOKE] %s: %d reps, rung=%s, cores=%d\n", RUN_ID, N_REPS, RUNGS_R12[1], cores))
}

# ---- cell grid: rung x rep-chunk (single Design C scenario) ----------
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(RUNGS_R12, function(r)
  data.frame(rung = r, chunk = seq_len(n_chunks), stringsAsFactors = FALSE)))
stopifnot(task >= 1L, task <= nrow(grid))
job   <- grid[task, ]
rung  <- job$rung; chunk <- job$chunk
learners <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[%s task %d] DesignC rung=%s chunk=%d reps=%d-%d cores=%d\n",
            RUN_ID, task, rung, chunk, rep_lo, rep_hi, cores))
cat("  Design C knobs: "); print(unlist(designC_knobs()[c("sigma2_Y","sigma2_A","base_n0","base_m","alpha_strat")]))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
fn <- file.path(OUT, sprintf("r12_DesignC_%s_chunk%03d%s.rds", rung, chunk, if (SMOKE) "_SMOKE" else ""))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build the Design C population + truth ONCE (POP_SEED) -----------
pop <- make_population_designC(pop_seed = POP_SEED, truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0, abs(pop$truth$gap_super_census)))

# ---- one replication: Design C sample -> estimators -> deff ----------
# Mirrors run_sim.R::one_rep exactly, but draws with draw_sample_designC
# (bigger within-PSU n). The estimator call, EIF DEFF computation, and the
# carried check columns are IDENTICAL to the canonical driver so the
# per-rep schema (and hence aggregation) matches results/sim_full_summary.csv.
one_rep <- function(i, pop, learners) tryCatch({
  obs <- draw_sample_designC(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_estimators(obs, learners = learners)
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow, deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk across cores (parallel-safe L'Ecuyer streams) ------
cl <- makeCluster(cores)
on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl, c("CODE", ".this_dir", "RUN_ID"), envir = environment())
clusterEvalQ(cl, {
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  source(file.path(.this_dir, "design_c.R"))   # exposes draw_sample_designC on workers
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners)
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored (e.g. degenerate learner fit)
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps in this chunk failed", RUN_ID, task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (private; never clobbers locked sim outputs) ----
out <- list(run_id = RUN_ID, scenario = "DesignC", rung = rung, learners = learners,
            chunk = chunk, reps = reps, Psi = Psi, truth = pop$truth, params = pop$params,
            designC = designC_knobs(), results = res, diagnostics = diags)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- realized-DEFF report (the headline number reviewers asked for) ---
fa_rows <- res[res$method == "Fully-Aware", ]
cat(sprintf("[%s task %d] REALIZED deff_clust (FA eif): mean=%.2f  median=%.2f  range=[%.2f, %.2f]  | icc_eif mean=%.3f\n",
            RUN_ID, task, mean(fa_rows$deff), median(fa_rows$deff),
            min(fa_rows$deff), max(fa_rows$deff), mean(fa_rows$icc_eif)))

# ---- in-task summary (matches aggregate_sim.R column conventions) ----
# Written so a single smoke task is self-describing; the full run is
# aggregated across chunks by aggregate.R (this folder).
z_or_t <- function(df) qt(0.975, pmax(1, df))
summ <- do.call(rbind, by(res, res$method, function(d) {
  crit <- z_or_t(d$df)
  cov  <- mean(abs(d$b - Psi) <= crit * d$se)
  data.frame(scenario = "DesignC", rung = rung, method = d$method[1],
             n_reps = nrow(d), Psi = Psi,
             bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
             se_ratio = mean(d$se) / sd(d$b),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
             deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif),
             stringsAsFactors = FALSE)
}))
rownames(summ) <- NULL
ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
summ <- summ[order(match(summ$method, ord_m)), ]
np <- sapply(summ, is.numeric); sp <- summ; sp[np] <- round(sp[np], 4)
cat(sprintf("\n[%s task %d] per-task summary (rung=%s):\n", RUN_ID, task, rung))
print(sp, row.names = FALSE)

# write/append this task's summary to the private CSV
csv <- file.path(RUN_SUMMARY, sprintf("%s_summary%s.csv", RUN_ID, if (SMOKE) "_SMOKE" else ""))
write.table(sp, csv, sep = ",", row.names = FALSE,
            col.names = !file.exists(csv), append = file.exists(csv))
cat(sprintf("[%s task %d] summary -> %s\n", RUN_ID, task, csv))

# ---- reproducibility manifest (one per task; reviewer-proofing) ------
manifest <- list(
  run_id = RUN_ID, scenario = "DesignC", rung = rung, learners = learners,
  chunk = chunk, reps = reps, designC = designC_knobs(),
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),
  realized_deff = list(mean = mean(fa_rows$deff), median = median(fa_rows$deff),
                       min = min(fa_rows$deff), max = max(fa_rows$deff),
                       icc_eif = mean(fa_rows$icc_eif)),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  smoke = SMOKE
)
mf <- file.path(RUN_MANIFEST, sprintf("manifest_DesignC_%s_chunk%03d%s.rds",
                                      rung, chunk, if (SMOKE) "_SMOKE" else ""))
saveRDS(manifest, mf)
cat(sprintf("[%s task %d] manifest %s\n", RUN_ID, task, mf))
