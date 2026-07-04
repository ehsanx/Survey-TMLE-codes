# =====================================================================
# run.R  —  R07_informative driver (SLURM array job)
#
# GOAL: stress-test the de-weighted out-of-fold (OOF) nuisance choice in the
# Fully-Aware-CF arm under INFORMATIVE-sampling-BEYOND-C. For each rho in the
# sweep we draw informatively-selected samples (dgp_infsamp.R) and run TWO
# cross-fitted arms (cf_arms.R) that differ ONLY in whether the OOF fits are
# weighted:
#     Fully-Aware-CF-unwt   (paper default; de-weighted OOF)
#     Fully-Aware-CF-wt     (weighted OOF; the defended alternative)
# Reported: BIAS and COVERAGE curves for both arms vs rho.
#
# One SLURM array task = one (rho x rung x rep-chunk) cell. Each task builds the
# fixed population ONCE (canonical make_population), computes the truth ONCE,
# attaches the per-PSU selection driver uS(rho) ONCE, runs its chunk of reps over
# the cores SLURM gave it, and writes ONE RDS. Aggregate across tasks at the end.
#
# Reuses the canonical engine by SOURCING it; defines NO new estimator/aggregation
# logic beyond the two helper files in THIS folder.
#
# Env vars (with local-test fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores for this task           (default detectCores-1)
#   SIM_N_REPS           total reps per (rho,rung)      (default 1000)
#   SIM_CHUNK            reps per array task            (default 100)
#   R07_OUT              per-task RDS dir   (default sim_output/arc_runs/R07_informative)
#   SIM_CODE             canonical code dir (default R)
#   R07_DIR              this run's folder  (default simulation/enhancements/R07_informative)
#   SMOKE=1              tiny run: 20 reps, rho in {0,0.9}, rung L2_smooth only
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))   # exposes .sl/.eif_from_tmle/.se_des/make_cf_folds
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

# this run's helpers (the ONLY new logic)
R07_DIR <- Sys.getenv("R07_DIR", file.path(CODE, "arc_runs", "R07_informative"))
source(file.path(R07_DIR, "dgp_infsamp.R"))
source(file.path(R07_DIR, "cf_arms.R"))

suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

# ---- run knobs --------------------------------------------------------------
SCENARIO <- "standard"                      # standard design per spec
LAMBDA   <- 0.8                             # informative-selection strength (weight spread);
                                            # tuned: mean(sumw/N)~1, w_cv~0.9 (realistic survey
                                            # weight spread), and the unweighted in-sample E[Y]
                                            # drift ~DOUBLES from rho=0 to rho=0.9 (selection bites)
RHOS     <- c(0, 0.3, 0.6, 0.9)             # informativeness sweep
RUNGS    <- c("L2_smooth", "L3_adaptive")   # the two libraries in the spec

if (SMOKE) {                                # tiny, ~1-3 min, one population build
  RHOS   <- c(0, 0.9)
  RUNGS  <- c("L2_smooth")
  N_REPS <- 20L; CHUNK <- 20L
} else {
  N_REPS <- geti("SIM_N_REPS", 1000L)
  CHUNK  <- geti("SIM_CHUNK", 100L)
}

task  <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
OUT   <- Sys.getenv("R07_OUT", file.path(DATA_ROOT, "arc_runs", "R07_informative"))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: rho x rung x rep-chunk --------------------------------------
cells <- expand.grid(rho = RHOS, rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job <- grid[task, ]
rho <- as.numeric(job$rho); rung <- job$rung; chunk <- as.integer(job$chunk)
learners <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[task %d] R07_informative rho=%.2f rung=%s chunk=%d reps=%d-%d cores=%d (SMOKE=%s)\n",
            task, rho, rung, chunk, rep_lo, rep_hi, cores, SMOKE))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
fn <- file.path(OUT, sprintf("r07_rho%03d_%s_chunk%03d.rds",
                             round(rho * 100), rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (canonical; deterministic from POP_SEED) --
pop <- make_population(SCENARIO, model_type = "complex", truth_M = 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f\n",
            task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0))

# ---- attach the per-PSU selection driver uS(rho) ONCE -----------------------
# deterministic in (POP_SEED, rho) so every chunk for this rho uses the SAME uS;
# reps differ only through the PPS draw (sample_seed) and stage-2 SRS.
uS_tab <- attach_uS(pop, rho = rho, seed = POP_SEED)
cat(sprintf("[task %d] uS attached: cor(uS, uY+uA)=%.3f over %d PSUs\n",
            task, suppressWarnings(stats::cor(uS_tab$uS, uS_tab$r_psu)), nrow(uS_tab)))

# ---- one replication: informative sample -> both CF arms -> diagnostics -----
one_rep <- function(i, pop, learners, rho, lambda, uS_tab) tryCatch({
  obs <- draw_sample_infsamp(pop, sample_seed = SAMPLE_SEED_BASE + i,
                             model_type = "complex", rho = rho, lambda = lambda,
                             uS_tab = uS_tab)
  pair <- run_cf_pair(obs, learners = learners)
  ch   <- attr(obs, "checks")
  dd   <- deff_clust(pair$eif_unwt, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, rho = rho, rung = NA_character_, pair$results,
                    deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, rho = rho, pair$diag,
                    deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- parallel over reps (mirror run_sim.R worker setup) ---------------------
run_seq <- function() lapply(reps, one_rep, pop = pop, learners = learners,
                             rho = rho, lambda = LAMBDA, uS_tab = uS_tab)
reps_out <- if (cores > 1L && length(reps) > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("R07_DIR"), envir = environment())
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    R07_DIR <- Sys.getenv("R07_DIR", file.path(CODE, "arc_runs", "R07_informative"))
    source(file.path(R07_DIR, "dgp_infsamp.R")); source(file.path(R07_DIR, "cf_arms.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
  parLapply(cl, reps, one_rep, pop = pop, learners = learners,
            rho = rho, lambda = LAMBDA, uS_tab = uS_tab)
} else {
  set.seed(SAMPLE_SEED_BASE + task)
  run_seq()
}

reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored (e.g. NULL-learner crash)
if (!length(reps_out)) stop(sprintf("[task %d] all %d reps in this chunk failed", task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
res$rung <- rung                            # fill the rung label (NA placeholder above)
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output --------------------------------------------------------
out <- list(run = "R07_informative", scenario = SCENARIO, rho = rho, rung = rung,
            lambda = LAMBDA, learners = learners, chunk = chunk, reps = reps,
            Psi = Psi, truth = pop$truth, params = pop$params,
            uS_cor = suppressWarnings(stats::cor(uS_tab$uS, uS_tab$r_psu)),
            results = res, diagnostics = diags)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s  (%d est rows, %d diag rows)\n",
            task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; reviewer-proofing) -------------
manifest <- list(
  run = "R07_informative", scenario = SCENARIO, rho = rho, rung = rung,
  lambda = LAMBDA, learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  uS_cor = out$uS_cor,
  pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]]
)
mdir <- file.path(OUT, "manifest")
if (!dir.exists(mdir)) dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
mf <- file.path(mdir, sprintf("manifest_rho%03d_%s_chunk%03d.rds",
                              round(rho * 100), rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[task %d] manifest %s\n", task, mf))
