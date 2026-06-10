# =====================================================================
# run_sim.R  —  cluster-ready simulation driver (SLURM array job)
#
# One SLURM array task = one (scenario x library-rung x rep-chunk) cell.
# Each task builds the fixed population ONCE, computes the truth ONCE, runs its
# chunk of replications across the cores SLURM gave it, and writes ONE RDS to
# the output dir. Aggregate afterwards with aggregate_sim.R.
#
# Reps are seeded by a GLOBAL rep index (chunk offset) so chunks never collide
# and the whole run is reproducible regardless of how it is split across tasks.
#
# Env vars (with local-test fallbacks):
#   SLURM_ARRAY_TASK_ID   1-based task index               (default 1)
#   SLURM_CPUS_PER_TASK   cores for this task              (default detectCores-1)
#   SIM_N_REPS            total reps per (scenario,rung)   (default 1000)
#   SIM_CHUNK             reps per array task              (default 100)
#   SIM_OUT              output dir for per-task RDS       (default DATA_INTERMEDIATE)
#   SIM_CONFIG           "ladder" (complex, full ladder) | "test" (subset)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "codes")            # set SIM_CODE to the code dir on ARC
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
OUT    <- Sys.getenv("SIM_OUT", DATA_INTERMEDIATE)
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x rung (ladder runs on the realistic 'complex' DGP) ----
RUNGS <- names(SL_LADDER)                            # trimmed in learners.R to the chosen ladder
cells <- expand.grid(scenario = c("standard", "R1"), rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job <- grid[task, ]
scenario <- job$scenario; rung <- job$rung; chunk <- job$chunk
learners <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[task %d] scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, cores))

# ---- build population + truth ONCE (deterministic from POP_SEED) ----
pop <- make_population(scenario, model_type = "complex", truth_M = 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0, abs(pop$truth$gap_super_census)))

one_rep <- function(i, pop, learners) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_estimators(obs, learners = learners)
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow, deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}

cl <- makeCluster(cores)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "codes")
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners)
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output: estimates + per-rep diagnostics + truth (full detail) ----
out <- list(scenario = scenario, rung = rung, learners = learners, chunk = chunk,
            reps = reps, Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags)
fn <- file.path(OUT, sprintf("sim_%s_%s_chunk%03d.rds", scenario, rung, chunk))
saveRDS(out, fn)
cat(sprintf("[task %d] saved %s  (%d est rows, %d diag rows)\n", task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; race-free; reviewer-proofing) ----
manifest <- list(
  scenario = scenario, rung = rung, learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),                 # realized ICC of Y/A/C, mean A/Y
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA),
  sysname = Sys.info()[["nodename"]]
)
mf <- file.path(MANIFEST_DIR, sprintf("manifest_%s_%s_chunk%03d.rds", scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[task %d] manifest %s\n", task, mf))
