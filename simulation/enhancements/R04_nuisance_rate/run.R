# =====================================================================
# run.R  —  R04_nuisance_rate driver  (NON-DESTRUCTIVE ARC enhancement run)
#
# GOAL: empirically MEASURE the cross-fit nuisance L2 errors per rung (L1-L4) and
# scenario, to DEMONSTRATE the TMLE product-rate condition
#   ||Qhat - Q0|| * ||ghat - g0||  =  o_P( m^{-1/2} )
# and show it is MET on the smooth/parametric rungs (L1-L3) but deliberately
# VIOLATED on L4 (interpolating deep RF). The diagnostic object is
#   product * sqrt(m)  ->  shrinks across L1->L3, blows up / stays large at L4.
#
# This driver SOURCES the canonical engine read-only and reuses run_estimators'
# building blocks via nuisance_rate_helpers.R. It NEVER edits R/*.R and writes
# only into sim_output/arc_runs/R04_nuisance_rate/ + results/.
#
# One SLURM array task = one (scenario x rung x rep-chunk) cell, mirroring run_sim.R.
#
# Env vars (with local fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index           (default 1)
#   SLURM_CPUS_PER_TASK  cores                         (default detectCores-1)
#   SIM_N_REPS           total reps per (scenario,rung)(default 200)
#   SIM_CHUNK            reps per array task           (default 50)
#   SIM_CODE             canonical code dir            (default R)
#   R04_OUT              per-task RDS out dir          (default DATA_ROOT/arc_runs/R04_nuisance_rate)
#   SMOKE                "1" -> tiny 1-cell, few-rep local validation
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))     # exposes .sl(), make_cf_folds(), etc.
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))       # SL_LADDER, custom RF wrappers
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# run-local helpers (truth extraction + L2 instrumentation); resolve next to this file
.this_dir <- Sys.getenv("R04_DIR",
  file.path(REPO_ROOT, "simulation", "enhancements", "R04_nuisance_rate"))
source(file.path(.this_dir, "nuisance_rate_helpers.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 200L)
CHUNK  <- geti("SIM_CHUNK", 50L)

# ---- output dirs (NEVER clobber locked results) ----
OUT <- Sys.getenv("R04_OUT", file.path(DATA_ROOT, "arc_runs", "R04_nuisance_rate"))
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
for (d in c(OUT, ARC_RESULTS)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x rung x chunk  (same shape as run_sim.R) ----
RUNGS <- names(SL_LADDER)
cells <- expand.grid(scenario = c("standard", "R1"), rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

if (SMOKE) {
  # tiny validation: ONE cheap cell (standard x L1_param), few reps, serial.
  cells  <- cells[cells$scenario == "standard" & cells$rung == "L1_param", , drop = FALSE]
  N_REPS <- 20L; CHUNK <- 20L; cores <- 1L; task <- 1L
  cat("[SMOKE] standard x L1_param, 20 reps, 1 core\n")
}

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

cat(sprintf("[task %d] scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d learners=%s\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, cores, paste(learners, collapse = "+")))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr04_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ----
pop <- make_population(scenario, model_type = "complex",
                       truth_M = if (SMOKE) 2e5L else 2e6L)
cat(sprintf("[task %d] Psi=%.5f Psi_N=%.5f gap=%.6f\n",
            task, pop$truth$psi, pop$truth$psi_N_Q0, abs(pop$truth$gap_super_census)))

# ---- run the chunk's reps (parallel over reps, like run_sim.R) ----
run_one <- function(i) tryCatch(
  one_rep_rate(i, pop, learners, rung, scenario, V_cf = 5L, g_oof_bound = 0.05),
  error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    .this_dir <- Sys.getenv("R04_DIR",
      file.path(REPO_ROOT, "simulation", "enhancements", "R04_nuisance_rate"))
    source(file.path(.this_dir, "nuisance_rate_helpers.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
  clusterExport(cl, c("pop", "learners", "rung", "scenario"), envir = environment())  # FIX: workers need `pop`
  rows <- parLapply(cl, reps, run_one)
} else {
  set.seed(SAMPLE_SEED_BASE + task)
  rows <- lapply(reps, run_one)
}
rows <- Filter(Negate(is.null), rows)   # drop reps that errored (e.g. degenerate learner fit)
if (!length(rows)) stop(sprintf("[task %d] all %d reps in this chunk failed", task, length(reps)))
res <- do.call(rbind, rows)

# ---- per-task RDS (full detail; never clobbers locked outputs) ----
out <- list(run = "R04_nuisance_rate", scenario = scenario, rung = rung,
            learners = learners, chunk = chunk, reps = reps,
            Psi = pop$truth$psi, truth = pop$truth, params = pop$params,
            per_rep = res)
saveRDS(out, fn)   # `tag`/`fn` were computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s  (%d rep rows)\n", task, fn, nrow(res)))

# ---- reproducibility manifest (one per task; mirrors run_sim.R) ----
manifest <- list(
  run = "R04_nuisance_rate",
  scenario = scenario, rung = rung, learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  truth_join = unique(res$truth_join))
mf <- file.path(OUT, sprintf("%smanifest_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[task %d] manifest %s  (truth_join=%s)\n",
            task, mf, paste(unique(res$truth_join), collapse = ",")))

# ---- SMOKE: aggregate inline + print the rate table so the gate is visible ----
if (SMOKE) {
  cat("\n[SMOKE] per-rung nuisance-rate summary (this single cell):\n")
  agg <- with(res, data.frame(
    scenario = scenario[1], rung = rung[1], n_reps = nrow(res),
    mean_eQ = mean(eQ_int), mean_eg = mean(eg_int),
    mean_prod = mean(prod_int),
    mean_prod_sqrtm = mean(prod_int_sqrtm),
    mean_prod_sqrtn = mean(prod_int_sqrtn),
    mean_m_psu = mean(m_psu)))
  print(format(agg, digits = 4), row.names = FALSE)
  cat("\n[SMOKE] OK if columns are finite & truth_join=OK. ",
      "Full-run decision rule is in NOTES.md.\n")
}
