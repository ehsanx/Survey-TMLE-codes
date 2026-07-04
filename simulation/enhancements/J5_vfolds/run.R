# =====================================================================
# J5_vfolds/run.R -- number-of-folds (V) sensitivity for the CF arm
#
# Closes the review finding: no V sensitivity; fold count is confounded with
# design geometry in the headline (V_cf = 5 everywhere in Design A).
#
# Cells (6): V in {2, 5, 10} x rung in {L2_smooth, L4_aggressive}, Design A
# (standard/complex), 1000 reps. V=5 reproduces the headline CF cell (an
# internal consistency anchor). FULL: 6 cells x 10 chunks -> --array=1-60.
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- nzchar(Sys.getenv("SMOKE"))
task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- if (SMOKE) 2L else geti("SIM_N_REPS", 1000L)
CHUNK  <- if (SMOKE) 2L else geti("SIM_CHUNK", 100L)
OUT    <- Sys.getenv("J_OUT", "sim_output/enhancements/J5_vfolds")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CELLS <- if (SMOKE) data.frame(V = 2L, rung = "L2_smooth", stringsAsFactors = FALSE) else
  expand.grid(V = c(2L, 5L, 10L), rung = c("L2_smooth", "L4_aggressive"),
              stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(CELLS)), function(i)
  cbind(CELLS[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
V <- grid$V[task]; rung <- grid$rung[task]; chunk <- grid$chunk[task]
learners <- SL_LADDER[[rung]]
reps <- ((chunk - 1L) * CHUNK + 1L):min(chunk * CHUNK, N_REPS)
cat(sprintf("[J5 task %d] V=%d rung=%s chunk=%d reps=%d-%d\n",
            task, V, rung, chunk, min(reps), max(reps)))

pop <- make_population("standard", model_type = "complex",
                       truth_M = if (SMOKE) 5e4L else 2e6L)
Psi <- pop$truth$psi

one_rep <- function(i, pop, learners, V) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- tryCatch(run_estimators(obs, learners = learners, V_cf = V),
                  error = function(e) NULL)
  if (is.null(est)) return(NULL)
  cbind(rep = i, V_cf = V, est$results)
}

cl <- makeCluster(cores); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "R")
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 920000L + task)
res <- do.call(rbind, Filter(Negate(is.null),
         parLapply(cl, reps, one_rep, pop = pop, learners = learners, V = V)))

out <- list(run = "J5_vfolds", V_cf = V, rung = rung, learners = learners,
            chunk = chunk, reps = reps, Psi = Psi, truth = pop$truth,
            params = pop$params, results = res)
fn <- file.path(OUT, sprintf("j5_V%02d_%s_chunk%03d.rds", V, rung, chunk))
saveRDS(out, fn)
cat(sprintf("[J5 task %d] saved %s (%d rows)\n", task, fn, nrow(res)))
if (SMOKE) cat("J5 SMOKE OK\n")
