# =====================================================================
# J4_overlap/run.R -- limited-overlap / near-positivity stress
#
# Closes the review finding: no simulated near-positivity-violation scenario
# with a KNOWN truth (webF_overlap is an NHANES diagnostic, not a simulation).
#
# Knob: alpha_g, the confounding strength in the treatment model
# (eta = alpha_g * (C1+C2+C3+C4) + uA; the intercept is re-solved so the
# marginal P(A=1) stays at the target, so increasing alpha_g stretches the
# propensity toward {0,1} WITHOUT changing the outcome law -- the
# super-population truth Psi is invariant to alpha_g by construction).
#   base  alpha_g = log(1.3)  (the headline value; moderate overlap)
#   mid   alpha_g = log(2.5)  (strained overlap)
#   high  alpha_g = log(5.0)  (near-positivity violation)
#
# Cells (5): L2_smooth x {base, mid, high} + L4_aggressive x {base, high}.
# FULL: 5 cells x 10 chunks -> --array=1-50.
# Watch: CF coverage + the g-floor interaction (g_oof_bound = 0.05 default);
# per-rep fitted-g diagnostics (min/max, epsilon) ship in the diag rows.
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
OUT    <- Sys.getenv("J_OUT", "sim_output/enhancements/J4_overlap")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ALPHAS <- c(base = log(1.3), mid = log(2.5), high = log(5.0))
CELLS  <- if (SMOKE) data.frame(rung = "L2_smooth", alpha = "base", stringsAsFactors = FALSE) else
  rbind(expand.grid(rung = "L2_smooth",     alpha = c("base", "mid", "high"),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE),
        expand.grid(rung = "L4_aggressive", alpha = c("base", "high"),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE))
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(CELLS)), function(i)
  cbind(CELLS[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
rung <- grid$rung[task]; atag <- grid$alpha[task]; chunk <- grid$chunk[task]
learners <- SL_LADDER[[rung]]
reps <- ((chunk - 1L) * CHUNK + 1L):min(chunk * CHUNK, N_REPS)
cat(sprintf("[J4 task %d] rung=%s alpha=%s(%.3f) chunk=%d reps=%d-%d\n",
            task, rung, atag, ALPHAS[[atag]], chunk, min(reps), max(reps)))

pop <- make_population("standard", model_type = "complex", alpha_g = ALPHAS[[atag]],
                       truth_M = if (SMOKE) 5e4L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[J4 task %d] Psi=%.5f (alpha_g invariance check: base headline ~0.209)\n", task, Psi))

one_rep <- function(i, pop, learners) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- tryCatch(run_estimators(obs, learners = learners), error = function(e) NULL)
  if (is.null(est)) return(NULL)
  list(results = cbind(rep = i, est$results),
       diag    = cbind(rep = i, est$diagnostics$drow))
}

cl <- makeCluster(cores); on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "R")
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 910000L + task)
reps_out <- Filter(Negate(is.null), parLapply(cl, reps, one_rep, pop = pop, learners = learners))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

out <- list(run = "J4_overlap", rung = rung, alpha_tag = atag, alpha_g = ALPHAS[[atag]],
            learners = learners, chunk = chunk, reps = reps, Psi = Psi,
            truth = pop$truth, params = pop$params, results = res, diagnostics = diags)
fn <- file.path(OUT, sprintf("j4_%s_%s_chunk%03d.rds", rung, atag, chunk))
saveRDS(out, fn)
cat(sprintf("[J4 task %d] saved %s (%d rows)\n", task, fn, nrow(res)))
if (SMOKE) cat("J4 SMOKE OK\n")
