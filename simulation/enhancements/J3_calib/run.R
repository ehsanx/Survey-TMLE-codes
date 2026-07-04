# =====================================================================
# J3_calib/run.R -- calibrated (raked) weights vs design weights
#
# Closes the review finding: the NHANES application treats CALIBRATED weights
# as inverse-inclusion-probability weights, but no simulation checks that the
# CF-TMLE + Eq-8 linearization variance stays calibrated under raking.
#
# Design: headline Design A population (standard/complex, POP_SEED fixed).
# Per rep: draw the two-stage sample (design weights = true HT weights), then
# RAKE those weights to two KNOWN population margins (stratum counts x L1
# tertile counts -- the practice-analogue of post-stratifying to demographic
# controls). Run the full 5-arm suite TWICE: design weights vs raked weights.
# If Eq-8 + CF-TMLE is robust to calibration (the paper's implicit NHANES
# assumption), raked-weight coverage should match design-weight coverage
# (raking may mildly shrink SEs -- classical calibration efficiency).
#
# Cells: rung in {L1_param, L2_smooth}. FULL: 2 rungs x 10 chunks -> --array=1-20
# Env: SLURM_ARRAY_TASK_ID, SIM_N_REPS (1000), SIM_CHUNK (100), SMOKE=1 for a
#      2-rep single-cell check. Same seed conventions as run_sim.R.
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
OUT    <- Sys.getenv("J_OUT", "sim_output/enhancements/J3_calib")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RUNGS <- if (SMOKE) "L1_param" else c("L1_param", "L2_smooth")
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- expand.grid(rung = RUNGS, chunk = seq_len(n_chunks),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
stopifnot(task >= 1L, task <= nrow(grid))
rung  <- grid$rung[task]; chunk <- grid$chunk[task]
learners <- SL_LADDER[[rung]]
reps <- ((chunk - 1L) * CHUNK + 1L):min(chunk * CHUNK, N_REPS)
cat(sprintf("[J3 task %d] rung=%s chunk=%d reps=%d-%d cores=%d\n",
            task, rung, chunk, min(reps), max(reps), cores))

pop <- make_population("standard", model_type = "complex",
                       truth_M = if (SMOKE) 5e4L else 2e6L)
Psi <- pop$truth$psi

# ---- KNOWN population margins (computable because the population is known) --
# Margin 1: stratum counts. Margin 2: tertiles of the OBSERVED covariate L1
# (apply_L is deterministic given the C's, so population L1 is available).
pop_L  <- apply_L(pop$pop$C1, pop$pop$C2, pop$pop$C3, pop$pop$C4, "complex")
l1_cut <- quantile(pop_L$L1, c(1/3, 2/3))
pop_l1ter <- cut(pop_L$L1, c(-Inf, l1_cut, Inf), labels = c("T1", "T2", "T3"))
MARG_STRATA <- as.data.frame(table(strata = pop$pop$strata))
MARG_L1TER  <- as.data.frame(table(L1ter  = pop_l1ter))

rake_weights <- function(obs) {
  obs$L1ter <- cut(obs$L1, c(-Inf, l1_cut, Inf), labels = c("T1", "T2", "T3"))
  des <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight,
                   data = obs, nest = FALSE)
  rk  <- tryCatch(
    rake(des, sample.margins = list(~strata, ~L1ter),
         population.margins = list(MARG_STRATA, MARG_L1TER),
         control = list(maxit = 100, epsilon = 1e-7)),
    error = function(e) NULL)
  if (is.null(rk)) return(NULL)
  as.numeric(weights(rk))
}

one_rep <- function(i, pop, learners) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  out <- list()
  for (scheme in c("design", "raked")) {
    o <- obs
    if (scheme == "raked") {
      wcal <- rake_weights(obs)
      if (is.null(wcal)) next                      # raking failed: drop scheme this rep
      o$weight <- wcal
    }
    est <- tryCatch(run_estimators(o, learners = learners), error = function(e) NULL)
    if (is.null(est)) next
    out[[scheme]] <- cbind(rep = i, weight_scheme = scheme, est$results,
                           w_cv = sd(o$weight) / mean(o$weight))
  }
  do.call(rbind, out)
}

cl <- makeCluster(cores); on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl, c("rake_weights", "l1_cut", "MARG_STRATA", "MARG_L1TER",
                    "SAMPLE_SEED_BASE"), envir = environment())
clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "R")
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 900000L + task)
res <- do.call(rbind, parLapply(cl, reps, one_rep, pop = pop, learners = learners))

out <- list(run = "J3_calib", rung = rung, learners = learners, chunk = chunk,
            reps = reps, Psi = Psi, truth = pop$truth, params = pop$params,
            margins = list(strata = MARG_STRATA, L1ter = MARG_L1TER, l1_cut = l1_cut),
            results = res)
fn <- file.path(OUT, sprintf("j3_%s_chunk%03d.rds", rung, chunk))
saveRDS(out, fn)
cat(sprintf("[J3 task %d] saved %s (%d rows)\n", task, fn, nrow(res)))
if (SMOKE) cat("J3 SMOKE OK\n")
