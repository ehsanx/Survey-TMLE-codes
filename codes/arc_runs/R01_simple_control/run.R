# =====================================================================
# R01_simple_control / run.R
#   Pipeline-correctness CONTROL for the survey-TMLE cross-fit bias story.
#
# WHY THIS RUN EXISTS
#   The headline run (codes/run_sim.R) builds the population with
#   model_type = 'complex' (Kang-Schafer transforms in dgp.R::apply_L), so the
#   parametric/smooth learners on L1-L3 are MIS-specified -- even L1's GLM sees
#   the transformed L1..L4, not the linear C1..C4 that actually drive A and Y.
#   That misspecification is the most plausible source of the +0.012..+0.021
#   cross-fit (and single-fit) bias seen at L1-L3 in results/sim_full_summary.csv.
#
#   This control re-runs the SAME engine, SAME five arms, but on the
#   model_type = 'simple' DGP (apply_L returns C1..C4 unchanged), where the
#   parametric learner IS correctly specified, so bias should collapse to ~0 and
#   coverage to nominal. We run L1 (+ L2 for context) under 'simple', and ALSO
#   re-run L1 under 'complex' so the matched complex/simple pair sits side by
#   side in one CSV. Every output row is tagged with a model_type column (the
#   canonical tables lack it).
#
# NON-DESTRUCTIVE: this driver SOURCES codes/*.R unchanged. The only new logic
#   is (a) looping model_type and building the population with the MATCHING
#   model_type, and (b) drawing samples with that model_type. No canonical file
#   is edited. Outputs go to a private arc_runs/ subtree; nothing clobbers the
#   locked results.
#
# ENV VARS (local-test fallbacks in brackets):
#   SLURM_ARRAY_TASK_ID  1-based task index                       [1]
#   SLURM_CPUS_PER_TASK  cores for this task                      [detectCores-1]
#   SIM_N_REPS           total reps per (scenario,model_type,rung)[1000]
#   SIM_CHUNK            reps per array task                      [100]
#   SIM_CODE             canonical code dir                       [codes]
#   ARC_OUT              per-task RDS dir   [sim_output/arc_runs/R01_simple_control]
#   SMOKE                "1" -> tiny run (20 reps, 1 cell) for local validation
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "codes")          # set SIM_CODE to the code dir on ARC
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

RUN_ID <- "R01_simple_control"
SMOKE  <- Sys.getenv("SMOKE", "0") == "1"

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)

# default per-task output dir lives UNDER the data root's sim_output tree, in a
# private arc_runs subfolder so it never collides with the locked intermediate/.
ARC_OUT <- Sys.getenv("ARC_OUT", file.path(dirname(DATA_INTERMEDIATE), "arc_runs", RUN_ID))
if (!dir.exists(ARC_OUT)) dir.create(ARC_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: (scenario) x (model_type, rung) ----------------------------
# The (model_type, rung) pairs are the substance of this control:
#   simple  / L1_param   correctly-specified parametric -> EXPECT bias ~0, cov ~.95
#   simple  / L2_smooth  correctly-specified smooth     -> context
#   complex / L1_param   misspecified parametric        -> reproduces headline bias
# (No complex/L2 here; the headline run already has it -- this run is the cheap
#  matched control, not a re-do of the full ladder.)
mt_rung <- data.frame(
  model_type = c("simple",   "simple",    "complex"),
  rung       = c("L1_param", "L2_smooth", "L1_param"),
  stringsAsFactors = FALSE
)
cells <- do.call(rbind, lapply(c("standard", "R1"), function(sc)
  cbind(scenario = sc, mt_rung, stringsAsFactors = FALSE)))

if (SMOKE) {
  # tiny: one cell (standard, simple, L1), 20 reps, single chunk -> 1-3 min.
  cells  <- cells[cells$scenario == "standard" &
                  cells$model_type == "simple" & cells$rung == "L1_param", , drop = FALSE]
  N_REPS <- 20L; CHUNK <- 20L
  cores  <- min(cores, 4L)
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job        <- grid[task, ]
scenario   <- job$scenario
model_type <- job$model_type
rung       <- job$rung
chunk      <- job$chunk
learners   <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[%s task %d] scenario=%s model_type=%s rung=%s chunk=%d reps=%d-%d cores=%d%s\n",
            RUN_ID, task, scenario, model_type, rung, chunk, rep_lo, rep_hi, cores,
            if (SMOKE) "  [SMOKE]" else ""))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
fn <- file.path(ARC_OUT, sprintf("r01_%s_%s_%s_chunk%03d.rds", scenario, model_type, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE, with the MATCHING model_type ------------
# This is the one substantive difference from run_sim.R, which hardcodes
# model_type = 'complex'. Under 'simple', apply_L is the identity so the same
# C1..C4 that drive A and Y are what the learners see -> correct specification.
pop <- make_population(scenario, model_type = model_type, truth_M = 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0,
            abs(pop$truth$gap_super_census)))

# one_rep mirrors run_sim.R exactly, except draw_sample uses the cell model_type.
one_rep <- function(i, pop, learners, model_type) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = model_type)
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

cl <- makeCluster(cores)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "codes")
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
# seed stream keyed by task (same convention as run_sim.R; distinct per cell/chunk)
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners,
                      model_type = model_type)
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored (e.g. degenerate learner fit)
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps failed", RUN_ID, task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# tag both tables with model_type (canonical tables lack it) ------------------
res$model_type   <- model_type
diags$model_type <- model_type

# ---- per-task output (private; never touches the locked intermediate/) ------
out <- list(run_id = RUN_ID, scenario = scenario, model_type = model_type,
            rung = rung, learners = learners, chunk = chunk, reps = reps,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; reviewer-proofing) -------------
manifest <- list(
  run_id = RUN_ID, scenario = scenario, model_type = model_type, rung = rung,
  learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]]
)
mf <- file.path(ARC_OUT, sprintf("manifest_%s_%s_%s_chunk%03d.rds",
                                 scenario, model_type, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[%s task %d] manifest %s\n", RUN_ID, task, mf))
