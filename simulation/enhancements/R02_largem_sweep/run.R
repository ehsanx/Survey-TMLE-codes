# =====================================================================
# run.R  —  R02_largem_sweep driver (SLURM array job, NON-DESTRUCTIVE)
#
# GOAL: settle whether the L4 cross-fit over-coverage (se_ratio ~ 1.44,
# coverage 0.985-0.992 in results/sim_full_summary.csv) shrinks as the number
# of sampled PSUs m grows (finite-m, conservative-consistent per Theorem 2) or
# is structural. Hold DGP + library FIXED at the 'standard' scenario and sweep
# m = sampled-PSUs-per-stratum (base_m) across {6, 12, 20, 30}
# (= baseline, 2x, ~3.3x, 5x). 1000 reps per m.
#
# Two library rungs per m:
#   * L4_aggressive = SL.ranger.deep   (the over-covering deep-RF rung)
#   * RF_shallow     = SL.ranger.t1     (default-depth RF contrast; its
#                                        se_ratio should approach 1 faster)
# We ALSO track bias at L1 (SL.glm) and L3 (SL_LADDER$L3_adaptive) across m to
# test the O(1/m) Hajek-ratio explanation for the R01 smooth-rung bias -- these
# two extra rungs are CHEAP, so they ride along in every task.
#
# The driver SOURCES the canonical engine (config/dgp/estimators/diagnostics/
# learners) and a run-local helpers.R. It NEVER edits R/*.R. The m knob is
# draw_sample()'s existing `base_m` arg (see helpers.R header).
#
# OUTPUTS (do NOT clobber locked results):
#   per-task RDS  -> sim_output/arc_runs/R02_largem_sweep/
#   manifest RDS  -> sim_output/arc_runs/R02_largem_sweep/manifest/
#   summary CSV   -> results/R02_largem_sweep_summary.csv (written by
#                    aggregate.R after all tasks finish; this driver writes
#                    only per-task RDS).
#
# Array layout: task = (m-level) x (rep-chunk).  ALL rungs run inside each task
# (so a task is one m at one chunk, across all 4 rungs). With 4 m-levels and
# N_REPS/CHUNK chunks, set --array=1-(4*n_chunks).
#
# Env vars (with local-test fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores for this task           (default detectCores-1)
#   SIM_N_REPS           reps per (m, rung)            (default 1000)
#   SIM_CHUNK            reps per array task           (default 100)
#   SIM_CODE             canonical code dir            (default R)
#   SMOKE                "1" -> tiny run (20 reps, largest m, L4+shallow only)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# run-local helpers (one_rep_m forwards base_m; summarise_sweep matches cols)
RUN_DIR <- file.path(REPO_ROOT, "simulation", "enhancements", "R02_largem_sweep")
# fall back to the directory of this script if REPO_ROOT layout differs
if (!file.exists(file.path(RUN_DIR, "helpers.R"))) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) RUN_DIR <- dirname(normalizePath(f))
}
source(file.path(RUN_DIR, "helpers.R"))

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE  <- identical(Sys.getenv("SMOKE"), "1")
task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)

# ---- run identity + output dirs (NON-DESTRUCTIVE; never touch locked dirs) ----
RUN_ID  <- "R02_largem_sweep"
OUT     <- file.path(DATA_ROOT, "arc_runs", RUN_ID)         # sim_output/arc_runs/<id>
MAN_OUT <- file.path(OUT, "manifest")
for (d in c(OUT, MAN_OUT)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

SCENARIO <- "standard"   # FIXED per spec (the L4 over-coverage is on standard)

# ---- the m sweep: sampled-PSUs-per-stratum (base_m). standard J/stratum=60. ---
# {baseline, 2x, ~3.3x, 5x}; baseline base_m=6 -> m_total=H*6=60 PSUs.
M_LEVELS <- c(6L, 12L, 20L, 30L)

# ---- the rungs run in EVERY task (L4 + shallow-RF contrast + L1/L3 bias track) -
# Names are the LABELS that land in the summary 'rung' column.
RUNGS <- list(
  L4_aggressive = SL_LADDER$L4_aggressive,   # "SL.ranger.deep"  (the over-coverer)
  RF_shallow    = "SL.ranger.t1",            # default-depth RF contrast
  L1_param      = SL_LADDER$L1_param,        # "SL.glm"          (bias track)
  L3_adaptive   = SL_LADDER$L3_adaptive      # ensemble + RF     (bias track)
)

# ---- SMOKE: fast pipeline validation (~3-4 min, validated on R 4.5) ----------
# Goal of the smoke is to validate the PIPELINE end-to-end at the largest m
# (correct df_design = H*(base_m-1) = 290, finite b/se, all columns present),
# NOT to estimate the L4 se_ratio (that needs 1000 reps on ARC). The default
# smoke rung is therefore the CHEAP L1_param (GLM in every arm); 20 reps at m=30
# on 2 cores took ~3.7 min locally. To eyeball the heavy deep-RF decision cell
# locally instead, set SMOKE_RUNG=L4_aggressive (expect several minutes PER REP
# -- this is the cell that gets a 6h ARC walltime; keep cores small to avoid OOM).
TRUTH_M <- 2e6L
if (SMOKE) {
  M_LEVELS <- c(30L)                                   # largest m (decision cell)
  sm_rung  <- Sys.getenv("SMOKE_RUNG", "L1_param")     # cheap GLM by default
  if (!sm_rung %in% names(RUNGS)) stop("SMOKE_RUNG must be one of: ",
                                       paste(names(RUNGS), collapse = ", "))
  RUNGS    <- RUNGS[sm_rung]
  N_REPS   <- if (identical(sm_rung, "L1_param")) 20L else 4L
  CHUNK    <- N_REPS
  cores    <- min(cores, 2L)                           # cap workers -> low RAM
  TRUTH_M  <- 2e5L                                      # faster population build
  cat(sprintf("[SMOKE] m=30, rung=%s, %d reps, cores=%d, truth_M=%d\n",
              sm_rung, N_REPS, cores, TRUTH_M))
}

# ---- build the (m-level x rep-chunk) task grid -------------------------------
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_along(M_LEVELS), function(mi)
  data.frame(m_idx = mi, base_m = M_LEVELS[mi], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job    <- grid[task, ]
base_m <- job$base_m
chunk  <- job$chunk
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[task %d] scenario=%s base_m=%d (m_total=%d) chunk=%d reps=%d-%d rungs=%s cores=%d\n",
            task, SCENARIO, base_m, 10L * base_m, chunk, rep_lo, rep_hi,
            paste(names(RUNGS), collapse = ","), cores))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%ssim_m%02d_chunk%03d.rds", tag, base_m, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) -------------
# Population geometry is FIXED (standard); only the SAMPLING base_m varies, so
# one population serves every m-level. H = pop$params$H = 10 for standard.
pop <- make_population(SCENARIO, model_type = "complex", truth_M = TRUTH_M)
Psi <- pop$truth$psi
Hstr <- pop$params$H
cat(sprintf("[task %d] Psi=%.5f (mc se %.6f) H=%d J/str=%d  (base_m=%d <= J/str OK: %s)\n",
            task, Psi, pop$truth$se_mc, Hstr, pop$params$J_per_stratum,
            base_m, base_m <= pop$params$J_per_stratum))
stopifnot(base_m <= pop$params$J_per_stratum)   # stage-1 SRS-WOR feasibility

# ---- parallel cluster (mirror run_sim.R worker setup) ------------------------
cl <- makeCluster(cores)
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "R")
  source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
  source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
})
clusterExport(cl, c("one_rep_m"), envir = environment())
# per-task RNG stream; offset by m_idx so different m / chunks never collide
clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 1000L * job$m_idx + task)

# ---- run all rungs for this (m, chunk) ---------------------------------------
all_res <- list(); all_diag <- list()
for (rl in names(RUNGS)) {
  learners <- RUNGS[[rl]]
  t0 <- Sys.time()
  reps_out <- parLapply(cl, reps, one_rep_m, pop = pop, learners = learners, base_m = base_m)
  res  <- do.call(rbind, lapply(reps_out, `[[`, "results"))
  diag <- do.call(rbind, lapply(reps_out, `[[`, "diag"))
  # tag with the swept knobs + rung label so aggregate.R / plots have them
  res$scenario   <- SCENARIO
  res$rung_label <- rl
  res$base_m     <- base_m
  res$m_total    <- Hstr * base_m
  res$Psi        <- Psi
  diag$rung_label <- rl; diag$base_m <- base_m
  all_res[[rl]]  <- res
  all_diag[[rl]] <- diag
  cat(sprintf("[task %d] rung=%-14s done (%d reps, %.1f min)\n",
              task, rl, length(reps), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res_all  <- do.call(rbind, all_res);  rownames(res_all)  <- NULL
diag_all <- do.call(rbind, all_diag); rownames(diag_all) <- NULL

# ---- per-task output (full detail; NON-DESTRUCTIVE path) ---------------------
out <- list(run_id = RUN_ID, scenario = SCENARIO, base_m = base_m,
            m_total = Hstr * base_m, m_idx = job$m_idx, chunk = chunk,
            reps = reps, rungs = RUNGS, Psi = Psi, truth = pop$truth,
            params = pop$params, results = res_all, diagnostics = diag_all,
            smoke = SMOKE)
saveRDS(out, fn)   # `tag`/`fn` were computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s  (%d est rows, %d diag rows)\n",
            task, fn, nrow(res_all), nrow(diag_all)))

# ---- reproducibility manifest (one per task; race-free) ----------------------
manifest <- list(
  run_id = RUN_ID, scenario = SCENARIO, base_m = base_m, m_total = Hstr * base_m,
  m_idx = job$m_idx, chunk = chunk, reps = reps, rungs = RUNGS,
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
mf <- file.path(MAN_OUT, sprintf("%smanifest_m%02d_chunk%03d.rds", tag, base_m, chunk))
saveRDS(manifest, mf)
cat(sprintf("[task %d] manifest %s\n", task, mf))
