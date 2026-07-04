# =====================================================================
# run.R  —  R19_rate_sweep driver  (NON-DESTRUCTIVE ARC enhancement run)
#
# GOAL (Writing/comments/phase4-arc-sim-specs.md, "√m-scaled product-error
# diagnostic across the F.5 m-sweep" — the re-review's sharpest evidence ask):
# show the realized cross-fit nuisance PRODUCT error scaled by sqrt(m) ACROSS
# the large-m sweep, to probe whether the assumed product rate (C1)
#   ||Qhat - Q0|| * ||ghat - g0||  =  o_P( m^{-1/2} )
# is plausibly met at each ladder rung. A single-m diagnostic (R04) cannot test
# a RATE; R02 swept base_m but recorded no nuisance errors. R19 CROSSES them:
# R04's truth-extraction/L2 machinery x R02's base_m knob.
#   Expected: the sqrt(m_total)-scaled product is FLAT/DECLINING across the
#   sweep at L2-L3 (rate plausibly met), GROWING at L4 (interpolating deep RF
#   -- rate fails), and PLAUSIBLY GROWING at L1 (SL.glm is misspecified for
#   both nuisances under the complex covariates -> positive approximation
#   floors; matches R02's locked L1 bias floor + coverage decay -- a finding,
#   NOT a bug). STOP-and-report outcomes are only the genuinely contradictory
#   ones: L4 flat/declining, or L2 clearly growing (see NOTES.md).
#
# This driver SOURCES the canonical engine read-only, sources R04's OWN helpers
# (codes/arc_runs/R04_nuisance_rate/nuisance_rate_helpers.R — reused, NOT
# copied), and adds only the thin base_m-forwarding wrapper in
# rate_sweep_helpers.R. It NEVER edits codes/*.R or any other arc_run, and
# writes only into sim_output/arc_runs/R19_rate_sweep/ + results/arc/.
#
# CELL GRID: scenario 'standard' (FIXED) x base_m {6,12,20,30} x rung L1-L4
#   = 16 cells; FULL = 200 reps, SIM_CHUNK=50 -> 4 chunks/cell
#   -> FULL --array = 1-64  (16 cells x 4 chunks; keep submit.slurm in sync).
# Grid order (cells row-major, chunk fastest): expand.grid puts base_m fastest
# within rung, so tasks 1-16 = L1_param (m6 c1-4, m12 c1-4, m20 c1-4, m30 c1-4),
# 17-32 = L2_smooth, 33-48 = L3_adaptive, 49-64 = L4_aggressive (heaviest).
#
# m CONVENTION (matches R04): m = total #sampled PSUs, m_total = H*base_m
#   = {60, 120, 200, 300} for standard (H=10). base_m AND m_total recorded per
#   row so results join the existing R02/R04 tables.
#
# Env vars (with local fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index             (default 1)
#   SLURM_CPUS_PER_TASK  cores                          (default detectCores-1)
#   SIM_N_REPS           total reps per cell            (default 200)
#   SIM_CHUNK            reps per array task            (default 50)
#   SIM_CODE             canonical code dir             (default "codes")
#   R19_DIR              this run's folder              (default REPO_ROOT/codes/arc_runs/R19_rate_sweep)
#   R19_OUT              per-task RDS out dir           (default DATA_ROOT/arc_runs/R19_rate_sweep)
#   R04_DIR              R04's folder (for its helpers) (default SIM_CODE/arc_runs/R04_nuisance_rate)
#   SMOKE                "1" -> tiny 1-cell, 10-rep serial local validation
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))     # exposes .sl(), make_cf_folds(), etc.
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))       # SL_LADDER, custom RF wrappers
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# ---- R04's truth-extraction + L2 instrumentation, sourced from R04's OWN
#      folder (REUSED read-only, NOT copied; provides attach_truth,
#      .g0_marginal, cf_nuisance_oof, .wL2, one_rep_rate) -------------------
.r04_dir <- Sys.getenv("R04_DIR", file.path(CODE, "arc_runs", "R04_nuisance_rate"))
source(file.path(.r04_dir, "nuisance_rate_helpers.R"))

# ---- run-local wrapper (one_rep_rate_m forwards base_m); resolve our folder --
.this_dir <- Sys.getenv("R19_DIR",
  file.path(REPO_ROOT, "simulation", "enhancements", "R19_rate_sweep"))
source(file.path(.this_dir, "rate_sweep_helpers.R"))

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 200L)
CHUNK  <- geti("SIM_CHUNK", 50L)

# ---- output dirs (NEVER clobber locked results) ----
OUT     <- Sys.getenv("R19_OUT", file.path(DATA_ROOT, "arc_runs", "R19_rate_sweep"))
MAN_OUT <- file.path(OUT, "manifest")
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
for (d in c(OUT, MAN_OUT, ARC_RESULTS))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: base_m x rung at FIXED scenario 'standard' -------------------
SCENARIO <- "standard"            # FIXED per spec (the F.5 m-sweep scenario)
M_LEVELS <- c(6L, 12L, 20L, 30L)  # = R02's sweep; m_total = 10*base_m = 60..300
RUNGS    <- names(SL_LADDER)      # L1_param, L2_smooth, L3_adaptive, L4_aggressive
cells <- expand.grid(base_m = M_LEVELS, rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

TRUTH_M <- 2e6L
if (SMOKE) {
  # tiny validation: ONE cheap cell (standard x L1_param x base_m=30), 10 reps,
  # serial, cheap truth integral. ~2-4 min locally.
  cells   <- cells[cells$rung == "L1_param" & cells$base_m == 30L, , drop = FALSE]
  N_REPS  <- 10L; CHUNK <- 10L; cores <- 1L; task <- 1L
  TRUTH_M <- 2e5L
  cat("[SMOKE] standard x L1_param x base_m=30, 10 reps, serial, truth_M=2e5\n")
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci)
  cbind(cells[ci, , drop = FALSE], cell_index = ci, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))   # FULL grid = 16*4 = 64 -> --array=1-64
job        <- grid[task, ]
base_m     <- job$base_m
rung       <- job$rung
chunk      <- job$chunk
cell_index <- job$cell_index
learners   <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi          # GLOBAL rep indices -> sample seeds match R02/R04

cat(sprintf("[task %d] scenario=%s base_m=%d rung=%s chunk=%d reps=%d-%d cores=%d learners=%s\n",
            task, SCENARIO, base_m, rung, chunk, rep_lo, rep_hi, cores,
            paste(learners, collapse = "+")))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr19_m%02d_%s_chunk%03d.rds", tag, base_m, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ----
# Population geometry is FIXED; only the SAMPLING base_m varies, so the SAME
# population serves every cell (and is identical to the locked headline run's).
pop <- make_population(SCENARIO, model_type = "complex", truth_M = TRUTH_M)
stopifnot(base_m <= pop$params$J_per_stratum)   # stage-1 SRS-WOR feasibility (60 for standard)
m_total <- pop$params$H * base_m
cat(sprintf("[task %d] Psi=%.5f Psi_N=%.5f gap=%.6f  H=%d -> m_total=%d\n",
            task, pop$truth$psi, pop$truth$psi_N_Q0,
            abs(pop$truth$gap_super_census), pop$params$H, m_total))

# ---- run the chunk's reps (parallel over reps, mirroring R04) ----------------
run_one <- function(i) tryCatch(
  one_rep_rate_m(i, pop, learners, rung, SCENARIO, base_m = base_m,
                 V_cf = 5L, g_oof_bound = 0.05),
  error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    # R04 helpers (attach_truth etc.) + R19 wrapper, same resolution as master
    .r04_dir <- Sys.getenv("R04_DIR", file.path(CODE, "arc_runs", "R04_nuisance_rate"))
    source(file.path(.r04_dir, "nuisance_rate_helpers.R"))
    .this_dir <- Sys.getenv("R19_DIR",
      file.path(REPO_ROOT, "simulation", "enhancements", "R19_rate_sweep"))
    source(file.path(.this_dir, "rate_sweep_helpers.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  # RNG stream offset: 1000L*cell_index + task. cell_index in 1..16 and
  # task in 1..64 < 1000, so no two (cell, chunk) tasks ever share a stream.
  # (Worker RNG only drives fold splits/SL CV; sample draws use scoped seeds.)
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 1000L * cell_index + task)
  # export EVERYTHING run_one's closure needs on the workers (R04 lesson: pop!)
  clusterExport(cl, c("pop", "learners", "rung", "SCENARIO", "base_m"),
                envir = environment())
  rows <- parLapply(cl, reps, run_one)
} else {
  set.seed(SAMPLE_SEED_BASE + 1000L * cell_index + task)   # same offset, serial
  rows <- lapply(reps, run_one)
}
rows <- Filter(Negate(is.null), rows)   # drop reps that errored
if (!length(rows)) stop(sprintf("[task %d] all %d reps in this chunk failed", task, length(reps)))
res <- do.call(rbind, rows)

# ---- per-task RDS (full detail; never clobbers locked outputs) ----
out <- list(run = "R19_rate_sweep", scenario = SCENARIO, rung = rung,
            base_m = base_m, m_total = m_total,
            learners = learners, chunk = chunk, reps = reps,
            Psi = pop$truth$psi, truth = pop$truth, params = pop$params,
            per_rep = res, smoke = SMOKE)
saveRDS(out, fn)   # `tag`/`fn` were computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s  (%d rep rows)\n", task, fn, nrow(res)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----
manifest <- list(
  run = "R19_rate_sweep",
  scenario = SCENARIO, rung = rung, base_m = base_m, m_total = m_total,
  learners = learners, chunk = chunk, reps = reps, cell_index = cell_index,
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED,
               rng_stream_iseed = SAMPLE_SEED_BASE + 1000L * cell_index + task),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  truth_join = unique(res$truth_join))
mf <- file.path(MAN_OUT, sprintf("%smanifest_m%02d_%s_chunk%03d.rds", tag, base_m, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[task %d] manifest %s  (truth_join=%s)\n",
            task, mf, paste(unique(res$truth_join), collapse = ",")))

# ---- SMOKE: aggregate inline + explicit gate verdict -------------------------
if (SMOKE) {
  cat("\n[SMOKE] single-cell rate summary (standard x L1_param x base_m=30):\n")
  agg <- data.frame(
    scenario = SCENARIO, rung = rung, base_m = base_m, m_total = m_total,
    n_reps = nrow(res), mean_n = mean(res$n), mean_m_psu = mean(res$m_psu),
    mean_eQ = mean(res$eQ_int), mean_eg = mean(res$eg_int),
    mean_prod = mean(res$prod_int),
    mean_prod_sqrtm = mean(res$prod_int_sqrtm),
    truth_join = paste(unique(res$truth_join), collapse = ","))
  print(format(agg, digits = 4), row.names = FALSE)

  tj_ok  <- all(res$truth_join == "OK")
  fin_ok <- all(is.finite(res$eQ_int)) && all(is.finite(res$eg_int)) &&
            all(is.finite(res$prod_int)) && all(is.finite(res$prod_int_sqrtm)) &&
            all(is.finite(res$eQ_real)) && all(is.finite(res$eg_real)) &&
            all(is.finite(res$prod_real)) && all(is.finite(res$prod_real_sqrtm))
  l1_ok  <- is.finite(mean(res$prod_int_sqrtm)) && mean(res$prod_int_sqrtm) < 1
  if (tj_ok && fin_ok && l1_ok) {
    cat(sprintf("[SMOKE-GATE] PASS: truth_join=OK, all error norms finite, L1 mean_prod_sqrtm=%.4f < 1\n",
                mean(res$prod_int_sqrtm)))
  } else {
    cat(sprintf("[SMOKE-GATE] STOP: truth_join_ok=%s finite_ok=%s l1_prod_sqrtm_lt1=%s (mean_prod_sqrtm=%.4f) -- investigate before submitting\n",
                tj_ok, fin_ok, l1_ok, mean(res$prod_int_sqrtm)))
  }
}
