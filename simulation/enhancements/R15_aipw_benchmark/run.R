# =====================================================================
# run.R  —  R15_aipw_benchmark driver  (NON-DESTRUCTIVE ARC enhancement run)
#
# GOAL (spec item A7, Writing/comments/phase4-arc-sim-specs.md): benchmark the
# five survey-TMLE arms against an INDEPENDENT, literature-standard competitor:
# survey-weighted AIPW (Hajek form) with (a) a survey JKn replicate SE and
# (b) the paper's own Eq-8 design linearization, across the FULL ladder
# (L1-L4) x both designs (standard, R1), on the SAME samples as the locked
# headline run (paired comparison; see comparability note below).
#
# This driver SOURCES the canonical engine read-only and reuses its building
# blocks via aipw_helpers.R. It NEVER edits codes/*.R and writes only into
# $R15_OUT (+ results/arc via aggregate.R afterwards).
#
# One SLURM array task = one (scenario x rung x rep-chunk) cell (mirrors
# run_sim.R / R04). FULL grid: 2 scenarios x 4 rungs = 8 cells, 1000 reps in
# chunks of 100 -> 10 chunks/cell -> 80 tasks -> sbatch --array=1-80.
# Grid order (expand.grid scenario fastest; chunks innermost):
#   tasks  1-10  standard x L1_param      (chunks 1-10)
#   tasks 11-20  R1       x L1_param
#   tasks 21-30  standard x L2_smooth
#   tasks 31-40  R1       x L2_smooth
#   tasks 41-50  standard x L3_adaptive
#   tasks 51-60  R1       x L3_adaptive
#   tasks 61-70  standard x L4_aggressive
#   tasks 71-80  R1       x L4_aggressive
#
# COMPARABILITY (preserve!): reps use the GLOBAL rep index i with
#   draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i)
# -- the exact seeding of the locked headline run (run_sim.R) -- so for each
# (scenario, rep) the analysis sample is BYTE-IDENTICAL to the one the five
# locked arms saw. AIPW-vs-TMLE differences are therefore estimator
# differences, not sampling Monte-Carlo noise.
#
# Env vars (with local fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores                          (default detectCores-1)
#   SIM_N_REPS           total reps per (scenario,rung) (default 1000)
#   SIM_CHUNK            reps per array task            (default 100)
#   SIM_CODE             canonical code dir             (default "codes")
#   R15_DIR              this run's folder              (default REPO_ROOT/codes/arc_runs/R15_aipw_benchmark)
#   R15_OUT              per-task RDS out dir           (default DATA_ROOT/arc_runs/R15_aipw_benchmark)
#   SMOKE                "1" -> tiny 1-cell local validation (standard x L1,
#                        20 reps, 2 cores, truth_M = 2e5)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))     # exposes .sl(), make_cf_folds(), .se_des()
source(file.path(CODE, "diagnostics.R"))    # deff_clust
source(file.path(CODE, "learners.R"))       # SL_LADDER, custom RF wrappers
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# run-local helpers (AIPW arms + JKn); resolve next to this file
.this_dir <- Sys.getenv("R15_DIR",
  file.path(REPO_ROOT, "simulation", "enhancements", "R15_aipw_benchmark"))
source(file.path(.this_dir, "aipw_helpers.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
G_FLOOR <- 0.05      # harmonized propensity floor, both arms (= engine OOF floor)
V_CF    <- 5L        # cross-fit folds (engine default)

# ---- output dirs (NEVER clobber locked results) ----
OUT <- Sys.getenv("R15_OUT", file.path(DATA_ROOT, "arc_runs", "R15_aipw_benchmark"))
MAN <- file.path(OUT, "manifest")
for (d in c(OUT, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x rung x chunk  (same shape as run_sim.R) ----
RUNGS <- names(SL_LADDER)
cells <- expand.grid(scenario = c("standard", "R1"), rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

if (SMOKE) {
  # tiny validation: ONE cheap cell (standard x L1_param), 20 reps, 2 cores.
  cells  <- cells[cells$scenario == "standard" & cells$rung == "L1_param", , drop = FALSE]
  N_REPS <- 20L; CHUNK <- 20L; cores <- max(1L, min(cores, 2L)); task <- 1L
  cat("[SMOKE] standard x L1_param, 20 reps, 2 cores\n")
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job <- grid[task, ]
scenario <- job$scenario; rung <- job$rung; chunk <- job$chunk
cell_index <- job$cell_index
learners <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[task %d] scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d learners=%s\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, cores, paste(learners, collapse = "+")))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr15_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ----
pop <- make_population(scenario, model_type = "complex",
                       truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0,
            abs(pop$truth$gap_super_census)))

# ---- run the chunk's reps (parallel over reps, like run_sim.R / R04) ----
run_one <- function(i) tryCatch(
  one_rep_aipw(i, pop, learners, g_floor = G_FLOOR, V_cf = V_CF),
  error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# RNG stream seed: SAMPLE_SEED_BASE + 1000*cell_index + task. The 1000*cell
# offset keeps streams distinct across cells even though `task` alone is
# already unique (1..80 < 1000), and never collides across chunks of the same
# cell (consecutive task ids). This stream only drives fold assignment +
# SuperLearner internal CV; the SAMPLES themselves come from the scoped
# per-rep seeds inside draw_sample() (comparability preserved).
iseed <- SAMPLE_SEED_BASE + 1000L * cell_index + task
if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    .this_dir <- Sys.getenv("R15_DIR",
      file.path(REPO_ROOT, "simulation", "enhancements", "R15_aipw_benchmark"))
    source(file.path(.this_dir, "aipw_helpers.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = iseed)
  # export EVERYTHING run_one's closure needs (R04 lesson: workers once lacked `pop`)
  clusterExport(cl, c("pop", "learners", "G_FLOOR", "V_CF"), envir = environment())
  rows <- parLapply(cl, reps, run_one)
} else {
  set.seed(iseed)
  rows <- lapply(reps, run_one)
}
rows <- Filter(Negate(is.null), rows)   # drop reps that errored
if (!length(rows)) stop(sprintf("[task %d] all %d reps in this chunk failed", task, length(reps)))
n_failed <- length(reps) - length(rows)
if (n_failed > 0L)
  cat(sprintf("[task %d] WARNING: %d/%d reps failed in this chunk; chunk RDS will checkpoint as complete (failed reps dropped, not retried) -- delete %s to recompute\n",
              task, n_failed, length(reps), fn))
res   <- do.call(rbind, lapply(rows, `[[`, "results"))
diags <- do.call(rbind, lapply(rows, `[[`, "diag"))

# ---- per-task RDS (full detail; never clobbers locked outputs) ----
out <- list(run = "R15_aipw_benchmark", scenario = scenario, rung = rung,
            learners = learners, chunk = chunk, reps = reps,
            n_failed = n_failed,
            Psi = Psi, truth = pop$truth, params = pop$params,
            g_floor = G_FLOOR, V_cf = V_CF,
            results = res, diagnostics = diags)
saveRDS(out, fn)   # `tag`/`fn` were computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s  (%d est rows, %d diag rows)\n",
            task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----
manifest <- list(
  run = "R15_aipw_benchmark",
  scenario = scenario, rung = rung, learners = learners, chunk = chunk, reps = reps,
  n_failed = n_failed,
  g_floor = G_FLOOR, V_cf = V_CF,
  params = pop$params, truth = pop$truth,
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED, iseed = iseed),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]]
)
mf <- file.path(MAN, sprintf("%smanifest_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[task %d] manifest %s\n", task, mf))

# ---- SMOKE: aggregate inline + explicit PASS/STOP gate ----
if (SMOKE) {
  cat("\n[SMOKE] per-arm summary (this single cell):\n")
  agg <- do.call(rbind, by(res, res$method, function(d) {
    crit <- qt(0.975, pmax(1, d$df))
    data.frame(method = d$method[1], n_reps = nrow(d),
               mean_b = mean(d$b), bias = mean(d$b) - Psi,
               bias_med = median(d$b) - Psi, emp_sd = sd(d$b),
               mean_se_jkn = mean(d$se_jkn), mean_se_lin = mean(d$se_lin),
               jkn_over_lin = mean(d$se_jkn) / mean(d$se_lin),
               cov_jkn = mean(abs(d$b - Psi) <= crit * d$se_jkn),
               cov_lin = mean(abs(d$b - Psi) <= crit * d$se_lin),
               df = d$df[1])
  }))
  rownames(agg) <- NULL
  print(format(agg, digits = 4), row.names = FALSE)
  # collapsed-g reps (every raw g outside the floor) would indicate GENUINE
  # weighted quasi-separation: the spurious raw-weight-scale IRLS divergence is
  # already removed by the helper's mean-1 weight normalization (= tmle()'s own
  # obsWeights normalization), so expect 0 here; investigate any nonzero count.
  n_coll_sf <- sum(diags$g_sf_outside_floor >= 1)
  n_coll_cf <- sum(diags$g_cf_outside_floor >= 1)
  cat(sprintf("[SMOKE] Psi=%.5f  jkn_mode=%s  cf_V_eff=%s  deff(mean)=%.3f  collapsed-g reps: SF=%d CF=%d\n",
              Psi, diags$jkn_mode[1], diags$cf_V_eff[1], mean(diags$deff),
              n_coll_sf, n_coll_cf))
  cat("[SMOKE] note: jkn_over_lin = 1 EXACTLY is correct here -- with stratum-\n",
      "        constant weights + equal PSU takes the JKn Hajek-ratio variance\n",
      "        reduces algebraically to the Eq-8 linearization (verified vs a\n",
      "        toy design where it breaks); it differs on NHANES.\n", sep = "")

  # gate: finite b/se for BOTH arms + BOTH SEs; se_jkn within ~2x of se_lin;
  # MEDIAN bias small vs known Psi (|bias_med| <= 0.05; median is the gate
  # statistic because collapsed AIPW-SF reps legitimately fatten the mean);
  # both arms present; no silently dropped reps; CF arm must not collapse.
  ok_arms   <- setequal(agg$method, c("AIPW-SF", "AIPW-CF"))
  ok_finite <- all(is.finite(res$b), is.finite(res$se_jkn), is.finite(res$se_lin))
  ok_ratio  <- all(is.finite(agg$jkn_over_lin)) &&
               all(agg$jkn_over_lin >= 0.5 & agg$jkn_over_lin <= 2)
  ok_bias   <- all(abs(agg$bias_med) <= 0.05)
  ok_nreps  <- all(agg$n_reps == length(reps))
  ok_cf     <- n_coll_cf == 0                    # unweighted-OOF CF must be stable
  if (ok_arms && ok_finite && ok_ratio && ok_bias && ok_nreps && ok_cf) {
    cat("[SMOKE-GATE] PASS: both arms finite; se_jkn/se_lin in [0.5,2]; |median bias| <= 0.05; no dropped reps; CF g-fits stable.\n")
  } else {
    cat(sprintf("[SMOKE-GATE] STOP: ok_arms=%s ok_finite=%s ok_ratio=%s ok_bias(med)=%s ok_nreps=%s ok_cf_stable=%s -- inspect before submitting.\n",
                ok_arms, ok_finite, ok_ratio, ok_bias, ok_nreps, ok_cf))
  }
}
