# =====================================================================
# run.R  —  R18_cvu driver  (NON-DESTRUCTIVE ARC enhancement run)
#
# SPEC A4 (Writing/comments/phase4-arc-sim-specs.md): unweighted-nuisance
# internal-CV simulation arm (CV-u). The application's CV arm is UNWEIGHTED
# (real-data path of run_estimators) but the simulation's CV foil is WEIGHTED
# (CV-w); this run adds the unweighted variant to the sim so the
# "internal CV != cross-fitting" contrast is shown for the object the
# application actually reports. ONE arm per rep ("Fully-Aware-CVu"; cheap).
#
# Cells = 2 scenarios (standard, R1) x 2 rungs (L2_smooth, L3_adaptive) = 4.
# FULL = SIM_N_REPS=1000, SIM_CHUNK=100 -> 4 cells x 10 chunks = 40 tasks
#        => sbatch --array=1-40   (grid is CELL-MAJOR; task order:
#        1-10  standard x L2_smooth   (chunks 1-10)
#        11-20 R1       x L2_smooth
#        21-30 standard x L3_adaptive
#        31-40 R1       x L3_adaptive)
#
# COMPARABILITY (key property — do not break): each rep i draws its sample via
# draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i) with i the GLOBAL rep
# index, exactly like run_sim.R — so the samples are BYTE-IDENTICAL to the
# locked headline run for the same scenario, and the CVu rows can be laid next
# to the locked Fully-Aware-CV / Fully-Aware-CF / Fully-Aware rows.
#
# This driver SOURCES the canonical engine read-only and adds ONLY cvu_arm.R.
# It NEVER edits codes/*.R and writes only into R18_OUT (+ results/arc via
# aggregate.R afterwards).
#
# Env vars (with local fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores                          (default detectCores-1)
#   SIM_N_REPS           total reps per cell            (default 1000)
#   SIM_CHUNK            reps per array task            (default 100)
#   SIM_CODE             canonical code dir             (default "codes")
#   R18_DIR              this run's folder              (default REPO_ROOT/codes/arc_runs/R18_cvu)
#   R18_OUT              per-task RDS out dir           (default DATA_ROOT/arc_runs/R18_cvu)
#   SMOKE                "1" -> tiny 1-cell local validation (standard x L2, 10 reps)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))     # exposes .eif_from_tmle, .se_des (reused read-only)
source(file.path(CODE, "diagnostics.R"))    # deff_clust, audit_population
source(file.path(CODE, "learners.R"))       # SL_LADDER
# run-local helper (the CV-u arm); resolve this folder via R18_DIR
.this_dir <- Sys.getenv("R18_DIR",
  file.path(REPO_ROOT, "simulation", "enhancements", "R18_cvu"))
source(file.path(.this_dir, "cvu_arm.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)

# ---- output dirs (NEVER clobber locked results) ----
OUT <- Sys.getenv("R18_OUT", file.path(DATA_ROOT, "arc_runs", "R18_cvu"))
MAN <- file.path(OUT, "manifest")
for (.d in c(OUT, MAN)) if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x rung x chunk (CV arm exists only on L2/L3) ----
SCEN  <- c("standard", "R1")
RUNGS <- c("L2_smooth", "L3_adaptive")

if (SMOKE) {
  # tiny validation: ONE cheap cell (standard x L2_smooth), 10 reps, 2 cores.
  SCEN <- "standard"; RUNGS <- "L2_smooth"
  N_REPS <- 10L; CHUNK <- 10L; cores <- min(cores, 2L); task <- 1L
  cat("[SMOKE] standard x L2_smooth, 10 reps, <=2 cores\n")
}

cells <- expand.grid(scenario = SCEN, rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))     # FULL grid = 4 x 10 = 40 (--array=1-40)
job <- grid[task, ]
scenario <- job$scenario; rung <- job$rung
chunk <- job$chunk; cell_index <- job$cell_index
learners <- SL_LADDER[[rung]]
# GLOBAL rep index so chunks never collide and sample seeds match run_sim.R exactly.
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[R18 task %d] scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d learners=%s%s\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, cores,
            paste(learners, collapse = "+"), if (SMOKE) " [SMOKE]" else ""))

# ---- per-chunk checkpoint: compute fn EARLY, skip instantly if already done ----
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr18_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ----
pop <- make_population(scenario, model_type = "complex",
                       truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[R18 task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0, abs(pop$truth$gap_super_census)))

# ---- one rep: sample (seed-identical to the locked run) + the CV-u arm ONLY ----
one_rep <- function(i, pop, learners) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_cvu(obs, learners = learners)     # inner_cv_folds=5, g_oof_bound=0.05 (engine defaults)
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_cvu, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow, deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    .this_dir <- Sys.getenv("R18_DIR",
      file.path(REPO_ROOT, "simulation", "enhancements", "R18_cvu"))
    source(file.path(.this_dir, "cvu_arm.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  # RNG-stream offset: iseed = SAMPLE_SEED_BASE + 1000L*cell_index + task.
  # task (1..40) is already unique across the grid, so streams are unique WITHIN
  # this run; the 1000L*cell_index term keeps them disjoint from the OLDER arc
  # runs' SAMPLE_SEED_BASE + task range (R01/R03/R04/R07/R11/R12). The sibling
  # batch runs (R13/R14/R15/R16/R19, and R02 via 1000L*m_idx + task) share this
  # same formula, so cross-run iseed coincidences DO occur -- harmless, because
  # the worker stream drives only learner-internal RNG while samples self-seed
  # via draw_sample(SAMPLE_SEED_BASE + global rep index) (see NOTE below).
  # NOTE: the SAMPLES are unaffected — draw_sample() seeds itself with
  # SAMPLE_SEED_BASE + i (global rep index) via .scoped_seed, preserving
  # byte-identity with the locked headline run; the worker stream governs only
  # learner-internal RNG (SuperLearner CV splits, folds.svy, ranger).
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 1000L * cell_index + task)
  # workers need NOTHING beyond the sourced files: `pop` and `learners` are
  # passed as explicit parLapply arguments (the R04 'missing pop' bug class);
  # export them anyway as belt-and-braces against future closure edits.
  clusterExport(cl, c("pop", "learners"), envir = environment())
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners)
} else {
  set.seed(SAMPLE_SEED_BASE + 1000L * cell_index + task)   # same offset, serial path
  reps_out <- lapply(reps, one_rep, pop = pop, learners = learners)
}
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored
if (!length(reps_out)) stop(sprintf("[R18 task %d] all %d reps in this chunk failed", task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task RDS (full detail; never clobbers locked outputs) ----
out <- list(run = "R18_cvu", scenario = scenario, rung = rung,
            learners = learners, chunk = chunk, reps = reps,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags, smoke = SMOKE)
saveRDS(out, fn)   # `tag`/`fn` were computed early for the checkpoint guard
cat(sprintf("[R18 task %d] saved %s  (%d est rows, %d diag rows)\n",
            task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task, in manifest/ subdir) ----
manifest <- list(
  run = "R18_cvu", scenario = scenario, rung = rung, learners = learners,
  chunk = chunk, reps = reps, cell_index = cell_index,
  params = pop$params, truth = pop$truth, pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED,
               rng_stream_iseed = SAMPLE_SEED_BASE + 1000L * cell_index + task),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]])
mf <- file.path(MAN, sprintf("%smanifest_r18_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[R18 task %d] manifest %s\n", task, mf))

# ---- SMOKE: aggregate inline + print the explicit smoke-gate verdict ----
if (SMOKE) {
  crit <- qt(0.975, pmax(1, res$df))
  cov  <- mean(abs(res$b - Psi) <= crit * res$se)
  agg  <- data.frame(scenario = scenario, rung = rung, method = res$method[1],
                     n_reps = nrow(res), Psi = Psi,
                     bias = mean(res$b) - Psi, emp_sd = sd(res$b),
                     mean_se = mean(res$se), se_ratio = mean(res$se) / sd(res$b),
                     coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(res)),
                     mean_df = mean(res$df), df_design = res$df_design[1],
                     g_cvu_min = mean(diags$g_cvu_min), eps_cvu = mean(diags$eps_cvu))
  cat("\n==== R18 SMOKE summary (", scenario, " x ", rung, ", ", nrow(res), " reps) ====\n", sep = "")
  num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
  print(ap, row.names = FALSE)

  # gate: finite b/se, df ~ design df, coverage at 10 reps in [0.6, 1]
  finite_ok <- all(is.finite(res$b)) && all(is.finite(res$se))
  df_ok     <- all(abs(res$df - res$df_design) <= 2)
  cov_ok    <- is.finite(cov) && cov >= 0.6 && cov <= 1
  if (finite_ok && df_ok && cov_ok) {
    cat(sprintf("\n[SMOKE-GATE] PASS: all b/se finite; df (mean %.1f) matches design df (%d); coverage %.2f in [0.6,1]. Submit the full --array=1-40 run.\n",
                mean(res$df), res$df_design[1], cov))
  } else {
    cat(sprintf("\n[SMOKE-GATE] STOP: finite_ok=%s df_ok=%s (mean df %.1f vs design %d) cov_ok=%s (coverage=%.2f). Investigate before submitting.\n",
                finite_ok, df_ok, mean(res$df), res$df_design[1], cov_ok, cov))
  }
}
