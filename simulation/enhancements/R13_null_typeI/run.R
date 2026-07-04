# =====================================================================
# run.R  —  R13_null_typeI driver  (NON-DESTRUCTIVE ARC enhancement run)
#
# SPEC: Writing/comments/phase4-arc-sim-specs.md item A5 — "Null / type-I
# error (psi=0) + small-effect power".
#
# GOAL: confirm type-I control at a TRUE NULL (te_log_odds = 0; the E3
# "apparent association dissolves" story is never simulated at psi = 0)
# and add a small-effect power point (te_log_odds = 0.3). Both designs
# (standard, R1), the full ladder L1-L4, all five arms via
# run_estimators() UNCHANGED. aggregate.R (this folder) computes the
# rejection rate of H0: psi = 0 (target 0.05 at the null; the
# Fully-Aware-CF rows are the audit target) next to the standard
# coverage-of-true-Psi columns.
#
# NON-DESTRUCTIVE: this driver SOURCES the canonical engine
# (config/dgp/estimators/diagnostics/learners) read-only. The ONLY new
# knob is te_log_odds, which is a documented make_population() argument.
# Outputs go to a private sim_output/arc_runs/R13_null_typeI/ folder
# (+ results/arc/ via aggregate.R); locked results are never touched.
#
# GRID: te in {0, 0.3} x scenario in {standard, R1} x rung L1-L4 = 16 cells,
# each x ceil(SIM_N_REPS/SIM_CHUNK) rep-chunks.
#   FULL: SIM_N_REPS=1000, SIM_CHUNK=100 -> 16 x 10 = 160 tasks -> --array=1-160
# Cells are te-major (all te=0 cells first), within te expand.grid(scenario,
# rung) with scenario fastest; chunks are fastest within a cell:
#   tasks   1- 10  te=0.0 standard L1_param   |  tasks  81- 90  te=0.3 standard L1_param
#   tasks  11- 20  te=0.0 R1       L1_param   |  tasks  91-100  te=0.3 R1       L1_param
#   tasks  21- 30  te=0.0 standard L2_smooth  |  ...
#   tasks  71- 80  te=0.0 R1   L4_aggressive  |  tasks 151-160  te=0.3 R1   L4_aggressive
#
# Env vars (local fallbacks in brackets):
#   SLURM_ARRAY_TASK_ID  1-based task index                  [1]
#   SLURM_CPUS_PER_TASK  cores for this task                 [detectCores-1]
#   SIM_N_REPS           total reps per cell                 [1000]
#   SIM_CHUNK            reps per array task                 [100]
#   SIM_CODE             canonical code dir                  [codes]
#   R13_DIR              this run folder    [REPO_ROOT/codes/arc_runs/R13_null_typeI]
#   R13_OUT              per-task RDS dir   [DATA_ROOT/arc_runs/R13_null_typeI]
#   SMOKE                "1" -> tiny 1-cell validation (te=0 x standard x L1, 20 reps)
# =====================================================================

RUN_ID <- "R13_null_typeI"

# ---- locate + source the canonical engine (read-only) ----------------
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))      # POP_SEED, SAMPLE_SEED_BASE, paths, RNGkind
source(file.path(CODE, "dgp.R"))         # make_population (te_log_odds knob), draw_sample
source(file.path(CODE, "estimators.R"))  # run_estimators (UNCHANGED, all five arms)
source(file.path(CODE, "diagnostics.R")) # deff_clust, audit_population
source(file.path(CODE, "learners.R"))    # SL_LADDER
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# run folder (no run-local helpers needed; resolved + recorded for the manifest)
.this_dir <- Sys.getenv("R13_DIR",
  file.path(REPO_ROOT, "simulation", "enhancements", RUN_ID))

# ---- env helpers ------------------------------------------------------
geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)

# ---- output dirs (private; NEVER the locked sim_output/intermediate) --
OUT <- Sys.getenv("R13_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
MAN_DIR <- file.path(OUT, "manifest")
for (d in c(OUT, MAN_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: te x scenario x rung (te-major; scenario fastest) -----
RUNGS <- names(SL_LADDER)                              # all 4 rungs L1-L4
TE_GRID <- c(0, 0.3)                                   # null + small-effect power point
cells <- do.call(rbind, lapply(TE_GRID, function(te)
  cbind(expand.grid(scenario = c("standard", "R1"), rung = RUNGS,
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE),
        te = te)))

if (SMOKE) {
  # tiny validation: ONE cheap cell (te=0 x standard x L1_param), 20 reps, 2 cores.
  cells  <- cells[cells$te == 0 & cells$scenario == "standard" &
                  cells$rung == "L1_param", , drop = FALSE]
  N_REPS <- 20L; CHUNK <- 20L
  cores  <- min(cores, 2L)
  task   <- 1L
  cat(sprintf("[SMOKE] %s: te=0 x standard x L1_param, %d reps, %d cores\n",
              RUN_ID, N_REPS, cores))
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))      # FULL grid = 16 cells x 10 chunks = 160 tasks
job        <- grid[task, ]
scenario   <- job$scenario
rung       <- job$rung
te         <- job$te
cell_index <- job$cell_index
chunk      <- job$chunk
learners   <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[%s task %d] te=%g scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d learners=%s%s\n",
            RUN_ID, task, te, scenario, rung, chunk, rep_lo, rep_hi, cores,
            paste(learners, collapse = "+"), if (SMOKE) "  [SMOKE]" else ""))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
# filename convention (spec): r13_te{000|030}_<scenario>_<rung>_chunk###.rds
tag <- if (SMOKE) "SMOKE_" else ""
te_tag <- sprintf("te%03d", as.integer(round(te * 100)))
fn <- file.path(OUT, sprintf("%sr13_%s_%s_%s_chunk%03d.rds", tag, te_tag, scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ------
# The ONLY change vs run_sim.R is the te_log_odds knob (headline uses 1.5).
pop <- make_population(scenario, model_type = "complex", te_log_odds = te,
                       pop_seed = POP_SEED, truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.6f (mc se %.6f) Psi_N=%.6f gap=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0,
            abs(pop$truth$gap_super_census)))

# ---- null truth check (te=0 only): |psi| <= 4*se_mc ; warn, do NOT stop ----
# At te=0 the GH blip gh_expit(m1)-gh_expit(m0) is identically 0 (m1 == m0
# bitwise when theta=0), so psi = se_mc = 0 exactly; max(., 1e-12) keeps the
# check well-defined in that degenerate case.
null_ok <- NA
if (te == 0) {
  null_ok <- abs(pop$truth$psi) <= max(4 * pop$truth$se_mc, 1e-12)
  cat(sprintf("[%s task %d] NULL CHECK: |psi|=%.3e vs 4*se_mc=%.3e -> null_ok=%s\n",
              RUN_ID, task, abs(pop$truth$psi), 4 * pop$truth$se_mc, null_ok))
  if (!isTRUE(null_ok))
    warning(sprintf("[%s task %d] null truth check FAILED: psi=%.6f, se_mc=%.6f",
                    RUN_ID, task, pop$truth$psi, pop$truth$se_mc))
}

# ---- one replication: mirrors run_sim.R::one_rep exactly + a te tag ----
# i is the GLOBAL rep index -> sample_seed = SAMPLE_SEED_BASE + i. draw_sample
# is independent of te, and in dgp.R theta enters ONLY the final Y draw (after
# the C and A RNG), so for a given (scenario, rep) the sampled units, C, A,
# weights and design are BYTE-IDENTICAL to the locked headline run (te=1.5);
# only the Y column differs. This preserves cross-run comparability.
one_rep <- function(i, pop, learners, te) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_estimators(obs, learners = learners)        # UNCHANGED, all five arms
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, te = te, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, te = te, est$diagnostics$drow, deff = dd$deff_clust,
                    icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk across cores (parallel-safe L'Ecuyer streams) -------
# Stream seed: SAMPLE_SEED_BASE + 1000L*cell_index + task. cell_index in 1..16
# and task in 1..160, so the 1000-per-cell offset guarantees distinct streams
# across cells/chunks WITHIN this run, and keeps them disjoint from the OLDER
# arc runs' SAMPLE_SEED_BASE + task range (R01/R03/R04/R07/R11/R12). The
# sibling batch runs (R14/R15/R16/R18/R19, and R02 via 1000L*m_idx + task)
# share this same formula, so cross-run iseed coincidences DO occur — harmless,
# because the worker stream only drives learner-internal randomness
# (SuperLearner CV folds / ranger / make_cf_folds); rep-level SAMPLES self-seed
# inside draw_sample(SAMPLE_SEED_BASE + global rep index) regardless of stream.
iseed <- SAMPLE_SEED_BASE + 1000L * cell_index + task
if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = iseed)
  # belt-and-suspenders: one_rep gets pop/learners/te as EXPLICIT parLapply args
  # (nothing is captured by closure), and we ALSO export them so any future
  # closure-style edit cannot ship broken (the R04 'workers lacked pop' lesson).
  clusterExport(cl, c("pop", "learners", "te", "scenario", "rung"), envir = environment())
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners, te = te)
} else {
  set.seed(iseed)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = learners, te = te)
}
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps in this chunk failed",
                                    RUN_ID, task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (full detail; private folder) ---------------------
out <- list(run_id = RUN_ID, te = te, scenario = scenario, rung = rung,
            learners = learners, chunk = chunk, reps = reps,
            Psi = Psi, truth = pop$truth, params = pop$params, null_ok = null_ok,
            results = res, diagnostics = diags)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----------
manifest <- list(
  run_id = RUN_ID, te = te, scenario = scenario, rung = rung,
  learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  null_ok = null_ok,                              # te=0 truth gate (NA at te=0.3)
  pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED, iseed_stream = iseed),
  run_dir = .this_dir,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  smoke = SMOKE
)
mf <- file.path(MAN_DIR, sprintf("%smanifest_%s_%s_%s_chunk%03d.rds",
                                 tag, te_tag, scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[%s task %d] manifest %s\n", RUN_ID, task, mf))

# ---- SMOKE: inline summary + explicit gate verdict ----------------------
if (SMOKE) {
  crit <- qt(0.975, pmax(1, res$df))
  res$reject <- abs(res$b) / res$se > crit          # test of H0: psi = 0
  res$cover  <- abs(res$b - Psi) <= crit * res$se   # CI covers true Psi
  summ <- do.call(rbind, by(res, res$method, function(d) data.frame(
    method = d$method[1], n_reps = nrow(d),
    mean_b = mean(d$b), emp_sd = sd(d$b), mean_se = mean(d$se),
    reject_rate = mean(d$reject), coverage = mean(d$cover),
    stringsAsFactors = FALSE)))
  rownames(summ) <- NULL
  ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
  summ <- summ[order(match(summ$method, ord_m)), ]
  np <- sapply(summ, is.numeric); sp <- summ; sp[np] <- round(sp[np], 4)
  cat(sprintf("\n[SMOKE] inline summary (te=%g %s %s; H0: psi=0; true Psi=%.6f):\n",
              te, scenario, rung, Psi))
  print(sp, row.names = FALSE)
  cat("(at te=0, coverage of Psi = 1 - reject_rate by construction)\n")

  finite_ok  <- all(is.finite(res$b)) && all(is.finite(res$se))
  n_meth     <- length(unique(res$method))
  methods_ok <- n_meth >= 4L && n_meth <= 5L       # L1_param has 4 (no CV arm on single-learner rungs)
  rr_dev     <- max(abs(summ$reject_rate - 0.05))
  rr_ok      <- rr_dev <= 0.25                     # lenient at 20 reps
  null_pass  <- isTRUE(null_ok)
  reasons <- paste(c(
    sprintf("b/se finite=%s", finite_ok),
    sprintf("method rows=%d (expect 4 at L1_param)", n_meth),
    sprintf("max|reject_rate-0.05|=%.3f (gate <=0.25 at 20 reps)", rr_dev),
    sprintf("null truth check=%s", null_pass)), collapse = "; ")
  if (finite_ok && methods_ok && rr_ok && null_pass) {
    cat(sprintf("\n[SMOKE-GATE] PASS: %s\n", reasons))
  } else {
    cat(sprintf("\n[SMOKE-GATE] STOP: %s\n", reasons))
  }
}
