# =====================================================================
# run.R  —  R03_isolation_2x2 driver (SLURM array job)
#
# Isolates CROSS-FITTING from DE-WEIGHTING by crossing the two factors the
# locked Fully-Aware-vs-CF contrast confounds, at a HARMONIZED 0.05 floor.
# Four arms per rep (see estimators_isolation.R):
#   SF-W single-fit/weighted | SF-U single-fit/UNWEIGHTED (dangerous tail)
#   CF-W cross-fit/weighted  | CF-U cross-fit/UNWEIGHTED
#
# Scenario: standard only. Rungs: L4_aggressive (the demonstration) +
# L1_param (placebo: all four arms should cover, since GLM is Donsker).
#
# Grid = expand.grid(rung in {L4_aggressive, L1_param}) x rep-chunks.
# One SLURM array task = one (rung x rep-chunk). Population built ONCE per task
# from POP_SEED; reps seeded by SAMPLE_SEED_BASE + GLOBAL rep index (identical
# to run_sim.R, so samples are byte-comparable to the locked run).
#
# NON-DESTRUCTIVE: sources the canonical engine read-only; writes per-task RDS
# to sim_output/arc_runs/R03_isolation_2x2/ and (driver-side, optional) never
# touches the locked results/sim_full_*.csv.
#
# Env vars (local-test fallbacks in parens):
#   SLURM_ARRAY_TASK_ID  1-based task index            (1)
#   SLURM_CPUS_PER_TASK  cores for this task           (detectCores-1)
#   SIM_N_REPS           total reps per rung           (1000)
#   SIM_CHUNK            reps per array task           (100)
#   SIM_CODE             canonical code dir            ("codes")
#   R03_OUT              per-task RDS dir              (REPO_ROOT/sim_output/arc_runs/R03_isolation_2x2)
#   R03_GFLOOR           harmonized propensity floor   (0.05)
#   SMOKE=1              tiny run: 50 reps, L4 only, 1 chunk (see NOTES.md)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "codes")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))      # exposes .sl/.eif_from_tmle/.se_des/make_cf_folds
source(file.path(CODE, "diagnostics.R"))     # deff_clust
source(file.path(CODE, "learners.R"))        # SL_LADDER
# this run's helper (lives in THIS folder; reuses the engine building blocks):
HERE <- Sys.getenv("R03_DIR",
                   file.path(REPO_ROOT, "codes", "arc_runs", "R03_isolation_2x2"))
source(file.path(HERE, "estimators_isolation.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
GFLOOR <- { x <- Sys.getenv("R03_GFLOOR"); if (nzchar(x)) as.numeric(x) else 0.05 }

OUT <- Sys.getenv("R03_OUT",
                  file.path(REPO_ROOT, "sim_output", "arc_runs", "R03_isolation_2x2"))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: STANDARD scenario, two rungs (demonstration + placebo) -------
# L4_aggressive first so the smoke / task-1 hits the decisive rung immediately.
RUNGS <- c("L4_aggressive", "L1_param")
SCEN  <- "standard"

if (SMOKE) {                                  # tiny, single-cell, finishes ~1-3 min
  N_REPS <- 50L; CHUNK <- 50L
  RUNGS  <- "L4_aggressive"                    # the dangerous rung only
  cores  <- min(cores, 6L)
}

cells <- expand.grid(scenario = SCEN, rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job <- grid[task, ]
scenario <- job$scenario; rung <- job$rung; chunk <- job$chunk
learners <- SL_LADDER[[rung]]
# GLOBAL rep index so chunks never collide and seeds match run_sim.R exactly.
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[R03 task %d] scenario=%s rung=%s chunk=%d reps=%d-%d g_floor=%.3f cores=%d%s\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, GFLOOR, cores,
            if (SMOKE) " [SMOKE]" else ""))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
fn <- file.path(OUT, sprintf("R03_%s_%s_chunk%03d%s.rds",
                             scenario, rung, chunk, if (SMOKE) "_SMOKE" else ""))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) -------------
pop <- make_population(scenario, model_type = "complex", truth_M = 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[R03 task %d] Psi=%.5f (mc se %.6f)\n", task, Psi, pop$truth$se_mc))

one_rep <- function(i, pop, learners, g_floor) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_isolation(obs, learners = learners, g_floor = g_floor)
  ch  <- attr(obs, "checks")
  # DEFF/ICC on the SF-W arm's EIF (analogue of run_sim's eif_fa audit)
  dd  <- deff_clust(est$diagnostics$eif_sfw, obs$strata, obs$cluster, obs$weight)
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
    CODE <- Sys.getenv("SIM_CODE", "codes")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    HERE <- Sys.getenv("R03_DIR",
                       file.path(REPO_ROOT, "codes", "arc_runs", "R03_isolation_2x2"))
    source(file.path(HERE, "estimators_isolation.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners, g_floor = GFLOOR)
} else {
  set.seed(SAMPLE_SEED_BASE + task)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = learners, g_floor = GFLOOR)
}
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored (e.g. degenerate learner fit)
if (!length(reps_out)) stop("all reps in this chunk failed")
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (non-clobbering; own folder) ----------------------------
out <- list(run = "R03_isolation_2x2", scenario = scenario, rung = rung,
            learners = learners, g_floor = GFLOOR, chunk = chunk, reps = reps,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags, smoke = SMOKE)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[R03 task %d] saved %s  (%d est rows, %d diag rows)\n",
            task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; reviewer-proofing) --------------
manifest <- list(
  run = "R03_isolation_2x2", scenario = scenario, rung = rung, learners = learners,
  g_floor = GFLOOR, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth, pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]])
mf <- file.path(OUT, sprintf("manifest_R03_%s_%s_chunk%03d%s.rds",
                             scenario, rung, chunk, if (SMOKE) "_SMOKE" else ""))
saveRDS(manifest, mf)
cat(sprintf("[R03 task %d] manifest %s\n", task, mf))

# ---- SMOKE: aggregate-on-the-spot + print the smoke-gate verdict -------------
# (full runs are aggregated separately by aggregate.R; see NOTES.md.)
if (SMOKE) {
  z_or_t <- function(df) qt(0.975, pmax(1, df))
  agg <- do.call(rbind, by(res, res$method, function(d) {
    crit <- z_or_t(d$df); cov <- mean(abs(d$b - Psi) <= crit * d$se)
    data.frame(method = d$method[1], n_reps = nrow(d),
               bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
               se_ratio = mean(d$se) / sd(d$b),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)))
  }))
  rownames(agg) <- NULL
  ord <- c("SF-W", "SF-U", "CF-W", "CF-U"); agg <- agg[order(match(agg$method, ord)), ]
  cat("\n==== R03 SMOKE summary (standard / ", rung, ", g_floor=", GFLOOR,
      ", ", nrow(res) / 4, " reps) ====\n", sep = "")
  num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
  print(ap, row.names = FALSE)
  sfu <- agg$coverage[agg$method == "SF-U"]
  cat(sprintf("\nSMOKE-GATE: SF-U (single-fit / UNWEIGHTED) L4 coverage = %.3f\n", sfu))
  if (length(sfu) && !is.na(sfu) && sfu >= 0.90) {
    cat("  >>> STOP-AND-REPORT: SF-U already covers ~nominal at L4. De-weighting,\n",
        "      not cross-fitting, may be rescuing coverage. Do NOT proceed to the\n",
        "      full run / appendix de-weighting justification before review.\n", sep = "")
  } else {
    cat("  >>> PROCEED: SF-U under-covers at L4 (as expected) -> cross-fitting is\n",
        "      the candidate active ingredient. Submit the full 1000-rep run.\n", sep = "")
  }
}
