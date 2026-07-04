# =====================================================================
# run.R  —  R20_cv_vs_cf_nondonsker driver  (NON-DESTRUCTIVE ARC run)
#
# SPEC: reviewer-driven follow-up to A4/A9a. The manuscript ladder makes
# the two NON-Donsker rungs SINGLE-learner (L1=SL.glm, L4=SL.ranger.deep)
# on purpose -- so SuperLearner's internal CV cannot hide overfitting
# behind a smooth ensemble member (see codes/learners.R header). A side
# effect: the Fully-Aware-CV arm (gated on length(learners) > 1 in
# estimators.R) is STRUCTURALLY ABSENT at L1 and L4, so the head-to-head
# "internal CV vs cross-fitting at a genuinely non-Donsker library" is
# never shown. This run closes that gap.
#
# GOAL: run ALL FIVE arms (incl. Fully-Aware-CV) at a MULTI-learner,
# NON-Donsker library
#   L5_nondonsker = c("SL.glm", "SL.earth", "SL.ranger.deep")
# (= L3's shape but with the AGGRESSIVE deep RF in place of the tame
# ranger.t1), on both designs, 1000 reps. Headline: does internal CV's
# ability to down-weight the overfitter substitute for cross-fitting, or
# does only PSU-level cross-fitting (CF) hold coverage at a non-Donsker
# library? Compare Fully-Aware-CV vs Fully-Aware-CF vs Fully-Aware
# (single-fit) coverage + se_ratio.
#
# NON-DESTRUCTIVE: SOURCES the canonical engine read-only. The library is
# defined HERE (not added to SL_LADDER), so no engine file is edited.
# run_estimators() is UNCHANGED; with a length-3 library it emits all five
# arms including the (weighted) internal-CV arm.
#
# GRID: scenario in {standard, R1} x library L5 = 2 cells, each x
# ceil(SIM_N_REPS/SIM_CHUNK) chunks.
#   FULL: SIM_N_REPS=1000, SIM_CHUNK=100 -> 2 x 10 = 20 tasks -> --array=1-20
#   tasks  1-10  standard L5_nondonsker  |  tasks 11-20  R1 L5_nondonsker
#
# Env vars (local fallbacks in brackets):
#   SLURM_ARRAY_TASK_ID  1-based task index                  [1]
#   SLURM_CPUS_PER_TASK  cores for this task                 [detectCores-1]
#   SIM_N_REPS           total reps per cell                 [1000]
#   SIM_CHUNK            reps per array task                 [100]
#   SIM_CODE             canonical code dir                  [codes]
#   R20_DIR              this run folder   [REPO_ROOT/codes/arc_runs/R20_cv_vs_cf_nondonsker]
#   R20_OUT              per-task RDS dir  [DATA_ROOT/arc_runs/R20_cv_vs_cf_nondonsker]
#   SMOKE                "1" -> tiny plumbing cell (standard, cheap glm+earth, 8 reps)
#   R20_SMOKE_DEEP       "1" -> smoke uses the REAL deep-RF L5 (slow; eyeball only)
# =====================================================================

RUN_ID <- "R20_cv_vs_cf_nondonsker"
RUNG_LAB <- "L5_nondonsker"

# ---- locate + source the canonical engine (read-only) ----------------
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))      # POP_SEED, SAMPLE_SEED_BASE, paths, RNGkind
source(file.path(CODE, "dgp.R"))         # make_population, draw_sample
source(file.path(CODE, "estimators.R"))  # run_estimators (UNCHANGED; CV arm fires at length>1)
source(file.path(CODE, "diagnostics.R")) # deff_clust, audit_population
source(file.path(CODE, "learners.R"))    # SL.ranger.deep wrapper + SL_LADDER
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

.this_dir <- Sys.getenv("R20_DIR", file.path(REPO_ROOT, "simulation", "enhancements", RUN_ID))

# ---- the non-Donsker MULTI-learner library (defined here, not in SL_LADDER) ----
L5 <- c("SL.glm", "SL.earth", "SL.ranger.deep")

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)

OUT <- Sys.getenv("R20_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
MAN_DIR <- file.path(OUT, "manifest")
for (d in c(OUT, MAN_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x library(L5) (scenario-major) ---------------
cells <- data.frame(scenario = c("standard", "R1"), stringsAsFactors = FALSE)
LEARNERS <- L5

if (SMOKE) {
  cells <- cells[cells$scenario == "standard", , drop = FALSE]
  # default smoke uses a CHEAP multi-learner library so the CV gate still
  # fires (length>1) and all five arms emit, but it finishes fast. The
  # deep-RF science path is validated by R02/R03/R16; set R20_SMOKE_DEEP=1
  # to eyeball the real L5 (slow: minutes per rep).
  LEARNERS <- if (identical(Sys.getenv("R20_SMOKE_DEEP"), "1")) L5
              else c("SL.glm", "SL.earth")
  N_REPS <- 8L; CHUNK <- 8L; cores <- min(cores, 2L); task <- 1L
  cat(sprintf("[SMOKE] %s: standard x {%s}, %d reps, %d cores\n",
              RUN_ID, paste(LEARNERS, collapse = "+"), N_REPS, cores))
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))        # FULL grid = 2 cells x 10 chunks = 20 tasks
job        <- grid[task, ]
scenario   <- job$scenario
cell_index <- job$cell_index
chunk      <- job$chunk
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[%s task %d] scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d learners=%s%s\n",
            RUN_ID, task, scenario, RUNG_LAB, chunk, rep_lo, rep_hi, cores,
            paste(LEARNERS, collapse = "+"), if (SMOKE) "  [SMOKE]" else ""))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn <- file.path(OUT, sprintf("%sr20_%s_chunk%03d.rds", tag, scenario, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic; headline te_log_odds=1.5) ----
pop <- make_population(scenario, model_type = "complex",
                       pop_seed = POP_SEED, truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.6f (mc se %.6f) Psi_N=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0))

# ---- one replication: mirrors run_sim.R::one_rep; all five arms at L5 ----
# i is the GLOBAL rep index -> sample_seed = SAMPLE_SEED_BASE + i, so for a
# given (scenario, rep) the sampled units/C/A/weights/design are BYTE-IDENTICAL
# to the locked headline run; only the LIBRARY differs (L5 vs the L1-L4 ladder).
one_rep <- function(i, pop, learners) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_estimators(obs, learners = learners)   # length>1 -> all five arms incl. CV
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow, deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk across cores (parallel-safe L'Ecuyer streams) -------
# Stream seed SAMPLE_SEED_BASE + 1000*cell_index + task: unique WITHIN this run
# (cell_index in 1..2, task in 1..20) and disjoint from the OLDER runs'
# SAMPLE_SEED_BASE + task range. Sibling batch runs share this formula
# (harmless: streams drive only learner-internal RNG/folds; samples self-seed
# in draw_sample from the global rep index).
iseed <- SAMPLE_SEED_BASE + 1000L * cell_index + task
if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger); library(earth) })
  })
  clusterSetRNGStream(cl, iseed = iseed)
  clusterExport(cl, c("pop", "scenario"), envir = environment())   # learners passed explicitly
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = LEARNERS)
} else {
  set.seed(iseed)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = LEARNERS)
}
reps_out <- Filter(Negate(is.null), reps_out)
n_failed <- length(reps) - length(reps_out)
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps in this chunk failed",
                                    RUN_ID, task, length(reps)))
if (n_failed > 0)
  cat(sprintf("[%s task %d] WARNING: %d/%d reps failed; chunk RDS will checkpoint as complete -- delete %s to recompute\n",
              RUN_ID, task, n_failed, length(reps), fn))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
res$rung <- RUNG_LAB
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (full detail; private folder) ---------------------
out <- list(run_id = RUN_ID, scenario = scenario, rung = RUNG_LAB,
            learners = LEARNERS, chunk = chunk, reps = reps, n_failed = n_failed,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags)
saveRDS(out, fn)
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----------
manifest <- list(
  run_id = RUN_ID, scenario = scenario, rung = RUNG_LAB, learners = LEARNERS,
  chunk = chunk, reps = reps, n_failed = n_failed,
  params = pop$params, truth = pop$truth,
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
mf <- file.path(MAN_DIR, sprintf("%smanifest_%s_chunk%03d.rds", tag, scenario, chunk))
saveRDS(manifest, mf)
cat(sprintf("[%s task %d] manifest %s\n", RUN_ID, task, mf))

# ---- SMOKE: inline summary + explicit gate verdict ----------------------
if (SMOKE) {
  crit <- qt(0.975, pmax(1, res$df))
  res$cover <- abs(res$b - Psi) <= crit * res$se
  summ <- do.call(rbind, by(res, res$method, function(d) data.frame(
    method = d$method[1], n_reps = nrow(d),
    bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
    coverage = mean(d$cover), stringsAsFactors = FALSE)))
  rownames(summ) <- NULL
  ord_m <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
  summ <- summ[order(match(summ$method, ord_m)), ]
  np <- sapply(summ, is.numeric); sp <- summ; sp[np] <- round(sp[np], 4)
  cat(sprintf("\n[SMOKE] inline summary (standard %s; true Psi=%.6f):\n",
              paste(LEARNERS, collapse = "+"), Psi))
  print(sp, row.names = FALSE)

  finite_ok <- all(is.finite(res$b)) && all(is.finite(res$se))
  meths     <- unique(res$method)
  cv_ok     <- "Fully-Aware-CV" %in% meths            # THE key check: CV arm fires at multi-learner
  cf_ok     <- "Fully-Aware-CF" %in% meths
  n_meth_ok <- length(meths) == 5L
  reasons <- paste(c(
    sprintf("b/se finite=%s", finite_ok),
    sprintf("5 arms present=%s", n_meth_ok),
    sprintf("CV arm present=%s (the point of this run)", cv_ok),
    sprintf("CF arm present=%s", cf_ok)), collapse = "; ")
  if (finite_ok && n_meth_ok && cv_ok && cf_ok)
    cat(sprintf("\n[SMOKE-GATE] PASS: %s\n", reasons))
  else
    cat(sprintf("\n[SMOKE-GATE] STOP: %s\n", reasons))
}
