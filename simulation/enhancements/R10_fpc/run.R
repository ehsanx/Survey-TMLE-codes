# =====================================================================
# R10_fpc/run.R  —  FPC / first-stage-fraction sensitivity (estimand/regime defense)
#
# NON-DESTRUCTIVE enhancement run. SOURCES the canonical engine (R/*.R) and
# REUSES its building blocks; writes only into this run's own output folders.
# Nothing in R/*.R is modified.
#
# Two analyses (selected by the SLURM array task index; see grid below):
#
#  PART A — FPC branch (over-coverage / m/M_pop->0 defense).
#    On the R1 design (H=50, J_per_stratum=8 PSUs, base_m=2 -> first-stage
#    fraction f = m/M = 2/8 = 0.25) and ONLY for rungs L1_param and L2_smooth,
#    recompute the Fully-Aware-CF design SE WITH a finite-population correction
#    (fpc = population PSUs-per-stratum = J_per_stratum). Compare CF coverage
#    WITH vs WITHOUT fpc. We DO NOT apply fpc at L4_aggressive: L4-CF already
#    sits at ~0.985 (over-coverage); the fpc shrinks the SE, which would tip
#    L4 into UNDER-coverage and is therefore the wrong place to defend the
#    finite-population claim. We also carry the Fully-Aware (single weighted
#    fit) arm fpc-vs-no-fpc for context (the eif_fa is returned by the engine).
#
#  PART B — first-stage-fraction sweep (m/M_pop -> 0 insensitivity).
#    Grow the finite population's PSUs-per-stratum (J_per_stratum) while holding
#    base_m = 2 fixed, so the realized first-stage fraction f = base_m / J shrinks
#    8/8=... actually 2/J: J in {8,16,32,64} -> f in {0.25,0.125,0.0625,0.03125}.
#    M_per_psu is held at the R1 default (300) so the population N grows with J.
#    For each J we run the CANONICAL (no-fpc) Fully-Aware-CF arm and check that
#    super-population coverage stays flat as f->0 (the m/M_pop->0 assumption is
#    innocuous because the super-pop estimand and its EIF do not depend on f).
#    Done at rung L1_param only (cheap; the regime point is about f, not the
#    learner). The truth Psi is recomputed per J (population changes), but the
#    super-pop truth is essentially J-invariant by construction.
#
# Engine contract reused (verified against R/*.R, 2026-06-07):
#   make_population(scenario,model_type,...,J_per_stratum=,M_per_psu=) -> pop
#   draw_sample(pop, sample_seed=, model_type=)  [base_m defaults: R1 -> 2]
#   run_estimators(obs, learners=) -> results(method,b,se,df) + diagnostics(eif_fa,...)
#   make_cf_folds(), .sl(), .eif_from_tmle(), .se_des(), deff_clust()  [building blocks]
#
# Env vars:
#   SIM_CODE              code dir (default 'R')
#   REPO_ROOT/DATA_ROOT   inherited via config.R
#   SLURM_ARRAY_TASK_ID   1-based task index (default 1)
#   SLURM_CPUS_PER_TASK   cores (default detectCores-1)
#   FPC_N_REPS            reps per cell (default 1000)
#   SMOKE                 '1' -> tiny validation run (20 reps, 1 part, 1 J level)
# =====================================================================

suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))

# run-local helpers (this folder)
RUN_DIR <- Sys.getenv("FPC_RUN_DIR",
                      file.path(REPO_ROOT, "simulation", "enhancements", "R10_fpc"))
source(file.path(RUN_DIR, "fpc_helpers.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE  <- Sys.getenv("SMOKE", "0") == "1"
task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("FPC_N_REPS", 1000L)

ID  <- "R10_fpc"
OUT <- file.path(DATA_ROOT, "arc_runs", ID)          # per-task RDS (NON-clobbering)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# TASK GRID
#   PART A: 1 task  (R1 design; rungs L1_param + L2_smooth; fpc vs no-fpc)
#   PART B: one task per J-sweep level (rung L1_param; canonical no-fpc CF)
# --------------------------------------------------------------------------
J_SWEEP <- c(8L, 16L, 32L, 64L)                      # base_m=2 -> f = 2/J
A_rungs <- c("L1_param", "L2_smooth")                # FPC applied here ONLY (NOT L4)

grid <- rbind(
  data.frame(part = "A", J = NA_integer_, rung = NA_character_, stringsAsFactors = FALSE),
  data.frame(part = "B", J = J_SWEEP,     rung = "L1_param",    stringsAsFactors = FALSE)
)

if (SMOKE) {
  N_REPS  <- 20L
  cores   <- min(cores, 4L)
  J_SWEEP <- 8L
  grid    <- grid[grid$part == "A", , drop = FALSE]   # smoke: just PART A, 20 reps
  cat("[SMOKE] PART A only, 20 reps, R1 / {L1_param,L2_smooth}, fpc vs no-fpc\n")
}
stopifnot(task >= 1L, task <= nrow(grid))
job  <- grid[task, ]
part <- job$part

cat(sprintf("[%s task %d] part=%s J=%s rung=%s N_REPS=%d cores=%d\n",
            ID, task, part, ifelse(is.na(job$J), "-", job$J),
            ifelse(is.na(job$rung), "-", job$rung), N_REPS, cores))

# ---- per-task output filename + checkpoint: skip instantly if already done ----
# Computed EARLY (before any population build / parLapply) so a re-submitted
# --array idempotently skips finished tasks. One file per part (Part A: one
# combined RDS; Part B: one RDS per J-sweep level).
if (part == "A") {
  fn <- file.path(OUT, sprintf("%s_partA%s.rds", ID, if (SMOKE) "_SMOKE" else ""))
} else {
  fn <- file.path(OUT, sprintf("%s_partB_J%03d%s.rds", ID, as.integer(job$J),
                               if (SMOKE) "_SMOKE" else ""))
}
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- parallel cluster bootstrap (workers re-source engine + helpers) ----------
mk_cluster <- function() {
  cl <- makeCluster(cores)
  clusterExport(cl, c("CODE", "RUN_DIR"), envir = .GlobalEnv)
  clusterEvalQ(cl, {
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    source(file.path(RUN_DIR, "fpc_helpers.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  cl
}

# ===========================================================================
# PART A — FPC vs no-FPC on R1 (L1_param + L2_smooth)
# ===========================================================================
# For each rep we run the canonical engine once per rung (so the no-fpc CF SE is
# byte-identical to the locked pipeline), then recompute the CF EIF locally and
# apply BOTH .se_des (no fpc, == engine) and .se_des_fpc (with fpc). We also carry
# the Fully-Aware single-fit EIF (engine's eif_fa) through both SE formulas.
one_rep_A <- function(i, pop, rungs) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  # population PSUs-per-stratum, broadcast to each observed row (constant w/in stratum):
  Jpop <- pop$params$J_per_stratum
  fpc1 <- rep(Jpop, nrow(obs))

  out <- list()
  for (rung in rungs) {
    learners <- SL_LADDER[[rung]]
    est <- run_estimators(obs, learners = learners)

    # --- Fully-Aware (single weighted fit): eif_fa comes straight from engine ---
    eif_fa <- est$diagnostics$eif_fa
    b_fa   <- est$results$b[est$results$method == "Fully-Aware"]
    se_fa_nofpc <- .se_des(eif_fa, obs$strata, obs$cluster, obs$weight, clustered = TRUE)
    se_fa_fpc   <- .se_des_fpc(eif_fa, obs$strata, obs$cluster, obs$weight, fpc1)

    # --- Fully-Aware-CF: recompute the CF EIF locally, apply both SE formulas ---
    cfr <- cf_eif_for_rep(obs, learners)
    se_cf_nofpc <- .se_des(cfr$eif, obs$strata, obs$cluster, obs$weight, clustered = TRUE)
    se_cf_fpc   <- .se_des_fpc(cfr$eif, obs$strata, obs$cluster, obs$weight, fpc1)

    # sanity: locally-recomputed no-fpc CF SE should match the engine's CF SE.
    se_cf_engine <- est$results$se[est$results$method == "Fully-Aware-CF"]

    out[[rung]] <- data.frame(
      rep = i, rung = rung,
      method = rep(c("Fully-Aware", "Fully-Aware-CF"), each = 2),
      fpc    = rep(c("no", "yes"), times = 2),
      b      = c(b_fa, b_fa, cfr$psi, cfr$psi),
      se     = c(se_fa_nofpc$se, se_fa_fpc$se, se_cf_nofpc$se, se_cf_fpc$se),
      df     = c(se_fa_nofpc$df, se_fa_fpc$df, se_cf_nofpc$df, se_cf_fpc$df),
      se_cf_engine_check = se_cf_engine,   # diagnostic: == se for CF/no row
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# ===========================================================================
# PART B — first-stage-fraction sweep (canonical no-fpc CF, rung L1_param)
# ===========================================================================
one_rep_B <- function(i, pop, learners) {
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_estimators(obs, learners = learners)
  ch  <- attr(obs, "checks")
  r   <- est$results[est$results$method == "Fully-Aware-CF", c("method", "b", "se", "df")]
  data.frame(rep = i, r, n = ch$n, df_design = ch$df_design,
             w_cv = ch$w_cv, stringsAsFactors = FALSE)
}

# ===========================================================================
# DRIVE
# ===========================================================================
manifest_common <- function(extra = list()) {
  c(list(
    id = ID, part = part, n_reps = N_REPS,
    seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, TRUTH_SEED = TRUTH_SEED),
    R_version = R.version.string,
    packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                      function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    git = arc_git_sha(),
    sysname = Sys.info()[["nodename"]]
  ), extra)
}

reps <- seq_len(N_REPS)

if (part == "A") {
  # ---- R1 population built ONCE (canonical geometry) ----
  pop <- make_population("R1", model_type = "complex", truth_M = 2e6L)
  Psi <- pop$truth$psi
  cat(sprintf("[%s] R1 Psi=%.5f (mc se %.6f)  J/str=%d base_m=2 -> f=%.4f\n",
              ID, Psi, pop$truth$se_mc, pop$params$J_per_stratum,
              2 / pop$params$J_per_stratum))

  cl <- mk_cluster(); on.exit(stopCluster(cl), add = TRUE)
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 9100L + task)
  rows <- do.call(rbind, parLapply(cl, reps, one_rep_A, pop = pop, rungs = A_rungs))

  # per-task RDS
  out <- list(id = ID, part = "A", scenario = "R1", Psi = Psi,
              truth = pop$truth, params = pop$params, per_rep = rows,
              first_stage_fraction = 2 / pop$params$J_per_stratum)
  saveRDS(out, fn)   # `fn` computed early for the checkpoint guard

  # summary CSV (one row per rung x method x fpc), run_sim conventions
  summ <- do.call(rbind, lapply(split(rows, list(rows$rung, rows$method, rows$fpc),
                                      drop = TRUE), function(d) {
    summarize_arm(d$b, d$se, d$df, Psi, "R1", d$rung[1], d$method[1], d$fpc[1])
  }))
  rownames(summ) <- NULL
  ord <- order(summ$rung, summ$method, summ$fpc)
  summ <- summ[ord, ]
  num  <- sapply(summ, is.numeric); summ_p <- summ; summ_p[num] <- round(summ_p[num], 4)

  arc_csv_dir <- file.path(RESULTS_DIR, "arc"); dir.create(arc_csv_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(arc_csv_dir, sprintf("%s_partA_summary%s.csv", ID, if (SMOKE) "_SMOKE" else ""))
  write.csv(summ_p, csv, row.names = FALSE)
  cat(sprintf("[%s] PART A saved %s\n  summary -> %s\n", ID, fn, csv))
  print(summ_p, row.names = FALSE)

  mf <- file.path(MANIFEST_DIR, sprintf("manifest_%s_partA%s.rds", ID, if (SMOKE) "_SMOKE" else ""))
  saveRDS(manifest_common(list(scenario = "R1", rungs = A_rungs, Psi = Psi,
                               params = pop$params, truth = pop$truth)), mf)
  cat(sprintf("[%s] manifest %s\n", ID, mf))

} else {                                                # part == "B"
  Jval     <- as.integer(job$J)
  learners <- SL_LADDER[[job$rung]]
  frac     <- 2 / Jval

  # ---- grow the finite population: override J_per_stratum (M_per_psu held at R1=300)
  #      so the first-stage fraction f = base_m/J shrinks toward 0 as J grows. ----
  pop <- make_population("R1", model_type = "complex",
                         J_per_stratum = Jval, M_per_psu = 300L, truth_M = 2e6L)
  Psi <- pop$truth$psi
  cat(sprintf("[%s] sweep J=%d (M/PSU=300, N=%d)  Psi=%.5f  f=base_m/J=2/%d=%.5f\n",
              ID, Jval, pop$params$N_pop, Psi, Jval, frac))

  cl <- mk_cluster(); on.exit(stopCluster(cl), add = TRUE)
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 9200L + task)
  rows <- do.call(rbind, parLapply(cl, reps, one_rep_B, pop = pop, learners = learners))
  rows$J <- Jval; rows$frac <- frac

  out <- list(id = ID, part = "B", scenario = "R1", J = Jval, frac = frac,
              rung = job$rung, Psi = Psi, truth = pop$truth, params = pop$params,
              per_rep = rows)
  saveRDS(out, fn)   # `fn` computed early for the checkpoint guard

  summ <- summarize_arm(rows$b, rows$se, rows$df, Psi, "R1", job$rung,
                        "Fully-Aware-CF", fpc = "no")
  summ$J <- Jval; summ$frac <- frac; summ$N_pop <- pop$params$N_pop
  num <- sapply(summ, is.numeric); summ_p <- summ; summ_p[num] <- round(summ_p[num], 4)
  cat(sprintf("[%s] PART B J=%d saved %s\n", ID, Jval, fn))
  print(summ_p, row.names = FALSE)

  mf <- file.path(MANIFEST_DIR, sprintf("manifest_%s_partB_J%03d%s.rds", ID, Jval, if (SMOKE) "_SMOKE" else ""))
  saveRDS(manifest_common(list(scenario = "R1", J = Jval, frac = frac, rung = job$rung,
                               Psi = Psi, params = pop$params, truth = pop$truth)), mf)
  cat(sprintf("[%s] manifest %s\n", ID, mf))
}

cat(sprintf("[%s task %d] DONE %s\n", ID, task, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
