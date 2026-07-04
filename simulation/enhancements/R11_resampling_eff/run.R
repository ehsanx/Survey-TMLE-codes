# =====================================================================
# run.R  —  ARC run R11_resampling_eff  (Web Appendix D)
#
# (a) RESAMPLING-SE CHECK: corroborate the design-linearization (Taylor) SE
#     that the engine returns against a replication-based survey JACKKNIFE
#     (JKn, delete-one-PSU) SE computed on the SAME influence function, for
#     the Fully-Aware (single-fit) and Fully-Aware-CF (cross-fitted) arms.
#     They should agree (ratio ~ 1) -> defends Eq 8.
# (b) EFFICIENCY BENCHMARK: a survey-weighted IPW / g-computation (AIPW)
#     baseline (svyglm nuisances) scored with the SAME design-EIF machinery,
#     so coverage and CI width are directly comparable to cross-fitted TMLE.
#
# RUN ORDER (per spec): L1_param FIRST (where linearization and jackknife
# MUST agree). The SLURM grid lists rungs in the order (L1_param, L3_adaptive)
# so the first array tasks are L1; only extend to L3 if L1 agrees.
#
# NON-DESTRUCTIVE: sources the canonical engine read-only; all new logic is in
# helpers_R11.R. Per-task RDS -> sim_output/arc_runs/R11_resampling_eff/.
# Summary CSV -> results/R11_resampling_eff_summary.csv. Nothing clobbers
# the locked results.
#
# Env vars (local-test fallbacks in brackets):
#   SIM_CODE              canonical code dir                [R]
#   SLURM_ARRAY_TASK_ID   1-based task index                [1]
#   SLURM_CPUS_PER_TASK   cores for this task               [detectCores-1]
#   R11_N_REPS            total reps per (scenario,rung)    [500]
#   R11_CHUNK             reps per array task               [100]
#   R11_OUT               per-task RDS dir   [REPO_ROOT/sim_output/arc_runs/R11_resampling_eff]
#   SMOKE                 "1" -> tiny run (20 reps, 1 cell, serial)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# run-specific helpers (jackknife-on-EIF + IPW/svyglm baseline)
HERE <- Sys.getenv("R11_DIR",
                   file.path(REPO_ROOT, "simulation", "enhancements", "R11_resampling_eff"))
source(file.path(HERE, "helpers_R11.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

# ---- run grid: ONLY the two rungs the spec calls for, L1 listed FIRST -------
# (L1_param where Taylor==jackknife MUST hold; L3_adaptive as the extension.)
R11_RUNGS <- c("L1_param", "L3_adaptive")
cells <- expand.grid(scenario = c("standard", "R1"), rung = R11_RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
# order rows so ALL L1 cells precede ALL L3 cells (rung-major), L1 first
cells <- cells[order(match(cells$rung, R11_RUNGS), cells$scenario), ]
rownames(cells) <- NULL

if (SMOKE) {
  # tiny: one cell (standard / L1_param), 20 reps, serial, finishes ~1-3 min
  N_REPS <- 20L; CHUNK <- 20L; cores <- 1L
  cells <- cells[cells$scenario == "standard" & cells$rung == "L1_param", , drop = FALSE]
  task  <- 1L
} else {
  N_REPS <- geti("R11_N_REPS", 500L)
  CHUNK  <- geti("R11_CHUNK", 100L)
  cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
  task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], chunk = seq_len(n_chunks))))
rownames(grid) <- NULL
stopifnot(task >= 1L, task <= nrow(grid))
job <- grid[task, ]
scenario <- job$scenario; rung <- job$rung; chunk <- job$chunk
learners <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

OUT <- Sys.getenv("R11_OUT",
                  file.path(REPO_ROOT, "sim_output", "arc_runs", "R11_resampling_eff"))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[R11 task %d] scenario=%s rung=%s chunk=%d reps=%d-%d cores=%d SMOKE=%s\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, cores, SMOKE))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
fn <- file.path(OUT, sprintf("R11_%s_%s_chunk%03d.rds", scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ------------
pop <- make_population(scenario, model_type = "complex", truth_M = 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[R11 task %d] Psi=%.5f (mc se %.6f)\n", task, Psi, pop$truth$se_mc))

# ---- one replication --------------------------------------------------------
# Reuses run_estimators() to get the engine's TMLE arms + their EIFs, then:
#   * jackknife SE on the Fully-Aware and Fully-Aware-CF EIFs (resampling check)
#   * IPW/svyglm AIPW baseline (efficiency benchmark)
# Returns engine-style rows: method, b, se (Taylor), df, plus se_jk/df_jk on
# the two TMLE arms whose EIF we jackknife.
one_rep <- function(i, pop, learners) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_estimators(obs, learners = learners)        # canonical TMLE arms
  ch  <- attr(obs, "checks")

  # EIF of the single-fit Fully-Aware arm is returned in diagnostics$eif_fa.
  # The CF arm's EIF is not returned by run_estimators(); recompute it here
  # from the SAME cross-fitted nuisances using the exposed building blocks so
  # the jackknife sees the identical influence function the engine scored.
  eif_fa <- est$diagnostics$eif_fa

  # --- reconstruct the Fully-Aware-CF EIF (mirrors estimators.R CF arm) ---
  W_cols <- c("L1", "L2", "L3", "L4")
  d  <- obs; n <- nrow(d); w <- d$weight
  W  <- d[, W_cols, drop = FALSE]
  XA <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)
  fold <- make_cf_folds(d$strata, d$cluster, V = 5L)
  Q0o <- Q1o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = NULL, learners)   # UNWEIGHTED (ignorability)
    gf <- .sl(d$A[tr], W[tr, ],  weights = NULL, learners)
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g1o <- pmin(pmax(g1o, 0.05), 1 - 0.05)
  Qoo <- pmin(pmax(cbind(Q0o, Q1o), 1e-3), 1 - 1e-3)
  cf  <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
              Q = Qoo, g1W = g1o)
  e_cf <- .eif_from_tmle(cf, d$Y, d$A, w)$eif

  # --- (a) jackknife (JKn) SE on each influence function ---
  jk_fa <- jk_se_on_eif(eif_fa, d$strata, d$cluster, w)
  jk_cf <- jk_se_on_eif(e_cf,   d$strata, d$cluster, w)

  # --- (b) IPW/svyglm AIPW baseline (efficiency benchmark) ---
  ipw <- ipw_svyglm_ate(obs, W_cols = W_cols, g_bound = 0.05)

  # engine TMLE rows (Taylor SE) + jackknife columns where applicable
  eng <- est$results
  eng$se_jk <- NA_real_; eng$df_jk <- NA_real_
  eng$se_jk[eng$method == "Fully-Aware"]    <- jk_fa$se_jk
  eng$df_jk[eng$method == "Fully-Aware"]    <- jk_fa$df_jk
  eng$se_jk[eng$method == "Fully-Aware-CF"] <- jk_cf$se_jk
  eng$df_jk[eng$method == "Fully-Aware-CF"] <- jk_cf$df_jk

  # IPW row (Taylor SE only; no jackknife requested for the baseline)
  ipw_row <- ipw$row; ipw_row$se_jk <- NA_real_; ipw_row$df_jk <- NA_real_

  results <- rbind(eng, ipw_row)

  cbind(rep = i, results,
        n = ch$n, df_design = ch$df_design, w_cv = ch$w_cv)
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk (parallel on ARC, serial in SMOKE) -----------------------
if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    HERE <- Sys.getenv("R11_DIR",
                       file.path(REPO_ROOT, "simulation", "enhancements", "R11_resampling_eff"))
    source(file.path(HERE, "helpers_R11.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + task)
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners)
} else {
  set.seed(SAMPLE_SEED_BASE + task)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = learners)
}
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored (e.g. degenerate learner fit)
if (!length(reps_out)) stop("all reps in this chunk failed")
res <- do.call(rbind, reps_out)

# ---- per-task output --------------------------------------------------------
out <- list(scenario = scenario, rung = rung, learners = learners, chunk = chunk,
            reps = reps, Psi = Psi, truth = pop$truth, params = pop$params,
            results = res)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[R11 task %d] saved %s  (%d rows)\n", task, fn, nrow(res)))

# ---- reproducibility manifest (one per task; reviewer-proofing) -------------
manifest <- list(
  run = "R11_resampling_eff", scenario = scenario, rung = rung,
  learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)),
                                         error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]]
)
mf <- file.path(OUT, sprintf("manifest_%s_%s_chunk%03d.rds", scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[R11 task %d] manifest %s\n", task, mf))
