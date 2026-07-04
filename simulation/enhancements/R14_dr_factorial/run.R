# =====================================================================
# run.R  —  R14_dr_factorial driver (SLURM array job)
#   Double-robustness factorial incl. the informative-sampling
#   unweighted-nuisance threat (spec item A6,
#   Writing/comments/phase4-arc-sim-specs.md).
#
# GOAL: the paper claims double robustness but only ever misspecifies BOTH
# nuisances jointly. This run crosses the misspecification PER NUISANCE:
#   q_spec in {C,W} x g_spec in {C,W}      (C = correct/raw frame, W = wrong/
#                                            Kang-Schafer frame; dr_helpers.R)
#   x sampling in {info, noninfo}          (draw_sample oversample = TRUE/FALSE)
#   x rung in {L1_param, L2_smooth}        (GLM vs smooth SL library)
# = 16 cells, three arms each (FA-w single-fit weighted; CF-u unweighted OOF,
# the paper default; CF-w weighted OOF on the SAME folds). Two headline reads:
#   (1) the 2x2 DR table at L1/noninfo  -- both-correct AND each one-correct
#       cell ~unbiased (double robustness), both-wrong biased;
#   (2) the THREAT view: in g-wrong / both-wrong x INFORMATIVE cells, the
#       CF-u vs CF-w bias gap (unweighted misspecified library converging to
#       the SAMPLE projection -- the Gemini mechanism); the gap should vanish
#       under noninfo (weights equal by construction, w_cv ~ 0).
#
# NON-DESTRUCTIVE: sources the canonical engine read-only; ALL new logic is in
# dr_helpers.R (this folder). Outputs go to a private arc_runs subtree +
# results/arc/; nothing clobbers the locked results.
#
# One SLURM array task = one (cell x rep-chunk).
#   FULL: 16 cells x ceil(1000/100)=10 chunks  ->  --array=1-160
#   (grid is CELL-MAJOR: tasks 1-10 = cell 1 chunks 1-10, ..., 151-160 =
#    cell 16; cells ordered by expand.grid with q_spec fastest, then g_spec,
#    sampling, rung -- see `cells` below.)
#
# Env vars (with local-test fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores for this task           (default detectCores-1)
#   SIM_N_REPS           total reps per cell           (default 1000)
#   SIM_CHUNK            reps per array task           (default 100)
#   SIM_CODE             canonical code dir            (default "codes")
#   R14_DIR              this run's folder  (default codes/arc_runs/R14_dr_factorial)
#   R14_OUT              per-task RDS dir   (default DATA_ROOT/arc_runs/R14_dr_factorial)
#   SMOKE                "1" -> tiny 1-cell local validation (20 reps, 2 cores)
# =====================================================================

RUN_ID <- "R14_dr_factorial"

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))        # POP_SEED, SAMPLE_SEED_BASE, paths, RNGkind
source(file.path(CODE, "dgp.R"))           # make_population, draw_sample, apply_L
source(file.path(CODE, "estimators.R"))    # .sl, .eif_from_tmle, .se_des, make_cf_folds
source(file.path(CODE, "diagnostics.R"))   # deff_clust, audit_population
source(file.path(CODE, "learners.R"))      # SL_LADDER
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

# run-local helpers (the ONLY new logic); resolve next to this file
.this_dir <- Sys.getenv("R14_DIR",
  file.path(REPO_ROOT, "simulation", "enhancements", "R14_dr_factorial"))
source(file.path(.this_dir, "dr_helpers.R"))

suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)

OUT <- Sys.getenv("R14_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: q_spec x g_spec x sampling x rung (16 cells) ----------------
# expand.grid order (q_spec fastest): cell 1 = (C,C,info,L1_param),
# 2 = (W,C,info,L1), 3 = (C,W,info,L1), 4 = (W,W,info,L1), 5-8 = same x noninfo,
# 9-16 = the same eight for L2_smooth.
cells <- expand.grid(q_spec   = c("C", "W"),
                     g_spec   = c("C", "W"),
                     sampling = c("info", "noninfo"),
                     rung     = c("L1_param", "L2_smooth"),
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
cells$cell_index <- seq_len(nrow(cells))   # stable 1..16 id (used in the RNG stream)

if (SMOKE) {
  # tiny validation: ONE cheap both-correct cell (q=C, g=C, info, L1_param),
  # 20 reps, 2 cores -> ~2-3 min. This is full-grid cell_index 1.
  cells <- cells[cells$q_spec == "C" & cells$g_spec == "C" &
                 cells$sampling == "info" & cells$rung == "L1_param", , drop = FALSE]
  N_REPS <- 20L; CHUNK <- 20L; cores <- min(cores, 2L); task <- 1L
  cat(sprintf("[SMOKE] %s: qC/gC/info/L1_param, 20 reps, %d cores\n", RUN_ID, cores))
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))
job        <- grid[task, ]
q_spec     <- job$q_spec
g_spec     <- job$g_spec
sampling   <- job$sampling
rung       <- job$rung
chunk      <- job$chunk
cell_index <- job$cell_index
oversample <- sampling == "info"           # draw_sample's informativeness toggle
learners   <- SL_LADDER[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[%s task %d] q=%s g=%s sampling=%s rung=%s chunk=%d reps=%d-%d cores=%d learners=%s\n",
            RUN_ID, task, q_spec, g_spec, sampling, rung, chunk, rep_lo, rep_hi,
            cores, paste(learners, collapse = "+")))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr14_q%s_g%s_%s_%s_chunk%03d.rds",
                              tag, q_spec, g_spec, sampling, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) ------------
# model_type = "simple": apply_L is the identity, so obs$L1..L4 carry the RAW
# C1..C4 (the misspecification is applied per-nuisance in dr_helpers.R). The
# population data and truth Psi are IDENTICAL across model_type (R01 established
# this: model_type only changes what the learners SEE, via apply_L at sampling).
pop <- make_population("standard", model_type = "simple",
                       truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0,
            abs(pop$truth$gap_super_census)))

# ---- one replication: sample -> three DR arms -> diagnostics ----------------
# sample_seed = SAMPLE_SEED_BASE + i with i the GLOBAL rep index, and
# draw_sample is internally .scoped_seed-ed -> for the info cells the sampled
# ROWS are byte-identical to the locked headline run's rep i (same standard
# population; model_type only changes the L columns via apply_L, not the
# selection), preserving cross-run comparability. noninfo cells differ only
# through oversample = FALSE (equal stratum sampling fractions -> equal
# weights, w_cv ~ 0: non-informative by construction).
one_rep <- function(i, pop, learners, q_spec, g_spec, oversample) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i,
                     model_type = "simple", oversample = oversample)
  est <- run_dr_arms(obs, learners = learners, q_spec = q_spec, g_spec = g_spec,
                     V_cf = 5L, g_floor = 0.05)
  ch  <- attr(obs, "checks")
  dd  <- deff_clust(est$diagnostics$eif_faw, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow, deff = dd$deff_clust,
                    icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk (parallel-safe L'Ecuyer streams) --------------------------
# RNG stream: iseed = SAMPLE_SEED_BASE + 1000L*cell_index + task.
#   cell_index in 1..16 -> offsets 1000..16000 are >= 1000 apart while
#   task <= 160, so no two (cell, chunk) tasks can ever share a stream (and the
#   offset also keeps this run's streams disjoint from the older runs'
#   SAMPLE_SEED_BASE + task convention). The stream only drives fold assignment
#   (make_cf_folds) + SuperLearner CV splits -- the sample draw itself is
#   .scoped_seed-ed by sample_seed, so estimates stay seed-comparable anyway.
# All worker inputs are passed as EXPLICIT parLapply arguments (pop, learners,
# q_spec, g_spec, oversample), not free closure variables -- this is the R04
# lesson (workers once lacked `pop`); nothing here relies on clusterExport.
iseed <- SAMPLE_SEED_BASE + 1000L * cell_index + task
if (cores > 1L && length(reps) > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    .this_dir <- Sys.getenv("R14_DIR",
      file.path(REPO_ROOT, "simulation", "enhancements", "R14_dr_factorial"))
    source(file.path(.this_dir, "dr_helpers.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner) })
  })
  clusterSetRNGStream(cl, iseed = iseed)
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners,
                        q_spec = q_spec, g_spec = g_spec, oversample = oversample)
} else {
  set.seed(iseed)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = learners,
                     q_spec = q_spec, g_spec = g_spec, oversample = oversample)
}

reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps in this chunk failed",
                                    RUN_ID, task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# tag both tables with the cell knobs (handy for per-rep analysis downstream)
res$q_spec   <- q_spec;   diags$q_spec   <- q_spec
res$g_spec   <- g_spec;   diags$g_spec   <- g_spec
res$sampling <- sampling; diags$sampling <- sampling
res$rung     <- rung;     diags$rung     <- rung

# ---- per-task RDS (full detail; private subtree, never clobbers locked) -----
out <- list(run_id = RUN_ID, scenario = "standard", model_type = "simple",
            q_spec = q_spec, g_spec = g_spec, sampling = sampling, rung = rung,
            oversample = oversample, learners = learners, chunk = chunk,
            reps = reps, Psi = Psi, truth = pop$truth, params = pop$params,
            knobs = list(g_floor = 0.05, Q_clamp = 1e-3, V_cf = 5L),
            results = res, diagnostics = diags)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ---------------
mdir <- file.path(OUT, "manifest")
if (!dir.exists(mdir)) dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
manifest <- list(
  run_id = RUN_ID, scenario = "standard", model_type = "simple",
  q_spec = q_spec, g_spec = g_spec, sampling = sampling, rung = rung,
  oversample = oversample, learners = learners, chunk = chunk, reps = reps,
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),
  knobs = list(g_floor = 0.05, Q_clamp = 1e-3, V_cf = 5L),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED, iseed = iseed, cell_index = cell_index),
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  smoke = SMOKE
)
mf <- file.path(mdir, sprintf("%smanifest_q%s_g%s_%s_%s_chunk%03d.rds",
                              tag, q_spec, g_spec, sampling, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[%s task %d] manifest %s\n", RUN_ID, task, mf))

# ---- SMOKE: inline summary + explicit gate verdict ---------------------------
if (SMOKE) {
  z_or_t <- function(df) qt(0.975, pmax(1, df))
  summ <- do.call(rbind, by(res, res$method, function(d) {
    ok <- is.finite(d$b) & is.finite(d$se)        # drop diverged reps
    bo <- d$b[ok]; so <- d$se[ok]; crit <- z_or_t(d$df[ok])
    cov <- mean(abs(bo - Psi) <= crit * so)
    data.frame(method = d$method[1], n_reps = sum(ok), n_diverged = sum(!ok),
               Psi = Psi, bias = mean(bo) - Psi, emp_sd = sd(bo),
               mean_se = mean(so), coverage = cov, stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL
  summ <- summ[order(match(summ$method, c("FA-w", "CF-u", "CF-w"))), ]
  np <- sapply(summ, is.numeric); sp <- summ; sp[np] <- round(sp[np], 4)
  cat(sprintf("\n[SMOKE] per-arm summary (qC/gC/info/%s; both-correct cell):\n", rung))
  print(sp, row.names = FALSE)

  vfe <- unique(stats::na.omit(diags$V_eff))
  cat(sprintf("[SMOKE] CF folds V_eff = %s (expect 5: standard design has 6 PSUs/stratum, V_cf=5)\n",
              paste(vfe, collapse = ",")))

  # ---- gate: (a) all three arms present with FINITE summary stats; the
  #            UNWEIGHTED CF-u arm must have n_diverged = 0, and the weighted
  #            arms (FA-w, CF-w) a convergent share >= 70% -- weighted GLM
  #            nuisance fits are KNOWN to diverge/collapse under these
  #            heavy-tailed weights (locked R03 standard/L1: SF-W 172/1000,
  #            CF-W 166/1000 diverged; base-R glm itself fails to converge on
  #            those draws), so zero-divergence would be a stricter bar than
  #            the validated engine behavior;
  #            (b) both-correct bias small: |bias| <= 2*emp_sd/sqrt(n_reps);
  #            (c) folds sane: a single V_eff in [2,5].
  reasons <- character(0)
  if (!all(is.finite(summ$bias)) || !all(is.finite(summ$mean_se)) ||
      !all(is.finite(summ$coverage)) || nrow(summ) != 3L)
    reasons <- c(reasons, "non-finite arm summary")
  cfu_div <- summ$n_diverged[summ$method == "CF-u"]
  if (length(cfu_div) != 1L || cfu_div != 0L)
    reasons <- c(reasons, "CF-u (unweighted, paper default) has diverged reps -- real bug")
  conv_share <- summ$n_reps / pmax(1, summ$n_reps + summ$n_diverged)
  if (!all(conv_share >= 0.7))
    reasons <- c(reasons, sprintf("convergent share < 70%% [%s]",
      paste(sprintf("%s %.2f", summ$method, conv_share)[conv_share < 0.7], collapse = "; ")))
  bias_lim <- 2 * summ$emp_sd / sqrt(pmax(1, summ$n_reps))
  if (!all(abs(summ$bias) <= bias_lim))
    reasons <- c(reasons, sprintf("both-correct bias exceeds 2*emp_sd/sqrt(n) [%s]",
      paste(sprintf("%s %.4f>%.4f", summ$method, abs(summ$bias), bias_lim)[abs(summ$bias) > bias_lim],
            collapse = "; ")))
  if (length(vfe) != 1L || vfe < 2L || vfe > 5L)
    reasons <- c(reasons, sprintf("V_eff not sane (%s)", paste(vfe, collapse = ",")))
  if (length(reasons)) {
    cat(sprintf("[SMOKE-GATE] STOP: %s\n", paste(reasons, collapse = " | ")))
  } else {
    cat(sprintf("[SMOKE-GATE] PASS: all 3 arms finite (CF-u clean; weighted-arm divergence within the R03-documented rate), both-correct |bias|<=2*emp_sd/sqrt(n), V_eff=%d\n",
                vfe))
  }
}
