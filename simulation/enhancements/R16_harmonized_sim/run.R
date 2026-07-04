# =====================================================================
# run.R  —  R16_harmonized_sim driver  (NON-DESTRUCTIVE ARC enhancement run)
#
# SPEC: Writing/comments/phase4-arc-sim-specs.md item A9a (harmonized-truncation
# headline simulation). Re-runs the headline 2x4 (scenario x rung) grid with the
# SAME propensity floor (g_floor, default 0.05) and Q-truncation (q_lo, default
# 1e-3) applied to ALL FIVE paper arms BEFORE targeting -- removing the
# floor/cross-fitting confound from Figure 1 / Table 1. Pattern: pre-fit
# nuisances at a harmonized floor -> tmle(Q=, g1W=) targeting, validated by
# R03_isolation_2x2 (its 2x2 isolation cells extended here to the five PAPER
# arms; see estimators_harmonized.R).
#
# GRID: 2 scenarios x 4 rungs = 8 cells; FULL = SIM_N_REPS=1000 reps in
# SIM_CHUNK=100 chunks -> n_chunks=10 -> FULL --array=1-80. Cell-major order:
#   tasks  1-10 standard/L1_param,   11-20 R1/L1_param,
#   tasks 21-30 standard/L2_smooth,  31-40 R1/L2_smooth,
#   tasks 41-50 standard/L3_adaptive,51-60 R1/L3_adaptive,
#   tasks 61-70 standard/L4_aggressive, 71-80 R1/L4_aggressive.
#
# COMPARABILITY: the population is built ONCE per task from POP_SEED and reps
# draw via draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i) with i the
# GLOBAL rep index -- so for a given scenario the samples are BYTE-IDENTICAL to
# the locked headline run (run_sim.R). Coverage differences vs
# results/sim_full_summary.csv are therefore attributable to the harmonized
# truncation alone, not to sampling noise.
#
# NON-DESTRUCTIVE: sources the canonical engine read-only; writes per-task RDS
# only under R16_OUT (+ results/arc via aggregate.R). No codes/*.R is edited.
#
# Env vars (local-test fallbacks in parens):
#   SLURM_ARRAY_TASK_ID  1-based task index            (1)
#   SLURM_CPUS_PER_TASK  cores for this task           (detectCores-1)
#   SIM_N_REPS           total reps per cell           (1000)
#   SIM_CHUNK            reps per array task           (100)
#   SIM_CODE             canonical code dir            ("codes")
#   R16_DIR              this run's folder             (REPO_ROOT/codes/arc_runs/R16_harmonized_sim)
#   R16_OUT              per-task RDS dir              (DATA_ROOT/arc_runs/R16_harmonized_sim)
#   R16_GFLOOR           harmonized propensity floor   (0.05)
#   R16_QLO              harmonized Q truncation       (1e-3)
#   SMOKE=1              tiny validation: standard x L1_param, 20 reps, 2 cores
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))      # exposes .sl/.eif_from_tmle/.se_des/make_cf_folds
source(file.path(CODE, "diagnostics.R"))     # deff_clust, audit_population
source(file.path(CODE, "learners.R"))        # SL_LADDER
# this run's helper (lives in THIS folder; reuses the engine building blocks):
HERE <- Sys.getenv("R16_DIR",
                   file.path(REPO_ROOT, "simulation", "enhancements", "R16_harmonized_sim"))
source(file.path(HERE, "estimators_harmonized.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getn <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.numeric(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
GFLOOR <- getn("R16_GFLOOR", 0.05)
QLO    <- getn("R16_QLO", 1e-3)

OUT <- Sys.getenv("R16_OUT", file.path(DATA_ROOT, "arc_runs", "R16_harmonized_sim"))
MAN <- file.path(OUT, "manifest")
for (d in c(OUT, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x rung (the full headline grid), x rep-chunks -------
RUNGS <- names(SL_LADDER)                            # L1_param..L4_aggressive
cells <- expand.grid(scenario = c("standard", "R1"), rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

if (SMOKE) {
  # tiny validation: ONE cheap cell (standard x L1_param), 20 reps, 2 cores.
  cells  <- cells[cells$scenario == "standard" & cells$rung == "L1_param", , drop = FALSE]
  N_REPS <- 20L; CHUNK <- 20L; cores <- min(cores, 2L); task <- 1L
  cat(sprintf("[SMOKE] standard x L1_param, 20 reps, %d core(s)\n", cores))
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))            # FULL --array=1-80 (see header)
job <- grid[task, ]
scenario <- job$scenario; rung <- job$rung; chunk <- job$chunk
cell_index <- job$cell
learners <- SL_LADDER[[rung]]
# GLOBAL rep index so chunks never collide and sample seeds match run_sim.R exactly.
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[R16 task %d] scenario=%s rung=%s chunk=%d reps=%d-%d g_floor=%.3f q_lo=%g cores=%d%s\n",
            task, scenario, rung, chunk, rep_lo, rep_hi, GFLOOR, QLO, cores,
            if (SMOKE) " [SMOKE]" else ""))

# ---- per-chunk checkpoint: compute fn EARLY; skip instantly if done ----------
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%sr16_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic from POP_SEED) -------------
pop <- make_population(scenario, model_type = "complex",
                       truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[R16 task %d] Psi=%.5f (mc se %.6f) Psi_N=%.5f gap=%.6f\n",
            task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0,
            abs(pop$truth$gap_super_census)))

one_rep <- function(i, pop, learners, g_floor, q_lo) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est <- run_harmonized(obs, learners = learners, g_floor = g_floor, q_lo = q_lo,
                        V_cf = 5L, inner_cv_folds = 5L)
  ch  <- attr(obs, "checks")
  # DEFF/ICC on the Fully-Aware-h arm's EIF (analogue of run_sim.R's eif_fa audit)
  dd  <- deff_clust(est$diagnostics$eif_fah, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, est$results, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow, deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk's reps (parallel over reps, mirrors run_sim.R / R03) ------
# RNG stream offset: iseed = SAMPLE_SEED_BASE + 1000L*cell_index + task. The
# 1000L*cell stride dominates the task index (max 80 < 1000), so no two
# (cell, chunk) combinations ever share a stream. (Sample REALIZATIONS are pinned
# by draw_sample's per-rep sample_seed regardless; the stream only governs
# SL-internal CV splits / ranger randomness.)
if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    HERE <- Sys.getenv("R16_DIR",
                       file.path(REPO_ROOT, "simulation", "enhancements", "R16_harmonized_sim"))
    source(file.path(HERE, "estimators_harmonized.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger) })
  })
  clusterSetRNGStream(cl, iseed = SAMPLE_SEED_BASE + 1000L * cell_index + task)
  # belt-and-braces export (R04 lesson: workers once lacked `pop`); the same
  # objects are ALSO passed as explicit parLapply args below.
  clusterExport(cl, c("pop", "learners", "GFLOOR", "QLO", "one_rep"),
                envir = environment())
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = learners,
                        g_floor = GFLOOR, q_lo = QLO)
} else {
  set.seed(SAMPLE_SEED_BASE + 1000L * cell_index + task)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = learners,
                     g_floor = GFLOOR, q_lo = QLO)
}
# account for hard-errored reps (one_rep tryCatch -> NULL): unlike the diverged-
# guard NAs (which keep their rows and are counted as n_diverged), these reps
# vanish from the chunk entirely -- record n_err so the loss is never silent.
n_err <- sum(vapply(reps_out, is.null, logical(1)))
reps_out <- Filter(Negate(is.null), reps_out)   # drop reps that errored
if (n_err > 0L)
  cat(sprintf("[R16 task %d] WARNING: %d/%d reps errored (tryCatch -> NULL) and are MISSING from this chunk's rows (see worker 'rep %%d failed' messages in the .err log)\n",
              task, n_err, length(reps)))
if (!length(reps_out)) stop(sprintf("[R16 task %d] all %d reps in this chunk failed",
                                    task, length(reps)))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (full detail; own folder; never clobbers locked files) --
out <- list(run = "R16_harmonized_sim", scenario = scenario, rung = rung,
            learners = learners, g_floor = GFLOOR, q_lo = QLO,
            chunk = chunk, reps = reps, n_err = n_err,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags, smoke = SMOKE)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[R16 task %d] saved %s  (%d est rows, %d diag rows)\n",
            task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----------------
manifest <- list(
  run = "R16_harmonized_sim", scenario = scenario, rung = rung, learners = learners,
  g_floor = GFLOOR, q_lo = QLO, chunk = chunk, reps = reps, n_err = n_err,
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
mf <- file.path(MAN, sprintf("%smanifest_r16_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
saveRDS(manifest, mf)
cat(sprintf("[R16 task %d] manifest %s\n", task, mf))

# ---- SMOKE: inline summary + explicit gate verdict ----------------------------
# (full runs are aggregated separately by aggregate.R; see NOTES.md.)
if (SMOKE) {
  z_or_t <- function(df) qt(0.975, pmax(1, df))
  agg <- do.call(rbind, by(res, res$method, function(d) {
    ok <- is.finite(d$b) & is.finite(d$se)
    bo <- d$b[ok]; so <- d$se[ok]; crit <- z_or_t(d$df[ok])
    cov <- if (sum(ok)) mean(abs(bo - Psi) <= crit * so) else NA_real_
    data.frame(method = d$method[1], n_reps = sum(ok), n_diverged = sum(!ok),
               bias = mean(bo) - Psi, emp_sd = sd(bo), mean_se = mean(so),
               se_ratio = mean(so) / sd(bo),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / max(1, sum(ok))))
  }))
  rownames(agg) <- NULL
  ord <- c("Fully-Aware-h", "Fully-Aware-CF-h", "Fully-Aware-CV-h",
           "Partially-Aware-h", "Non-Aware-h")
  agg <- agg[order(match(agg$method, ord)), ]
  cat(sprintf("\n==== R16 SMOKE summary (%s / %s, g_floor=%.3f, q_lo=%g, %d reps) ====\n",
              scenario, rung, GFLOOR, QLO, length(unique(res$rep))))
  num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
  print(ap, row.names = FALSE)

  # ---- the gate: arms 1-3,5 present (no CV at single-learner rung); all b/se
  # finite EXCEPT a bounded share of guard-flagged non-convergent reps in the
  # WEIGHTED-nuisance arms; FA-h vs PA-h psi IDENTICAL per rep (same fit; only
  # the SE differs).
  # NOTE on "all finite": the spec's original gate assumed the GLM rung is clean,
  # but the locked R03_isolation_2x2 FULL run already showed the weighted
  # single-fit-nuisance + pooled-targeting path diverging at L1_param (SF-W:
  # 172/1000 reps, results/arc/R03_isolation_2x2_summary.csv); FA-h/PA-h are that
  # same path, so the gate requires (a) the UNWEIGHTED-nuisance arms (Non-Aware-h,
  # Fully-Aware-CF-h) strictly all-finite, (b) FA-h divergence bounded (<= 0.25,
  # vs R03's 0.172 reference) and (c) every NA exactly matching a diverged flag.
  need <- c("Fully-Aware-h", "Partially-Aware-h", "Non-Aware-h", "Fully-Aware-CF-h")
  have <- unique(res$method)
  missing_arms <- setdiff(need, have)
  cv_expected  <- length(learners) > 1
  cv_present   <- "Fully-Aware-CV-h" %in% have
  arm_finite <- function(m) { d <- res[res$method == m, ]
                              nrow(d) > 0 && all(is.finite(d$b) & is.finite(d$se)) }
  unw_finite <- arm_finite("Non-Aware-h") && arm_finite("Fully-Aware-CF-h")
  guard_consistent <- all(res$diverged == !(is.finite(res$b) & is.finite(res$se)))
  div_frac_fa <- mean(res$diverged[res$method == "Fully-Aware-h"])
  mfa <- res[res$method == "Fully-Aware-h",     c("rep", "b")]
  mpa <- res[res$method == "Partially-Aware-h", c("rep", "b")]
  mm  <- merge(mfa, mpa, by = "rep", suffixes = c("_fa", "_pa"))
  psi_ident <- nrow(mm) > 0 &&
    all((mm$b_fa == mm$b_pa) | (is.na(mm$b_fa) & is.na(mm$b_pa)))
  stopifnot(psi_ident)   # hard invariant: arms 1+2 share ONE fit

  fails <- character(0)
  if (length(missing_arms))
    fails <- c(fails, paste("missing arms:", paste(missing_arms, collapse = ",")))
  if (!cv_expected && cv_present)
    fails <- c(fails, "CV arm present at a single-learner rung")
  if (cv_expected && !cv_present)
    fails <- c(fails, "CV arm missing at a multi-learner rung")
  if (!unw_finite)
    fails <- c(fails, "non-finite b/se in an UNWEIGHTED-nuisance arm (NA-h/CF-h must never diverge; R03: 1000/1000 finite)")
  if (!guard_consistent)
    fails <- c(fails, "divergence guard inconsistent (NA b/se without diverged flag, or flagged rep with finite b/se)")
  if (is.finite(div_frac_fa) && div_frac_fa > 0.25)
    fails <- c(fails, sprintf("FA-h diverged fraction %.2f > 0.25 (R03 L1 SF-W reference: 0.172)", div_frac_fa))
  if (!psi_ident)
    fails <- c(fails, "FA-h vs PA-h psi differ (must share one fit)")

  if (length(fails)) {
    cat(sprintf("\n[SMOKE-GATE] STOP: %s\n", paste(fails, collapse = "; ")))
  } else {
    cat(sprintf(paste0(
      "\n[SMOKE-GATE] PASS: arms 1-3,5 present (no CV at the single-learner rung); ",
      "unweighted arms all finite; FA-h guarded divergence %.2f within the R03-known ",
      "band (<= 0.25); FA-h == PA-h psi identical per rep. Submit the full --array=1-80.\n"),
      div_frac_fa))
  }
}
