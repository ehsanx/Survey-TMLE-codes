# =====================================================================
# run.R  —  R21_deployable_cvcf driver  (NON-DESTRUCTIVE ARC run)
#
# SPEC: Writing/JASA/R21_L6_ARC_SPEC.md  (DESIGN ONLY; >10 min/cell -> ARC).
#
# GOAL. Two gaps closed in one run, at a DEPLOYABLE non-Donsker library L6:
#   1) The lone-forest optic. L5 (R20) = c(SL.glm, SL.earth, SL.ranger.deep) is
#      a 3-member library whose one aggressive member is still a near-lone
#      interpolating forest. L6 is the library one would ACTUALLY deploy: a
#      richer CV-SuperLearner with TUNED xgboost carrying the aggressive load
#      and the deep RF demoted to one of 5 members.
#        L6 = c('SL.glm','SL.earth','SL.glmnet','SL.xgboost.tuned','SL.ranger.deep')
#      (optional 'SL.hal9001' appended behind R21_WITH_HAL=1; OFF by default --
#       slow / can fail to install on some ARC images; never make the headline
#       depend on it).
#   2) The missing AIPW-CV cell. The AIPW benchmark (R15) has only AIPW-SF and
#      AIPW-CF; the TMLE engine has a cluster-aware Fully-Aware-CV arm but its
#      AIPW counterpart was never built. R21 adds AIPW-CV (r21_helpers.R) so the
#      "cluster-CV != CF" question is answered ESTIMATOR-AGNOSTICALLY on a
#      deployable class.
#
# ARM MATRIX — 6 HEADLINE arms per cell (+2 inherited TMLE arms):
#   TMLE-{single-fit, cluster-CV, CF}   = Fully-Aware / Fully-Aware-CV /
#                                          Fully-Aware-CF  (emitted by
#                                          run_estimators when length(L6)>1)
#   AIPW-{single-fit, cluster-CV, CF}   = AIPW-SF / AIPW-CV / AIPW-CF
#                                          (aipw_arms_all in r21_helpers.R)
#   + Partially-Aware + Non-Aware       = inherited UNCHANGED from run_estimators
#                                          (i.i.d.-variance and unweighted arms;
#                                          the variance/weight-naivety axes, kept
#                                          in the CSV as a deployable-L6 counterpart
#                                          to the Figure-1 set, but NOT in the
#                                          headline 3x2 DECISION VIEW).
#
# NON-DESTRUCTIVE: SOURCES the canonical engine + R15's aipw_helpers.R
# read-only. L6 + SL.xgboost.tuned + the AIPW-CV arm are defined HERE (in
# r21_helpers.R), so NO engine file is edited. run_estimators() is UNCHANGED;
# with a length-5 library it emits FIVE TMLE arms: the three fold-protocol arms
# above PLUS Partially-Aware and Non-Aware. All five are kept; only the AIPW side
# is limited to the three fold-protocol arms (no AIPW variance/weight-naivety
# analogue is built or needed -- that axis is settled by the main sim + R03).
#
# GRID: scenario in {standard, R1} x library L6 = 2 cells, each x
#   ceil(SIM_N_REPS/SIM_CHUNK) chunks.
#   FULL: SIM_N_REPS=1000, SIM_CHUNK=100 -> 2 x 10 = 20 tasks -> --array=1-20
#     tasks  1-10  standard L6   |   tasks 11-20  R1 L6
#   OPTIONAL m-sweep behind R21_MSWEEP=1 (base_m {6,20}) -> doubles to 40 tasks
#     (scenario-major within each base_m block; see grid build below).
#
# SEEDING (byte-identical to locked R15/R20): per-rep sample seed
#   draw_sample(sample_seed = SAMPLE_SEED_BASE + i); parallel stream
#   iseed = SAMPLE_SEED_BASE + 1000*cell_index + task. The samples self-seed
#   from the GLOBAL rep index i, so for each (scenario, rep) the analysis sample
#   is byte-identical to the locked headline / R15 / R20 runs; only the library
#   and the added AIPW-CV arm differ.
#
# Env vars (local fallbacks in brackets):
#   SLURM_ARRAY_TASK_ID  1-based task index                  [1]
#   SLURM_CPUS_PER_TASK  cores for this task                 [detectCores-1]
#   SIM_N_REPS           total reps per cell                 [1000]
#   SIM_CHUNK            reps per array task                 [100]
#   SIM_CODE             canonical code dir                  [codes]
#   R21_DIR              this run folder    [REPO_ROOT/codes/arc_runs/R21_deployable_cvcf]
#   R21_OUT              per-task RDS dir   [DATA_ROOT/arc_runs/R21_deployable_cvcf]
#   R21_WITH_HAL         "1" -> append SL.hal9001 to L6 (OFF by default)
#   R21_MSWEEP           "1" -> base_m sweep {6,20} (doubles the array)
#   R21_XGB_TUNE         "1" -> genuinely CV-tune xgboost depth {3,6,8}
#                               (default: fixed max_depth=6; see r21_helpers.R)
#   R21_XGB_MAXDEPTH     fixed xgboost depth when not tuning  [6]
#   SMOKE                "1" -> tiny plumbing cell (standard, cheap
#                               glm+earth+glmnet, 8 reps; asserts ALL 6 arms,
#                               incl. BOTH *-CV arms, emit finite b/se)
# =====================================================================

WITH_HAL <- identical(Sys.getenv("R21_WITH_HAL"), "1")   # run HAL as a DISTINCT comparable library
RUN_ID   <- "R21_deployable_cvcf"
RUNG_LAB <- if (WITH_HAL) "L7_hal" else "L6_deployable"
# distinct output / manifest filenames so an L6 run and an L6+HAL run never collide and the
# per-chunk checkpoint never cross-skips; aggregate.R reads BOTH families into one summary.
RBASE  <- if (WITH_HAL) "r21hal"      else "r21"
MFBASE <- if (WITH_HAL) "halmanifest" else "manifest"

# ---- locate + source the canonical engine (read-only) ----------------
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) {
  cand  <- c("R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))      # POP_SEED, SAMPLE_SEED_BASE, paths, RNGkind
source(file.path(CODE, "dgp.R"))         # make_population, draw_sample
source(file.path(CODE, "estimators.R"))  # run_estimators, .sl, .se_des, make_cf_folds
source(file.path(CODE, "diagnostics.R")) # deff_clust, audit_population
source(file.path(CODE, "learners.R"))    # SL.ranger.deep wrapper + SL_LADDER
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

.this_dir <- Sys.getenv("R21_DIR", file.path(REPO_ROOT, "simulation", "enhancements", RUN_ID))
# R15's AIPW machinery (aipw_arms = SF+CF, .aipw_jkn, .AIPW_Q_CLAMP, one_rep_aipw)
R15_DIR <- Sys.getenv("R15_DIR", file.path(REPO_ROOT, "simulation", "enhancements", "R15_aipw_benchmark"))
source(file.path(R15_DIR, "aipw_helpers.R"))
# R21-local: SL.xgboost.tuned wrapper + AIPW-CV arm + aipw_arms_all (SF+CF+CV)
source(file.path(.this_dir, "r21_helpers.R"))
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# ---- the DEPLOYABLE non-Donsker MULTI-learner library L6 (defined here) ----
L6 <- c("SL.glm", "SL.earth", "SL.glmnet", "SL.xgboost.tuned", "SL.ranger.deep")
if (WITH_HAL) L6 <- c(L6, "SL.hal9001")   # L7_hal = L6 + Highly Adaptive Lasso (van der Laan 2017)

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE   <- identical(Sys.getenv("SMOKE"), "1")
MSWEEP  <- identical(Sys.getenv("R21_MSWEEP"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
G_FLOOR  <- 0.05      # harmonized propensity floor (= engine OOF floor; R15)
V_CF     <- 5L        # cross-fit folds (engine default)
INNER_CV <- 5L        # internal-CV folds for the CV arms (engine default)

OUT <- Sys.getenv("R21_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
MAN_DIR <- file.path(OUT, "manifest")
for (d in c(OUT, MAN_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x library(L6) [x base_m if MSWEEP] -----------
# core grid: scenario-major, 1 rung (L6). base_m defaults are scenario-specific
# inside draw_sample (standard 6, R1 2); MSWEEP overrides with an explicit
# {6,20} block (scenario fastest within each base_m).
LEARNERS <- L6
if (MSWEEP) {
  cells <- expand.grid(scenario = c("standard", "R1"), base_m = c(6L, 20L),
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
} else {
  cells <- data.frame(scenario = c("standard", "R1"), base_m = NA_integer_,
                      stringsAsFactors = FALSE)
}

if (SMOKE) {
  # cheap multi-learner library so the CV gate still fires (length>1) and all
  # six arms emit, but it finishes fast. The deep-RF / xgboost science path is
  # validated by R15/R20; this is a plumbing gate only.
  cells    <- data.frame(scenario = "standard", base_m = NA_integer_,
                         stringsAsFactors = FALSE)
  LEARNERS <- c("SL.glm", "SL.earth", "SL.glmnet")
  MSWEEP   <- FALSE
  N_REPS <- 8L; CHUNK <- 8L; cores <- min(cores, 2L); task <- 1L
  cat(sprintf("[SMOKE] %s: standard x {%s}, %d reps, %d cores\n",
              RUN_ID, paste(LEARNERS, collapse = "+"), N_REPS, cores))
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))   # FULL = 2 cells x 10 chunks = 20 tasks (40 if MSWEEP)
job        <- grid[task, ]
scenario   <- job$scenario
base_m     <- job$base_m                    # NA -> draw_sample uses scenario default
cell_index <- job$cell_index
chunk      <- job$chunk
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

mlab <- if (is.na(base_m)) "default" else as.character(base_m)
cat(sprintf("[%s task %d] scenario=%s rung=%s base_m=%s chunk=%d reps=%d-%d cores=%d learners=%s%s\n",
            RUN_ID, task, scenario, RUNG_LAB, mlab, chunk, rep_lo, rep_hi, cores,
            paste(LEARNERS, collapse = "+"), if (SMOKE) "  [SMOKE]" else ""))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
mtag <- if (is.na(base_m)) "" else sprintf("_m%02d", base_m)
fn <- file.path(OUT, sprintf("%s%s_%s%s_chunk%03d.rds", tag, RBASE, scenario, mtag, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic; headline te_log_odds=1.5) ----
pop <- make_population(scenario, model_type = "complex",
                       pop_seed = POP_SEED, truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.6f (mc se %.6f) Psi_N=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0))

# ---- one replication: all SIX arms at L6 ------------------------------
# i is the GLOBAL rep index -> sample_seed = SAMPLE_SEED_BASE + i, so for a
# given (scenario, rep) the sampled units/C/A/weights/design are BYTE-IDENTICAL
# to the locked headline run; only the LIBRARY + the AIPW-CV arm differ.
# TMLE side: run_estimators(length(L6)>1) emits Fully-Aware / -CV / -CF PLUS the
#   inherited Partially-Aware + Non-Aware (all five kept; the latter two headline
#   the variance/weight-naivety axes, not the CV-vs-CF axis R21 targets).
# AIPW side: aipw_arms_all() emits AIPW-SF / AIPW-CV / AIPW-CF on the SAME obs.
one_rep <- function(i, pop, learners, base_m) tryCatch({
  obs <- if (is.na(base_m))
    draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  else
    draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex",
                base_m = base_m)

  # --- TMLE arms (engine, unchanged): Fully-Aware / -CV / -CF ---
  est <- run_estimators(obs, learners = learners,
                        V_cf = V_CF, inner_cv_folds = INNER_CV, g_oof_bound = G_FLOOR)
  tmle_res <- est$results
  tmle_res$estimator <- "TMLE"

  # --- AIPW arms (R15 + R21-local): AIPW-SF / AIPW-CV / AIPW-CF ---
  ai <- aipw_arms_all(obs, learners = learners, g_floor = G_FLOOR, V_cf = V_CF,
                      inner_cv_folds = INNER_CV)
  # carry only the linearized (Eq-8) SE into the unified `se` column so the
  # TMLE/AIPW DECISION VIEW is apples-to-apples (both EIF-based); the JKn SE is
  # kept in the AIPW detail table (ai$results) for the appendix.
  aipw_res <- data.frame(method = ai$results$method, b = ai$results$b,
                         se = ai$results$se_lin, df = ai$results$df,
                         se_jkn = ai$results$se_jkn, estimator = "AIPW",
                         stringsAsFactors = FALSE)
  # align columns: TMLE rows have no se_jkn -> NA
  tmle_res$se_jkn <- NA_real_
  res <- rbind(tmle_res[, c("method", "b", "se", "df", "se_jkn", "estimator")],
               aipw_res[, c("method", "b", "se", "df", "se_jkn", "estimator")])

  ch <- attr(obs, "checks")
  dd <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, res, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, est$diagnostics$drow,
                    ai$diagnostics$drow, ai$diagnostics$drow_cv,
                    deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk across cores (parallel-safe L'Ecuyer streams) -------
# Stream seed SAMPLE_SEED_BASE + 1000*cell_index + task: unique WITHIN this run
# and disjoint from older runs' task range. Samples self-seed in draw_sample
# from the global rep index, so streams drive only learner-internal RNG/folds.
iseed <- SAMPLE_SEED_BASE + 1000L * cell_index + task
if (cores > 1L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    if (!file.exists(file.path(CODE, "config.R"))) CODE <- "R"
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R"))
    R15_DIR <- Sys.getenv("R15_DIR", file.path(REPO_ROOT, "simulation", "enhancements", "R15_aipw_benchmark"))
    source(file.path(R15_DIR, "aipw_helpers.R"))
    R21_DIR <- Sys.getenv("R21_DIR", file.path(REPO_ROOT, "simulation", "enhancements", "R21_deployable_cvcf"))
    source(file.path(R21_DIR, "r21_helpers.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner)
                       library(ranger); library(earth); library(glmnet)
                       library(xgboost); library(surveyCV) })
  })
  clusterSetRNGStream(cl, iseed = iseed)
  clusterExport(cl, c("pop", "scenario", "V_CF", "INNER_CV", "G_FLOOR"), envir = environment())
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = LEARNERS, base_m = base_m)
} else {
  set.seed(iseed)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = LEARNERS, base_m = base_m)
}
reps_out <- Filter(Negate(is.null), reps_out)
n_failed <- length(reps) - length(reps_out)
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps in this chunk failed",
                                    RUN_ID, task, length(reps)))
if (n_failed > 0)
  cat(sprintf("[%s task %d] WARNING: %d/%d reps failed; chunk RDS will checkpoint as complete -- delete %s to recompute\n",
              RUN_ID, task, n_failed, length(reps), fn))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
res$rung   <- RUNG_LAB
res$base_m <- if (is.na(base_m)) NA_integer_ else base_m
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (full detail; private folder) ---------------------
out <- list(run_id = RUN_ID, scenario = scenario, rung = RUNG_LAB, base_m = base_m,
            learners = LEARNERS, chunk = chunk, reps = reps, n_failed = n_failed,
            g_floor = G_FLOOR, V_cf = V_CF, inner_cv = INNER_CV,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags)
saveRDS(out, fn)
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----------
manifest <- list(
  run_id = RUN_ID, scenario = scenario, rung = RUNG_LAB, base_m = base_m,
  learners = LEARNERS, chunk = chunk, reps = reps, n_failed = n_failed,
  g_floor = G_FLOOR, V_cf = V_CF, inner_cv = INNER_CV,
  with_hal = identical(Sys.getenv("R21_WITH_HAL"), "1"),
  xgb_tune = identical(Sys.getenv("R21_XGB_TUNE"), "1"),
  xgb_maxdepth = geti("R21_XGB_MAXDEPTH", 6L),
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED, iseed_stream = iseed),
  run_dir = .this_dir,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam",
                      "glmnet","ranger","xgboost","hal9001"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  smoke = SMOKE
)
mf <- file.path(MAN_DIR, sprintf("%s%s_%s%s_chunk%03d.rds", tag, MFBASE, scenario, mtag, chunk))
saveRDS(manifest, mf)
cat(sprintf("[%s task %d] manifest %s\n", RUN_ID, task, mf))

# ---- SMOKE: inline summary + explicit gate verdict ----------------------
# THE key check: ALL SIX arms present, incl. BOTH *-CV arms, all finite b/se.
if (SMOKE) {
  crit <- qt(0.975, pmax(1, res$df))
  res$cover <- abs(res$b - Psi) <= crit * res$se
  summ <- do.call(rbind, by(res, res$method, function(d) data.frame(
    estimator = d$estimator[1], method = d$method[1], n_reps = nrow(d),
    bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
    coverage = mean(d$cover), stringsAsFactors = FALSE)))
  rownames(summ) <- NULL
  ord_m <- c("Fully-Aware", "Fully-Aware-CV", "Fully-Aware-CF",
             "AIPW-SF", "AIPW-CV", "AIPW-CF")
  summ <- summ[order(match(summ$method, ord_m)), ]
  np <- sapply(summ, is.numeric); sp <- summ; sp[np] <- round(sp[np], 4)
  cat(sprintf("\n[SMOKE] inline summary (standard %s; true Psi=%.6f):\n",
              paste(LEARNERS, collapse = "+"), Psi))
  print(sp, row.names = FALSE)

  finite_ok <- all(is.finite(res$b)) && all(is.finite(res$se))
  meths     <- unique(res$method)
  want      <- c("Fully-Aware", "Fully-Aware-CV", "Fully-Aware-CF",
                 "AIPW-SF", "AIPW-CV", "AIPW-CF")
  # run_estimators ALSO emits Partially-Aware + Non-Aware (estimators.R:168,174),
  # inherited free, kept in the per-task RDS + full summary CSV and excluded only
  # from the headline 3x2 DECISION VIEW. So the gate checks the six REQUIRED arms
  # are PRESENT (subset), NOT that the method set equals six.
  six_ok    <- all(want %in% meths)
  extra     <- setdiff(meths, want)            # expected: Partially-Aware, Non-Aware
  tmle_cv   <- "Fully-Aware-CV" %in% meths
  aipw_cv   <- "AIPW-CV" %in% meths      # the NEW estimator-agnostic cell
  reasons <- paste(c(
    sprintf("b/se finite=%s", finite_ok),
    sprintf("6 required arms present=%s", six_ok),
    sprintf("TMLE-CV present=%s", tmle_cv),
    sprintf("AIPW-CV present=%s (the new cell)", aipw_cv),
    sprintf("inherited extras={%s}", paste(extra, collapse = ","))), collapse = "; ")
  if (finite_ok && six_ok && tmle_cv && aipw_cv)
    cat(sprintf("\n[SMOKE-GATE] PASS: %s\n", reasons))
  else
    cat(sprintf("\n[SMOKE-GATE] STOP: %s\n", reasons))
}
