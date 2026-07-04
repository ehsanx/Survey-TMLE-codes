# =====================================================================
# run.R  —  R22_aipw_ladder_complete driver  (NON-DESTRUCTIVE ARC run)
#
# SPEC: Writing/comments/R22_AIPW_COMPLETE_SPEC.md.  README: README.md.
#
# GOAL (minimal "alpha" completion). Add ONLY the genuinely-missing AIPW cells
# of the unified L1-L7 x {TMLE,AIPW} x {single-fit, cluster-CV, CF} appendix
# figure -- nothing already computed elsewhere is recomputed:
#       AIPW-{SF,CV,CF} @ L5          (R20 ran TMLE only at L5)
#       AIPW-CV         @ L2, L3      (R15 ran AIPW-SF/CF at L2/L3, never the CV arm)
# NO TMLE arms (they exist in sim_full @ L1-L4 and R20 @ L5); NO AIPW-SF/CF at
# L2/L3 (they exist in R15). `L1/L4 internal-CV` is N/A by construction (a lone
# learner has nothing to weight -> internal CV degenerates to single-fit).
#
# So the per-rung arm set is:
#       L2_smooth     -> AIPW-CV                       (aipw_cv_arm)
#       L3_adaptive   -> AIPW-CV                       (aipw_cv_arm)
#       L5_nondonsker -> AIPW-SF, AIPW-CV, AIPW-CF     (aipw_arms_all)
#
# The manuscript continues to SOURCE every locked cell from sim_full/R15/R20/R21;
# R22 adds only the new cells, so NO locked number changes (audit only extends).
# Trust in the new cells comes from CODE IDENTITY: the AIPW machinery is the SAME
# (R15 aipw_arms / R21 aipw_cv_arm), already validated at L1-L4 (R15) and L6/L7
# (R21); R22 only applies it to the L2/L3/L5 libraries on the SAME samples.
#
# LIBRARIES (published / R20 definitions -- never edits SL_LADDER):
#   L2_smooth     = c('SL.glm','SL.gam','SL.earth')              (SL_LADDER)
#   L3_adaptive   = c('SL.glm','SL.earth','SL.ranger.t1')        (SL_LADDER)
#   L5_nondonsker = c('SL.glm','SL.earth','SL.ranger.deep')      (= R20 run.R:67)
#
# NON-DESTRUCTIVE: SOURCES the engine + R15's aipw_helpers.R + R21's r21_helpers.R
#   (aipw_cv_arm / aipw_arms_all) ALL read-only. No xgboost / HAL is used here, so
#   none of R21's HAL install fragility.
#
# GRID: (scenario in {standard, R1}) x (rung in {L2,L3,L5}) = 6 cells,
#   each x ceil(SIM_N_REPS/SIM_CHUNK) chunks.
#   FULL: SIM_N_REPS=1000, SIM_CHUNK=100 -> 6 x 10 = 60 tasks -> --array=1-60
#     cell order (scenario fastest): 1 std-L2 2 R1-L2 3 std-L3 4 R1-L3 5 std-L5 6 R1-L5
#     => tasks  1-10 std-L2 | 11-20 R1-L2 | 21-30 std-L3 | 31-40 R1-L3 |
#               41-50 std-L5 | 51-60 R1-L5    (L5 = the expensive block; probe task 41)
#
# SEEDING (byte-identical samples to locked sim_full/R15/R20/R21): per-rep
#   draw_sample(sample_seed = SAMPLE_SEED_BASE + i); parallel stream
#   iseed = SAMPLE_SEED_BASE + 1000*cell_index + task.
#
# Env vars (local fallbacks in brackets):
#   SLURM_ARRAY_TASK_ID [1] · SLURM_CPUS_PER_TASK [detectCores-1] ·
#   SIM_N_REPS [1000] · SIM_CHUNK [100] · SIM_CODE [codes] ·
#   R22_DIR [REPO_ROOT/codes/arc_runs/R22_aipw_ladder_complete] ·
#   R22_OUT [DATA_ROOT/arc_runs/R22_aipw_ladder_complete] ·
#   R15_DIR [.../R15_aipw_benchmark] · R21_DIR [.../R21_deployable_cvcf] ·
#   SMOKE "1" -> tiny plumbing cell (standard, cheap glm+earth, "all" arms, 8
#                reps; asserts AIPW-SF/CV/CF emit finite b/se)
# =====================================================================

RUN_ID <- "R22_aipw_ladder_complete"

# ---- locate + source the canonical engine (read-only) ----------------
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) {
  cand  <- c("R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))      # POP_SEED, SAMPLE_SEED_BASE, paths, RNGkind
source(file.path(CODE, "dgp.R"))         # make_population, draw_sample (sets attr(obs,'checks'))
source(file.path(CODE, "estimators.R"))  # .se_des, make_cf_folds (used by aipw helpers)
source(file.path(CODE, "diagnostics.R")) # deff_clust, audit_population
source(file.path(CODE, "learners.R"))    # SL_LADDER (L2/L3), SL.ranger.t1/.deep wrappers
source(file.path(CODE, "arc_runs", "_checkpoint.R"))   # arc_skip_if_done(), arc_git_sha()

.this_dir <- Sys.getenv("R22_DIR", file.path(REPO_ROOT, "simulation", "enhancements", RUN_ID))
R15_DIR <- Sys.getenv("R15_DIR", file.path(REPO_ROOT, "simulation", "enhancements", "R15_aipw_benchmark"))
R21_DIR <- Sys.getenv("R21_DIR", file.path(REPO_ROOT, "simulation", "enhancements", "R21_deployable_cvcf"))
source(file.path(R15_DIR, "aipw_helpers.R"))  # aipw_arms (SF+CF), .aipw_jkn, .AIPW_Q_CLAMP
source(file.path(R21_DIR, "r21_helpers.R"))   # aipw_cv_arm + aipw_arms_all (SF+CF+CV)
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# ---- the libraries to complete + the arms wanted at each (minimal scope) ----
LIBS <- list(
  L2_smooth     = SL_LADDER$L2_smooth,                       # c('SL.glm','SL.gam','SL.earth')
  L3_adaptive   = SL_LADDER$L3_adaptive,                     # c('SL.glm','SL.earth','SL.ranger.t1')
  L5_nondonsker = c("SL.glm", "SL.earth", "SL.ranger.deep")  # = R20 run.R:67
)
ARMS <- list(L2_smooth = "cv", L3_adaptive = "cv", L5_nondonsker = "all")  # cv = AIPW-CV only; all = SF+CV+CF

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE  <- identical(Sys.getenv("SMOKE"), "1")

task   <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 1000L)
CHUNK  <- geti("SIM_CHUNK", 100L)
G_FLOOR  <- 0.05      # harmonized propensity floor (= engine OOF floor; R15/R21)
V_CF     <- 5L        # cross-fit folds (engine default)
INNER_CV <- 5L        # internal-CV folds for the CV arms (engine default)

OUT <- Sys.getenv("R22_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
MAN_DIR <- file.path(OUT, "manifest")
for (d in c(OUT, MAN_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: scenario x rung -------------------------------------
cells <- expand.grid(scenario = c("standard", "R1"), rung = names(LIBS),
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

if (SMOKE) {
  cells   <- data.frame(scenario = "standard", rung = "SMOKE_cheap",
                        stringsAsFactors = FALSE)
  LIBS    <- c(LIBS, list(SMOKE_cheap = c("SL.glm", "SL.earth")))
  ARMS    <- c(ARMS, list(SMOKE_cheap = "all"))
  N_REPS <- 8L; CHUNK <- 8L; cores <- min(cores, 2L); task <- 1L
  cat(sprintf("[SMOKE] %s: standard x {%s}, AIPW 'all', %d reps, %d cores\n",
              RUN_ID, paste(LIBS$SMOKE_cheap, collapse = "+"), N_REPS, cores))
}

n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
  cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))   # FULL = 6 cells x 10 chunks = 60 tasks
job        <- grid[task, ]
scenario   <- job$scenario
rung       <- job$rung
cell_index <- job$cell_index
chunk      <- job$chunk
LEARNERS   <- LIBS[[rung]]
WANT       <- ARMS[[rung]]
rep_lo <- (chunk - 1L) * CHUNK + 1L
rep_hi <- min(chunk * CHUNK, N_REPS)
reps   <- rep_lo:rep_hi

cat(sprintf("[%s task %d] scenario=%s rung=%s arms=%s chunk=%d reps=%d-%d cores=%d learners=%s%s\n",
            RUN_ID, task, scenario, rung, WANT, chunk, rep_lo, rep_hi, cores,
            paste(LEARNERS, collapse = "+"), if (SMOKE) "  [SMOKE]" else ""))

# ---- per-chunk checkpoint: skip instantly if this chunk already finished ----
tag <- if (SMOKE) "SMOKE_" else ""
fn <- file.path(OUT, sprintf("%sr22_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- build population + truth ONCE (deterministic; headline te_log_odds=1.5) ----
pop <- make_population(scenario, model_type = "complex",
                       pop_seed = POP_SEED, truth_M = if (SMOKE) 2e5L else 2e6L)
Psi <- pop$truth$psi
cat(sprintf("[%s task %d] Psi=%.6f (mc se %.6f) Psi_N=%.6f\n",
            RUN_ID, task, Psi, pop$truth$se_mc, pop$truth$psi_N_Q0))

# ---- one replication: ONLY the wanted AIPW arms at `learners` ----------
# i is the GLOBAL rep index -> sample_seed = SAMPLE_SEED_BASE + i, so for a given
# (scenario, rep) the sampled units/C/A/weights/design are BYTE-IDENTICAL to the
# locked runs; only the library + the AIPW arms differ.
#   want='all' -> aipw_arms_all (AIPW-SF, AIPW-CF, AIPW-CV)
#   want='cv'  -> aipw_cv_arm   (AIPW-CV only)
# deff_clust / icc_eif are computed from the AIPW cluster-CV influence function
# (present in every cell), not from a TMLE EIF (no TMLE is run here).
one_rep <- function(i, pop, learners, want) tryCatch({
  obs <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")

  if (identical(want, "all")) {
    ai   <- aipw_arms_all(obs, learners = learners, g_floor = G_FLOOR, V_cf = V_CF,
                          inner_cv_folds = INNER_CV)
    rdf  <- ai$results                 # method, b, se_jkn, se_lin, df  (SF, CF, CV)
    Dcv  <- ai$diagnostics$D_cv        # AIPW cluster-CV influence fn (for deff)
    drow <- cbind(ai$diagnostics$drow, ai$diagnostics$drow_cv)
  } else {                              # "cv"
    cv   <- aipw_cv_arm(obs, learners = learners, g_floor = G_FLOOR,
                        inner_cv_folds = INNER_CV)
    rdf  <- cv$results                 # AIPW-CV
    Dcv  <- cv$D_centered
    drow <- cv$drow
  }

  res <- data.frame(method = rdf$method, b = rdf$b, se = rdf$se_lin, df = rdf$df,
                    se_jkn = rdf$se_jkn, estimator = "AIPW", stringsAsFactors = FALSE)

  ch <- attr(obs, "checks")
  dd <- deff_clust(Dcv, obs$strata, obs$cluster, obs$weight)
  list(
    results = cbind(rep = i, res, deff = dd$deff_clust, icc_eif = dd$icc_eif,
                    n = ch$n, df_design = ch$df_design, sumw_over_N = ch$sumw_over_N,
                    min_psu_str = ch$min_psu_str, w_cv = ch$w_cv),
    diag    = cbind(rep = i, drow, deff = dd$deff_clust, icc_eif = dd$icc_eif)
  )
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# ---- run the chunk across cores (parallel-safe L'Ecuyer streams) -------
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
                       library(ranger); library(earth); library(gam); library(surveyCV) })
  })
  clusterSetRNGStream(cl, iseed = iseed)
  clusterExport(cl, c("pop", "scenario", "V_CF", "INNER_CV", "G_FLOOR"), envir = environment())
  reps_out <- parLapply(cl, reps, one_rep, pop = pop, learners = LEARNERS, want = WANT)
} else {
  set.seed(iseed)
  reps_out <- lapply(reps, one_rep, pop = pop, learners = LEARNERS, want = WANT)
}
reps_out <- Filter(Negate(is.null), reps_out)
n_failed <- length(reps) - length(reps_out)
if (!length(reps_out)) stop(sprintf("[%s task %d] all %d reps in this chunk failed",
                                    RUN_ID, task, length(reps)))
if (n_failed > 0)
  cat(sprintf("[%s task %d] WARNING: %d/%d reps failed; chunk RDS checkpoints as complete -- delete %s to recompute\n",
              RUN_ID, task, n_failed, length(reps), fn))
res   <- do.call(rbind, lapply(reps_out, `[[`, "results"))
res$rung   <- rung
res$base_m <- NA_integer_
diags <- do.call(rbind, lapply(reps_out, `[[`, "diag"))

# ---- per-task output (full detail; private folder) ---------------------
out <- list(run_id = RUN_ID, scenario = scenario, rung = rung, base_m = NA_integer_,
            learners = LEARNERS, arms = WANT, chunk = chunk, reps = reps, n_failed = n_failed,
            g_floor = G_FLOOR, V_cf = V_CF, inner_cv = INNER_CV,
            Psi = Psi, truth = pop$truth, params = pop$params,
            results = res, diagnostics = diags)
saveRDS(out, fn)
cat(sprintf("[%s task %d] saved %s  (%d est rows, %d diag rows)\n",
            RUN_ID, task, fn, nrow(res), nrow(diags)))

# ---- reproducibility manifest (one per task; manifest/ subdir) ----------
manifest <- list(
  run_id = RUN_ID, scenario = scenario, rung = rung, arms = WANT, base_m = NA_integer_,
  learners = LEARNERS, chunk = chunk, reps = reps, n_failed = n_failed,
  g_floor = G_FLOOR, V_cf = V_CF, inner_cv = INNER_CV,
  params = pop$params, truth = pop$truth,
  pop_audit = audit_population(pop),
  seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE,
               TRUTH_SEED = TRUTH_SEED, iseed_stream = iseed),
  run_dir = .this_dir,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]],
  smoke = SMOKE
)
mf <- file.path(MAN_DIR, sprintf("%smanifest_%s_%s_chunk%03d.rds", tag, scenario, rung, chunk))
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
  ord_m <- c("AIPW-SF", "AIPW-CV", "AIPW-CF")
  summ <- summ[order(match(summ$method, ord_m)), ]
  np <- sapply(summ, is.numeric); sp <- summ; sp[np] <- round(sp[np], 4)
  cat(sprintf("\n[SMOKE] inline summary (standard %s, AIPW 'all'; true Psi=%.6f):\n",
              paste(LIBS$SMOKE_cheap, collapse = "+"), Psi))
  print(sp, row.names = FALSE)

  finite_ok <- all(is.finite(res$b)) && all(is.finite(res$se))
  meths     <- unique(res$method)
  want_m    <- c("AIPW-SF", "AIPW-CV", "AIPW-CF")
  arms_ok   <- all(want_m %in% meths)
  aipw_cv   <- "AIPW-CV" %in% meths
  reasons <- paste(c(
    sprintf("b/se finite=%s", finite_ok),
    sprintf("AIPW-SF/CV/CF present=%s", arms_ok),
    sprintf("AIPW-CV present=%s", aipw_cv),
    sprintf("methods={%s}", paste(meths, collapse = ","))), collapse = "; ")
  if (finite_ok && arms_ok && aipw_cv)
    cat(sprintf("\n[SMOKE-GATE] PASS: %s\n", reasons))
  else
    cat(sprintf("\n[SMOKE-GATE] STOP: %s\n", reasons))
}
