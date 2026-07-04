# =====================================================================
# A21_rr_or_continuous/run.R  —  RR / OR + continuous-outcome DEMONSTRATION
# (NON-DESTRUCTIVE ARC-style run).
#
# Validates that the delta-method estimands achieve nominal DESIGN coverage,
# i.e. the executable evidence for the Web-Appendix RR/OR corollary:
#   binary outcome     -> risk ratio (RR) + odds ratio (OR)
#   continuous outcome -> risk difference (RD) + ratio-of-means (RR)
# for the single-fit (Fully-Aware) and cross-fit (Fully-Aware-CF) arms, on the
# standard design, across rungs {L1_param, L2_smooth}.
#
# NON-DESTRUCTIVE: sources the engine read-only; the only new pieces are
# estimands_rr_or.R::delta_rr_or() and the per-arm components exposed by
# run_estimators() in diagnostics$arms.
#
# TWO RUN MODES (one file):
#   * LOCAL  : A21_ALL=1  -> loop ALL (outcome x rung x chunk) units in-process,
#              building each population ONCE; per-chunk RDS = resumable progress.
#   * ARC    : SLURM_ARRAY_TASK_ID=k -> the k-th (outcome x rung x chunk) unit.
#
# Env (local fallbacks in brackets):
#   SIM_N_REPS [800]  SIM_CHUNK [50]  SIM_CODE [codes]  SLURM_CPUS_PER_TASK [cores-1]
#   A21_OUT [DATA_ROOT/arc_runs/A21_rr_or_continuous]   SMOKE=1 -> tiny plumbing run
# =====================================================================

RUN_ID <- "A21_rr_or_continuous"

CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE"); CODE <- hit[1]
}
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
source(file.path(CODE, "estimands_rr_or.R"))
source(file.path(CODE, "arc_runs", "_checkpoint.R"))
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")
A21_ALL <- identical(Sys.getenv("A21_ALL"), "1")
cores  <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS <- geti("SIM_N_REPS", 800L)
CHUNK  <- geti("SIM_CHUNK", 50L)
OUT    <- Sys.getenv("A21_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
MAN_DIR <- file.path(OUT, "manifest")
for (d in c(OUT, MAN_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- gamma0: keep the CONTINUOUS outcome positive so ratio-of-means RR is
#      interpretable; binary keeps the headline intercept. -----------------------
g0_of <- function(ot) if (ot == "continuous") 5 else -2
fam_of <- function(ot) if (ot == "continuous") "gaussian" else "binomial"

# ---- cell grid: outcome_type x rung (scenario fixed = standard) ---------------
RUNGS <- strsplit(Sys.getenv("A21_RUNGS", "L1_param,L2_smooth"), ",")[[1]]
cells <- expand.grid(scenario = "standard",
                     outcome_type = c("binary", "continuous"),
                     rung = RUNGS, stringsAsFactors = FALSE)
if (SMOKE) { N_REPS <- 8L; CHUNK <- 8L; cores <- min(cores, 2L)
             cells <- cells[cells$rung == "L1_param", , drop = FALSE] }
n_chunks <- ceiling(N_REPS / CHUNK)

# =====================================================================
# one replication -> per-(method,estimand) point/CI/cover rows
# =====================================================================
one_rep <- function(i, pop, learners, fam, outcome_type, truth) tryCatch({
  obs  <- draw_sample(pop, sample_seed = SAMPLE_SEED_BASE + i, model_type = "complex")
  est  <- run_estimators(obs, learners = learners, family = fam)
  arms <- est$diagnostics$arms                 # Fully-Aware + Fully-Aware-CF
  bin  <- (outcome_type == "binary")
  rr_rows <- list()
  for (m in names(arms)) {
    a  <- arms[[m]]
    rr <- delta_rr_or(a$psi1, a$psi0, a$eif1, a$eif0,
                      obs$strata, obs$cluster, obs$weight, binary = bin)
    rr$method <- m; rr$rep <- i
    rr$truth  <- ifelse(rr$estimand == "RR", truth$RR, truth$OR)
    rr$cover  <- as.integer(rr$lo <= rr$truth & rr$truth <= rr$hi)
    rr$width  <- rr$hi - rr$lo
    rr_rows[[m]] <- rr[, c("rep","method","estimand","est","log_est","se_log",
                           "lo","hi","width","truth","cover")]
  }
  rd <- est$results[, c("method","b","se","df")]
  rd$estimand <- "RD"; rd$truth <- truth$RD; rd$rep <- i
  crit <- qt(0.975, pmax(1, rd$df))
  rd$cover <- as.integer(abs(rd$b - truth$RD) <= crit * rd$se)
  rd$lo <- rd$b - crit*rd$se; rd$hi <- rd$b + crit*rd$se; rd$width <- rd$hi - rd$lo
  list(rr = do.call(rbind, rr_rows),
       rd = rd[, c("rep","method","estimand","b","se","lo","hi","width","truth","cover")])
}, error = function(e) { message(sprintf("rep %d failed: %s", i, conditionMessage(e))); NULL })

# =====================================================================
# a worker cluster with the engine sourced ONCE (reused across all chunks
# so we don't pay makeCluster + source overhead per chunk)
# =====================================================================
make_cl <- function(cores) {
  cl <- makeCluster(cores)
  clusterEvalQ(cl, {
    CODE <- Sys.getenv("SIM_CODE", "R")
    if (!file.exists(file.path(CODE, "config.R"))) CODE <- "R"
    source(file.path(CODE, "config.R")); source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R")); source(file.path(CODE, "diagnostics.R"))
    source(file.path(CODE, "learners.R")); source(file.path(CODE, "estimands_rr_or.R"))
    suppressMessages({ library(survey); library(tmle); library(SuperLearner)
                       library(earth); library(gam); library(ranger) })
  })
  cl
}

# =====================================================================
# run one (outcome x rung x chunk) unit: parallel reps -> per-chunk RDS
# =====================================================================
run_unit <- function(scenario, outcome_type, rung, chunk, pop, truth, cell_index, cl = NULL) {
  tag <- if (SMOKE) "SMOKE_" else ""
  fn  <- file.path(OUT, sprintf("%sa21_%s_%s_%s_chunk%03d.rds",
                                tag, scenario, outcome_type, rung, chunk))
  if (!SMOKE && file.exists(fn) && file.size(fn) > 200) {
    cat(sprintf("[%s] SKIP (done): %s\n", RUN_ID, basename(fn))); return(invisible(NULL))
  }
  learners <- SL_LADDER[[rung]]; fam <- fam_of(outcome_type)
  rep_lo <- (chunk - 1L) * CHUNK + 1L; rep_hi <- min(chunk * CHUNK, N_REPS); reps <- rep_lo:rep_hi
  iseed  <- SAMPLE_SEED_BASE + 1000L * cell_index + chunk
  cat(sprintf("[%s] %s/%s/%s chunk %d reps %d-%d (%d cores)\n",
              RUN_ID, scenario, outcome_type, rung, chunk, rep_lo, rep_hi, cores))
  own_cl <- FALSE
  if (is.null(cl) && cores > 1L) { cl <- make_cl(cores); own_cl <- TRUE }
  if (!is.null(cl)) {
    clusterSetRNGStream(cl, iseed = iseed)
    ro <- parLapply(cl, reps, one_rep, pop = pop, learners = learners, fam = fam,
                    outcome_type = outcome_type, truth = truth)
    if (own_cl) stopCluster(cl)
  } else {
    set.seed(iseed)
    ro <- lapply(reps, one_rep, pop = pop, learners = learners, fam = fam,
                 outcome_type = outcome_type, truth = truth)
  }
  ro <- Filter(Negate(is.null), ro)
  n_failed <- length(reps) - length(ro)
  if (!length(ro)) { warning(sprintf("[%s] all reps failed in %s", RUN_ID, basename(fn))); return(invisible(NULL)) }
  out <- list(run_id = RUN_ID, scenario = scenario, outcome_type = outcome_type, rung = rung,
              learners = learners, chunk = chunk, reps = reps, n_failed = n_failed, truth = truth,
              rr = do.call(rbind, lapply(ro, `[[`, "rr")),
              rd = do.call(rbind, lapply(ro, `[[`, "rd")))
  saveRDS(out, fn)
  manifest <- list(run_id = RUN_ID, scenario = scenario, outcome_type = outcome_type, rung = rung,
                   learners = learners, chunk = chunk, reps = reps, truth = truth,
                   seeds = list(POP_SEED = POP_SEED, SAMPLE_SEED_BASE = SAMPLE_SEED_BASE, iseed = iseed),
                   R_version = R.version.string,
                   packages = sapply(c("SuperLearner","tmle","survey","earth","gam","ranger"),
                                     function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
                   timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), git = arc_git_sha(),
                   sysname = Sys.info()[["nodename"]], smoke = SMOKE)
  saveRDS(manifest, file.path(MAN_DIR, sprintf("%smanifest_%s_%s_%s_chunk%03d.rds",
                                               tag, scenario, outcome_type, rung, chunk)))
  cat(sprintf("[%s] saved %s (%d rr, %d rd rows, %d failed)\n",
              RUN_ID, basename(fn), nrow(out$rr), nrow(out$rd), n_failed))
  if (SMOKE) {
    cov <- aggregate(cover ~ method + estimand, out$rr, mean)
    cat("[SMOKE] RR/OR coverage:\n"); print(cov, row.names = FALSE)
    cat(sprintf("[SMOKE-GATE] %s: finite=%s arms={%s}\n",
                if (all(is.finite(out$rr$est))) "PASS" else "STOP",
                all(is.finite(out$rr$est)), paste(unique(out$rr$method), collapse = ",")))
  }
  invisible(out)
}

# ---- truth (per population, built once per outcome_type) ---------------------
truth_of <- function(pop, outcome_type) {
  tr <- pop$truth
  list(RD = tr$psi, psi1 = tr$psi1, psi0 = tr$psi0, RR = tr$psi1 / tr$psi0,
       OR = if (outcome_type == "binary") (tr$psi1 / (1 - tr$psi1)) / (tr$psi0 / (1 - tr$psi0)) else NA_real_)
}

# =====================================================================
# driver: LOCAL (loop all units, pop built once per outcome) or ARC (one task)
# =====================================================================
if (A21_ALL || SMOKE) {
  CL <- if (cores > 1L) make_cl(cores) else NULL
  for (ot in unique(cells$outcome_type)) {
    pop <- make_population("standard", model_type = "complex", outcome_type = ot,
                           gamma0 = g0_of(ot), pop_seed = POP_SEED,
                           truth_M = if (SMOKE) 2e5L else 2e6L)
    tr  <- truth_of(pop, ot)
    cat(sprintf("\n[%s] === outcome=%s  RD=%.4f psi1=%.4f psi0=%.4f RR=%.4f OR=%s ===\n",
                RUN_ID, ot, tr$RD, tr$psi1, tr$psi0, tr$RR,
                if (is.na(tr$OR)) "NA" else sprintf("%.4f", tr$OR)))
    cot <- cells[cells$outcome_type == ot, , drop = FALSE]
    for (r in seq_len(nrow(cot))) {
      rung <- cot$rung[r]; ci <- which(cells$outcome_type == ot & cells$rung == rung)
      for (ch in seq_len(n_chunks)) run_unit("standard", ot, rung, ch, pop, tr, ci, cl = CL)
    }
  }
  if (!is.null(CL)) stopCluster(CL)
} else {
  grid <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i)
    cbind(cells[i, , drop = FALSE], cell_index = i, chunk = seq_len(n_chunks))))
  task <- geti("SLURM_ARRAY_TASK_ID", 1L); stopifnot(task >= 1L, task <= nrow(grid))
  job  <- grid[task, ]
  pop  <- make_population("standard", model_type = "complex", outcome_type = job$outcome_type,
                          gamma0 = g0_of(job$outcome_type), pop_seed = POP_SEED, truth_M = 2e6L)
  tr   <- truth_of(pop, job$outcome_type)
  fn   <- file.path(OUT, sprintf("a21_standard_%s_%s_chunk%03d.rds", job$outcome_type, job$rung, job$chunk))
  if (file.exists(fn) && file.size(fn) > 200) arc_skip_if_done(fn, task)
  run_unit("standard", job$outcome_type, job$rung, job$chunk, pop, tr, job$cell_index)
}
cat(sprintf("\n[%s] done.\n", RUN_ID))
