# =====================================================================
# config.R  —  central paths + parallel settings for the survey-TMLE sims
# Sourced by dgp.R, estimators.R, run_smoke.R, and the refactored sim*.R.
# =====================================================================
#
# Storage split:
#   * CODE + small/organized outputs (summary CSVs, figures) -> the repository
#   * LARGE per-rep + aggregated RDS                          -> ./sim_output (gitignored)
# The locked summary CSVs the paper quotes are committed under results/; the large
# per-rep outputs are regenerated locally into ./sim_output.

# ---- roots (env-overridable so the same code runs on Windows + Linux/HPC) -------
# Defaults assume you run from the repository root. All paths are env-overridable;
# on a cluster, sim.slurm sets REPO_ROOT / DATA_ROOT / SIM_CODE via the environment.
REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
DATA_ROOT <- Sys.getenv("DATA_ROOT", file.path(REPO_ROOT, "outputs", "sim_output"))

# ---- locked summary tables (committed) + generated figures (gitignored) ----
RESULTS_DIR <- file.path(REPO_ROOT, "results")              # locked summary CSVs (committed)
IMAGES_DIR  <- file.path(REPO_ROOT, "outputs", "figures")   # generated figures (gitignored)

# ---- full simulation outputs (organized under DATA_ROOT = outputs/sim_output) ----
# Reviewer-proofing: keep ALL internal detail so later "please report X"
# requests can be answered without re-running. Layout:
#   DATA_ROOT/intermediate  per-task RDS: per-rep estimates, SEs, df, AND per-rep
#                           diagnostics (epsilon, g-ranges, DEFF, ICC, design checks)
#   DATA_ROOT/results       aggregated summary CSV + combined RDS
#   DATA_ROOT/manifest      run config, seeds, sessionInfo, package versions,
#                           population audits (ICC, alpha0, sigmas, truth)
DATA_INTERMEDIATE <- file.path(DATA_ROOT, "intermediate")
DATA_RESULTS      <- file.path(DATA_ROOT, "results")
MANIFEST_DIR      <- file.path(DATA_ROOT, "manifest")

# ---- ensure the directory tree exists --------------------------------
.ensure_dirs <- function() {
  for (d in c(RESULTS_DIR, IMAGES_DIR, DATA_INTERMEDIATE, DATA_RESULTS, MANIFEST_DIR)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}
.ensure_dirs()

# ---- parallel settings -----------------------------------------------
# Laptop: 6 physical cores / 12 logical, 15.8 GB RAM.
# Cap workers to protect RAM during SuperLearner-heavy runs (each worker holds
# a copy of the population + its SL fits). 6 is a safe default; raise for the
# cheap GLM smoke, lower if memory pressure appears.
N_CORES_SMOKE <- 6L
N_CORES_FULL  <- 6L
default_cores <- function(cap = N_CORES_FULL) {
  max(1L, min(cap, parallel::detectCores() - 1L))
}

# ---- reproducibility -------------------------------------------------
# Two independent seed roots (see plan-phase2.md B3): population vs sample.
POP_SEED         <- 20260606L   # fixed -> identical population across workers/reps
SAMPLE_SEED_BASE <- 1000L       # per-rep sample seed = SAMPLE_SEED_BASE + rep index
TRUTH_SEED       <- 424242L     # huge-N Monte-Carlo truth integral
RNGkind("L'Ecuyer-CMRG")        # parallel-safe streams

invisible(NULL)
