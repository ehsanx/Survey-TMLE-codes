# =====================================================================
# Nhanes/arc_runs/R17_e4_floor/run.R
#   R17: E4 floor sensitivity + share-at-floor for N6   (spec item A9b,
#        Writing/comments/phase4-arc-sim-specs.md)
#
# GOAL. For E4 (6% exposure; single-fit min ghat -> ~0.003 on the N1-corrected
# data) the locked 0.05 OOF propensity floor BINDS. Re-run ONLY the primary
# primary estimator (Fully-Aware-CF at the primary 3-learner library) with the
# floor swept over {0.05, 0.025, 0.01} on E4, and record the share-at-floor /
# weighted-exposure-mass-below-0.05 per split (RAW pre-clipping OOF ghat). The
# E2/E3 cells at the locked 0.05 floor give the cross-example N6 share table.
#
# CELL GRID (one SLURM array task = one cell; FULL --array=1-5):
#   task 1: E2, floor 0.05     (N6 share row)
#   task 2: E3, floor 0.05     (N6 share row)
#   task 3: E4, floor 0.05     (N6 share row + sensitivity anchor = locked floor)
#   task 4: E4, floor 0.025    (sensitivity)
#   task 5: E4, floor 0.01     (sensitivity)
# Per cell: B = NHANES_B (default 20) split-repeats of the primary arm with
# split seeds 20260607 + 1:B (the locked nhanes_arc.R convention).
#
# NON-DESTRUCTIVE: sources codes/{estimators,diagnostics}.R + _checkpoint.R
# read-only; ALL new logic is in this folder's estimators_r17.R. Writes ONLY to
#   Nhanes/nhanes_output/arc_runs/R17/   nh17_<example>_f<###>.rds (+ manifest/)
#   results/arc/                         (written by aggregate_R17.R, not here)
# Never touches the locked Nhanes/nhanes_output/{intermediate,results}.
#
# PRIMARY library == Nhanes/R/03_run_estimators.R:
#   LIB = c("SL.glm", "SL.earth", "SL.glmnet")
#
# PROVENANCE: Nhanes/analytic/E*_imputed.rds are the N1-CORRECTED frames (post
# RHD180 fix, commit dcce692; pre-fix copies live in Nhanes/analytic_backup_N1).
#
# REPRODUCIBILITY / COMPARABILITY:
#   * The analytic frames are FIXED data (no draw_sample); the inputs are
#     byte-identical to the locked headline ladder + the ITEM-0 E4 re-run.
#   * Split seeds 20260607 + 1:B match the locked B-split convention. set.seed()
#     is called per split INSIDE the worker, and make_cf_folds is the first RNG
#     consumer in run_cf_floor, so split b is fully deterministic regardless of
#     scheduling AND the raw OOF ghat is IDENTICAL across the three E4 floor
#     cells for the same split (the floor enters only post-hoc) -> exact floor
#     attribution. NOTE: realized folds differ from the locked all-arm run for
#     the same seed (there, FA/NA/CV consume RNG before the CF fold draw).
#   * No clusterSetRNGStream: the per-split set.seed overrides any stream (same
#     as the R05/nhanes_arc.R exemplars), so cells cannot collide by design.
#
# Env vars (mirror R05_harmonized_floor; local-test fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores                         (default detectCores-1)
#   NHANES_B             split repeats                  (default 20; SMOKE -> 2)
#   SIM_CODE             engine dir (codes/)            (default codes)
#   NH_ANA               analytic dir (imputed rds)     (default Nhanes/analytic)
#   R17_DIR              this run's folder              (default REPO_ROOT/Nhanes/arc_runs/R17_e4_floor)
#   R17_OUT              per-cell output dir            (default .../nhanes_output/arc_runs/R17)
#   R17_MANIFEST         manifest dir                   (default <R17_OUT>/manifest)
#   R17_LIB_OVERRIDE     comma-separated SL library override (wiring checks only)
#   SMOKE                "1" -> E4/f050 cell, B=2, GLM-only library, <8 min gate
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "R")
HERE <- Sys.getenv("R17_DIR",
                   file.path(Sys.getenv("REPO_ROOT", "."),
                             "nhanes", "arc_runs", "R17_e4_floor"))
source(file.path(CODE, "estimators.R"))               # building blocks (read-only)
source(file.path(CODE, "diagnostics.R"))              # deff_clust
source(file.path(HERE, "estimators_r17.R"))           # run_cf_floor (NEW logic)
source(file.path(CODE, "arc_runs", "_checkpoint.R"))  # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# one-way ANOVA ICC (deff_clust needs it; defined here so we need not source
# dgp.R) -- byte-identical to the helper in Nhanes/R/nhanes_arc.R.
icc_anova <- function(y, cluster) {
  y <- as.numeric(y); cl <- as.factor(cluster)
  ok <- !is.na(y) & !is.na(cl); y <- y[ok]; cl <- droplevels(cl[ok])
  if (nlevels(cl) < 2L) return(NA_real_)
  ni <- as.numeric(table(cl)); k <- length(ni); N <- length(y); gm <- mean(y)
  mu <- tapply(y, cl, mean)
  msb <- sum(ni * (mu - gm)^2) / (k - 1)
  msw <- sum((y - mu[cl])^2) / (N - k)
  n0  <- (N - sum(ni^2) / N) / (k - 1)
  vb  <- (msb - msw) / n0
  max(0, vb) / (max(0, vb) + msw)
}

geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
SMOKE <- identical(Sys.getenv("SMOKE"), "1")
task  <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
B     <- geti("NHANES_B", 20L)
ANA   <- Sys.getenv("NH_ANA", "nhanes/analytic")
OUT   <- Sys.getenv("R17_OUT",
                    file.path(Sys.getenv("REPO_ROOT", "."),
                              "nhanes", "nhanes_output", "arc_runs", "R17"))
MAN   <- Sys.getenv("R17_MANIFEST", file.path(OUT, "manifest"))
for (d in c(OUT, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- PRIMARY library (== Nhanes/R/03_run_estimators.R LIB) ---------
# R17_LIB_OVERRIDE exists for cheap wiring checks ONLY; the full LIB runs on ARC.
LIB_CERT <- c("SL.glm", "SL.earth", "SL.glmnet")
ov  <- Sys.getenv("R17_LIB_OVERRIDE", "")
learners <- if (nzchar(ov)) trimws(strsplit(ov, ",")[[1]]) else LIB_CERT
V_CF <- 5L                                   # locked cross-fit V (V_eff capped by min PSUs/stratum)
SPLIT_SEED_BASE <- 20260607L                 # locked B-split seed base (nhanes_arc.R)

# ---- cell grid: (example, floor)  [FULL --array=1-5] -------------------------
CELLS <- data.frame(example = c("E2", "E3", "E4",  "E4",  "E4"),
                    floor   = c(0.05, 0.05, 0.05, 0.025, 0.01),
                    stringsAsFactors = FALSE)

if (SMOKE) {
  # wiring-check cell: E4 / floor 0.05 (the locked floor on the binding example),
  # B=2, <=2 cores, GLM-only library unless an explicit override was given.
  task  <- 3L
  B     <- 2L
  cores <- max(1L, min(cores, 2L))
  if (!nzchar(ov)) learners <- "SL.glm"
  cat("[SMOKE] E4/f050 single cell, B=2, learners={", paste(learners, collapse = ","),
      "} (wiring check, NOT the primary library)\n", sep = "")
}

stopifnot(task >= 1L, task <= nrow(CELLS))
ex_id   <- CELLS$example[task]
g_floor <- CELLS$floor[task]
ftag    <- sprintf("f%03d", as.integer(round(g_floor * 1000)))  # 0.05->f050, 0.025->f025, 0.01->f010

# ---- per-cell checkpoint: skip instantly if this cell already finished -------
tag <- if (SMOKE) "SMOKE_" else ""
fn  <- file.path(OUT, sprintf("%snh17_%s_%s.rds", tag, ex_id, ftag))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- load the (already imputed, N1-corrected) analytic data; encode covs -----
# IDENTICAL encoding to Nhanes/R/nhanes_arc.R / 03_run_estimators.R so the design
# matrix + W_cols match the locked pipeline byte-for-byte.
obs  <- readRDS(file.path(ANA, paste0(ex_id, "_imputed.rds")))
covs <- attr(obs, "covs"); label <- attr(obs, "example"); sub <- which(obs$inpop)
covdf <- droplevels(obs[sub, covs, drop = FALSE])
mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
stopifnot(nrow(mm) == length(sub))
wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
for (j in seq_along(wnames)) { obs[[wnames[j]]] <- NA_real_; obs[[wnames[j]]][sub] <- mm[, j] }
obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED

cat(sprintf("[task %d] %s (%s) | floor=%.3f (%s) | learners={%s} | n_domain=%d | A_prev=%.4f | B=%d | cores=%d\n",
            task, ex_id, label, g_floor, ftag, paste(learners, collapse = ","),
            length(sub), mean(obs$A[sub]), B, cores))

# ---- B split-repeats of the primary arm across cores -----------------------
# Per-split tryCatch -> NULL on error; a cell fails only if ALL splits fail.
one_run <- function(b) tryCatch({
  seed <- SPLIT_SEED_BASE + b               # locked split-seed convention
  set.seed(seed)                            # first RNG consumer = make_cf_folds
  r <- run_cf_floor(obs, learners = learners, V_cf = V_CF, g_floor = g_floor,
                    W_cols = wnames, nest = TRUE, inpop = obs$inpop)
  r$row <- cbind(split = b, seed = seed, r$row)
  r
}, error = function(e) {
  message(sprintf("split %d failed: %s", b, conditionMessage(e))); NULL
})

ncl <- max(1L, min(cores, B))
cl  <- makeCluster(ncl); on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "R")
  HERE <- Sys.getenv("R17_DIR",
                     file.path(Sys.getenv("REPO_ROOT", "."),
                               "nhanes", "arc_runs", "R17_e4_floor"))
  source(file.path(CODE, "estimators.R"))
  source(file.path(HERE, "estimators_r17.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger);
                     library(earth); library(gam); library(glmnet) })
  TRUE
}))
# Export EVERYTHING the worker closure needs (one_run + run_cf_floor inputs):
clusterExport(cl, c("obs", "learners", "wnames", "g_floor", "V_CF",
                    "SPLIT_SEED_BASE", "run_cf_floor", "one_run"),
              envir = environment())
runs <- parLapply(cl, seq_len(B), one_run)
runs <- Filter(Negate(is.null), runs)
if (!length(runs)) stop(sprintf("[task %d] all %d splits in this cell failed", task, B))
B_eff <- length(runs)
if (B_eff < B) cat(sprintf("[task %d] WARNING: %d/%d splits failed (B_eff=%d)\n",
                           task, B - B_eff, B, B_eff))

# ---- summarize across the B splits (locked column conventions) ---------------
per_split <- do.call(rbind, lapply(runs, `[[`, "row"))
summ <- data.frame(
  example = ex_id, label = label, floor = g_floor, method = "Fully-Aware-CF",
  B = B_eff,
  b = mean(per_split$b), b_split_sd = sd(per_split$b),
  se = mean(per_split$se), df = per_split$df[1],
  share_clip_w   = mean(per_split$share_clip_w),
  share_clip_unw = mean(per_split$share_clip_unw),
  mass05_w       = mean(per_split$mass05_w),
  expmass05_w    = mean(per_split$expmass05_w),
  g_raw_min      = mean(per_split$g_raw_min),
  g_raw_max      = mean(per_split$g_raw_max))
crit <- qt(0.975, pmax(1, summ$df))
summ$lcl <- summ$b - crit * summ$se; summ$ucl <- summ$b + crit * summ$se

# ---- rich diagnostics (from the first successful split) ----------------------
r1 <- runs[[1]]
ucl_sub <- paste(obs$strata[sub], obs$cluster[sub], sep = "_")  # globally-unique PSU id
dd <- deff_clust(r1$eif, obs$strata[sub], ucl_sub, obs$weight[sub])
diagnostics <- list(
  deff_clust = dd$deff_clust, icc_eif = dd$icc_eif,
  g_raw = r1$g_raw,                                             # full vector for overlap/N6 figures
  g_raw_q = quantile(r1$g_raw, c(.001, .01, .05, .25, .5, .75, .95, .99, .999)),
  n_domain = length(sub), A_prev = mean(obs$A[sub]), Y_prev = mean(obs$Y[sub]),
  n_psu = length(unique(paste(obs$strata[sub], obs$cluster[sub]))),
  n_strata = length(unique(obs$strata[sub])),
  design_df = summ$df[1],
  min_psu_stratum = min(table(unique(data.frame(s = obs$strata[sub], c = obs$cluster[sub]))$s)),
  weight_cv = sd(obs$weight[sub]) / mean(obs$weight[sub]),
  g_floor = g_floor, q_bound = 1e-3, v_cf = V_CF)

out <- list(run = "R17_e4_floor", example = ex_id, label = label,
            floor = g_floor, floor_tag = ftag, learners = learners,
            B = B, B_eff = B_eff, seeds = SPLIT_SEED_BASE + seq_len(B),
            summary = summ, per_split = per_split, diagnostics = diagnostics)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s\n", task, fn))
print(summ[, c("example", "floor", "method", "b", "b_split_sd", "se", "df", "lcl", "ucl")],
      row.names = FALSE)
cat(sprintf("[task %d] share_clip_w=%.4f share_clip_unw=%.4f mass05_w=%.4f expmass05_w=%.4f | g_raw=[%.4f,%.4f] used=[%.4f,%.4f] | DEFF(eif)=%.2f\n",
            task, summ$share_clip_w, summ$share_clip_unw, summ$mass05_w, summ$expmass05_w,
            min(per_split$g_raw_min), max(per_split$g_raw_max),
            min(per_split$g_used_min), max(per_split$g_used_max), dd$deff_clust))

# ---- reproducibility manifest (mirrors nhanes_arc.R / R05) --------------------
manifest <- list(
  run_id = "R17_e4_floor",
  example = ex_id, label = label, floor = g_floor, floor_tag = ftag,
  learners = learners, lib_primary = LIB_CERT, lib_override = nzchar(ov),
  B = B, B_eff = B_eff, V_cf = V_CF, q_bound = 1e-3,
  n_domain = length(sub), covariates = covs, dummy_cols = wnames,
  seeds = list(SPLIT_SEED_BASE = SPLIT_SEED_BASE, split_seeds = SPLIT_SEED_BASE + seq_len(B)),
  data_provenance = "Nhanes/analytic/*_imputed.rds, N1-corrected (post RHD180 fix)",
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]])
saveRDS(manifest, file.path(MAN, sprintf("%smanifest_%s_%s.rds", tag, ex_id, ftag)))

# ---- SMOKE: inline summary + explicit gate verdict ----------------------------
if (SMOKE) {
  cat("\n[SMOKE] per-split records (E4/f050, wiring library):\n")
  show <- per_split[, c("split", "b", "se", "df", "share_clip_w", "share_clip_unw",
                        "mass05_w", "expmass05_w", "g_raw_min", "g_raw_max",
                        "g_used_min", "eps_cf", "cf_V_eff")]
  print(format(show, digits = 4), row.names = FALSE)
  cat("\n[SMOKE] cell summary:\n")
  print(format(summ, digits = 4), row.names = FALSE)

  checks <- c(
    "finite b"                 = all(is.finite(per_split$b)),
    "finite se > 0"            = all(is.finite(per_split$se) & per_split$se > 0),
    "df > 0"                   = all(per_split$df > 0),
    "share_clip_w in [0,1]"    = all(per_split$share_clip_w >= 0 & per_split$share_clip_w <= 1),
    "share_clip_unw in [0,1]"  = all(per_split$share_clip_unw >= 0 & per_split$share_clip_unw <= 1),
    "mass05_w in [0,1]"        = all(per_split$mass05_w >= 0 & per_split$mass05_w <= 1),
    "expmass05_w in [0,1]"     = all(per_split$expmass05_w >= 0 & per_split$expmass05_w <= 1),
    "g_raw range in (0,1)"     = all(per_split$g_raw_min > 0 & per_split$g_raw_max < 1 &
                                     per_split$g_raw_min < per_split$g_raw_max),
    "floor binds (raw<floor)"  = all(per_split$g_raw_min < g_floor),   # E4: floor must bind
    # NOTE: the g_used line is a cheap sanity invariant on OUR clipping only --
    # tmle 2.1.1 returns the user-supplied g1W unmodified (internal bounding hits
    # an unexposed g1W.total), so it canNOT detect a passthrough failure; floor
    # correctness rests on the source-inspection argument (default gbound
    # 5/sqrt(n)/log(n) ~ 0.0051 at E4's n < all three floors => .bound is a no-op).
    "g used >= floor (clip sanity)" = all(per_split$g_used_min >= g_floor - 1e-12),
    "eps_cf not diverged"      = all(is.finite(per_split$eps_cf) & per_split$eps_cf < 20),
    "no split failed"          = (B_eff == B))
  cat("\n[SMOKE] gate checks:\n")
  for (i in seq_along(checks))
    cat(sprintf("  %-38s %s\n", names(checks)[i], if (checks[i]) "OK" else "FAIL"))
  cat(sprintf("[SMOKE] share numbers: share_clip_w=%.4f share_clip_unw=%.4f mass05_w=%.4f expmass05_w=%.4f\n",
              summ$share_clip_w, summ$share_clip_unw, summ$mass05_w, summ$expmass05_w))
  if (all(checks)) {
    cat("[SMOKE-GATE] PASS: wiring OK (finite b/se, shares in [0,1], floor binds on raw ghat; clip sanity holds -- floor passthrough itself rests on the source-inspection argument, not this gate). Full primary-LIB runs on ARC only.\n")
  } else {
    cat(sprintf("[SMOKE-GATE] STOP: failed checks: %s\n",
                paste(names(checks)[!checks], collapse = "; ")))
    quit(save = "no", status = 1)
  }
}
