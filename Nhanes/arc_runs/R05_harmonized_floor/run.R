# =====================================================================
# Nhanes/arc_runs/R05_harmonized_floor/run.R
#   R05: NHANES harmonized-floor + weighted-OOF-CF ladder (E3/E4)
#
# Mirrors Nhanes/R/nhanes_arc.R (split-repeat B + rich diagnostics) but swaps the
# canonical run_estimators() for THIS folder's run_estimators_iso() (helper:
# estimators_nhanes_iso.R), which re-emits the locked Fully-Aware (1e-3) and
# Fully-Aware-CF (unweighted-OOF, 0.05) arms UNCHANGED and ADDS two de-confounding
# arms: Fully-Aware-h05 (FA targeting, EIF re-bounded at 0.05) and
# Fully-Aware-CF-wOOF (cross-fit with WEIGHTED per-fold nuisances).
#
# NON-DESTRUCTIVE: sources codes/estimators.R etc. read-only; writes ONLY to
#   Nhanes/nhanes_output/arc_runs/R05/   (per-cell RDS)
#   results/arc/R05_summary.csv          (aggregated, written by the aggregator)
# Never touches the locked Nhanes/nhanes_output/{intermediate,results}.
#
# GRID: example {E3,E4} x rung {L1_param,L2_smooth,L3_adaptive,L4_aggressive} = 8.
#   array task t (1..8): example = c("E3","E4")[((t-1) %% 2) + 1]
#                        rung    = RUNGS[((t-1) %/% 2) + 1]
#   (So tasks 1,2 = L1 E3/E4; 3,4 = L2; 5,6 = L3; 7,8 = L4.) Restrict via --array
#   to e.g. 5-8 for L3,L4 only.
#
# Env vars (mirror Nhanes/nhanes.slurm; local-test fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores                         (default detectCores-1)
#   NHANES_B             cross-fit-split repeats        (default 20; SMOKE -> 2)
#   SIM_CODE             engine dir (codes/)            (default codes)
#   NH_ANA               analytic dir (imputed rds)     (default Nhanes/analytic)
#   NH_OUT               per-task output dir            (default .../arc_runs/R05)
#   NH_MANIFEST          manifest dir                   (default .../arc_runs/R05/manifest)
#   SMOKE                "1" -> tiny run (E4/L4 single cell, B=2)  [decision gate]
# =====================================================================

CODE     <- Sys.getenv("SIM_CODE", "codes")
HERE     <- Sys.getenv("R05_DIR",                                  # this run's folder
                       file.path(Sys.getenv("REPO_ROOT", "."),
                                 "Nhanes", "arc_runs", "R05_harmonized_floor"))
source(file.path(CODE, "estimators.R"))                            # building blocks (read-only)
source(file.path(CODE, "diagnostics.R"))                           # deff_clust
source(file.path(CODE, "learners.R"))                              # SL_LADDER, custom wrappers
source(file.path(HERE, "estimators_nhanes_iso.R"))                 # run_estimators_iso (NEW arms)
source(file.path(CODE, "arc_runs", "_checkpoint.R"))               # arc_skip_if_done(), arc_git_sha()
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# one-way ANOVA ICC (deff_clust needs it; defined here so we need not source dgp.R)
# -- byte-identical to the helper in Nhanes/R/nhanes_arc.R.
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
ANA   <- Sys.getenv("NH_ANA", "Nhanes/analytic")
OUT   <- Sys.getenv("NH_OUT",
                    file.path(Sys.getenv("REPO_ROOT", "."),
                              "Nhanes", "nhanes_output", "arc_runs", "R05"))
MAN   <- Sys.getenv("NH_MANIFEST", file.path(OUT, "manifest"))
FA_GBOUND   <- as.numeric(Sys.getenv("R05_FA_GBOUND",  "0.05"))    # harmonized FA floor
G_OOF_BOUND <- as.numeric(Sys.getenv("R05_OOF_BOUND",  "0.05"))    # OOF propensity floor (matches locked)
for (d in c(OUT, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: {E3,E4} x ladder (FULL L1..L4 examples) ---------------------
EXAMPLES_ID <- c("E3", "E4")
RUNGS <- names(SL_LADDER)                                           # L1..L4 in ladder order
grid  <- expand.grid(example = EXAMPLES_ID, rung = RUNGS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

if (SMOKE) {
  # decision-gate cell: E4 / L4_aggressive (the deepest non-Donsker rung, where
  # the single-fit-vs-CF divergence is largest in the locked results).
  grid  <- grid[grid$example == "E4" & grid$rung == "L4_aggressive", , drop = FALSE]
  task  <- 1L
  B     <- 2L
  cores <- max(1L, min(cores, 2L))
  cat("[SMOKE] E4/L4_aggressive single cell, B=2\n")
}

stopifnot(task >= 1L, task <= nrow(grid))
ex_id <- grid$example[task]; rung <- grid$rung[task]; learners <- SL_LADDER[[rung]]

# ---- per-cell checkpoint: skip instantly if this cell already finished -------
fn <- file.path(OUT, sprintf("nh_%s_%s.rds", ex_id, rung))
if (!SMOKE) arc_skip_if_done(fn, task)

# ---- load the (already imputed) analytic data; encode covariates -----------
# IDENTICAL encoding to Nhanes/R/nhanes_arc.R so the design matrix + W_cols match.
obs  <- readRDS(file.path(ANA, paste0(ex_id, "_imputed.rds")))
covs <- attr(obs, "covs"); label <- attr(obs, "example"); sub <- which(obs$inpop)
covdf <- droplevels(obs[sub, covs, drop = FALSE])
mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
stopifnot(nrow(mm) == length(sub))
wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
for (j in seq_along(wnames)) { obs[[wnames[j]]] <- NA_real_; obs[[wnames[j]]][sub] <- mm[, j] }
obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED

cat(sprintf("[task %d] %s (%s) | rung=%s | learners={%s} | n_domain=%d | B=%d | cores=%d | fa_floor=%.3f oof_floor=%.3f\n",
            task, ex_id, label, rung, paste(learners, collapse = ","),
            length(sub), B, cores, FA_GBOUND, G_OOF_BOUND))

# ---- repeated cross-fit-split runs across cores (mirror nhanes_arc.R) -------
one_run <- function(seed) {
  set.seed(seed)
  r <- run_estimators_iso(obs, learners = learners, V_cf = 5L,
                          g_oof_bound = G_OOF_BOUND, fa_gbound = FA_GBOUND,
                          W_cols = wnames, nest = TRUE, inpop = obs$inpop)
  list(results = r$results, drow = r$diagnostics$drow,
       eif_fa = r$diagnostics$eif_fa, eif_fa_h05 = r$diagnostics$eif_fa_h05,
       g_fa = r$diagnostics$g_fa, g_cf = r$diagnostics$g_cf, g_cfw = r$diagnostics$g_cfw)
}
seeds <- 20260607L + seq_len(B)                                    # SAME seed base as nhanes_arc.R
ncl <- max(1L, min(cores, B))
cl <- makeCluster(ncl); on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "codes")
  HERE <- Sys.getenv("R05_DIR",
                     file.path(Sys.getenv("REPO_ROOT", "."),
                               "Nhanes", "arc_runs", "R05_harmonized_floor"))
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "learners.R"))
  source(file.path(HERE, "estimators_nhanes_iso.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger);
                     library(earth); library(gam); library(glmnet) })
  TRUE
}))
clusterExport(cl, c("obs", "learners", "wnames", "FA_GBOUND", "G_OOF_BOUND"),
              envir = environment())
runs <- parLapply(cl, seeds, one_run)

# ---- summarize across the B splits: per-arm mean estimate + split-SD ---------
# SAME columns/shape as Nhanes/R/nhanes_arc.R so aggregate_nhanes.R works as-is.
allres <- do.call(rbind, Map(function(r, b) cbind(split = b, r$results), runs, seq_along(runs)))
arms <- unique(allres$method)
summ <- do.call(rbind, lapply(arms, function(m) {
  s <- allres[allres$method == m, ]
  data.frame(example = ex_id, label = label, rung = rung, method = m,
             b = mean(s$b), b_split_sd = sd(s$b), se = mean(s$se), df = s$df[1])
}))
crit <- qt(0.975, pmax(1, summ$df)); summ$lcl <- summ$b - crit * summ$se; summ$ucl <- summ$b + crit * summ$se

# ---- rich diagnostics (from split 1) ----------------------------------------
r1 <- runs[[1]]
ucl_sub <- paste(obs$strata[sub], obs$cluster[sub], sep = "_")   # globally-unique PSU id
dd  <- deff_clust(r1$eif_fa,     obs$strata[sub], ucl_sub, obs$weight[sub])  # locked FA EIF
ddh <- deff_clust(r1$eif_fa_h05, obs$strata[sub], ucl_sub, obs$weight[sub])  # harmonized EIF
gb_fa <- FA_GBOUND
diagnostics <- list(
  drow = r1$drow,
  deff_clust = dd$deff_clust, icc_eif = dd$icc_eif,
  deff_clust_h05 = ddh$deff_clust, icc_eif_h05 = ddh$icc_eif,   # harmonized-floor sensitivity
  g_fa = r1$g_fa, g_cf = r1$g_cf, g_cfw = r1$g_cfw,             # full vectors for overlap plots
  g_fa_q  = quantile(r1$g_fa,  c(.01,.05,.25,.5,.75,.95,.99)),
  g_cf_q  = quantile(r1$g_cf,  c(.01,.05,.25,.5,.75,.95,.99)),
  g_cfw_q = quantile(r1$g_cfw, c(.01,.05,.25,.5,.75,.95,.99)),
  # floor-SENSITIVE near-bound mass (reported, NOT the decision statistic):
  g_fa_near_bound = mean(r1$g_fa < gb_fa | r1$g_fa > 1 - gb_fa),
  n_domain = length(sub), A_prev = mean(obs$A[sub]), Y_prev = mean(obs$Y[sub]),
  n_psu = length(unique(paste(obs$strata[sub], obs$cluster[sub]))),
  n_strata = length(unique(obs$strata[sub])),
  design_df = summ$df[summ$method == "Fully-Aware-CF"][1],
  min_psu_stratum = min(table(unique(data.frame(s = obs$strata[sub], c = obs$cluster[sub]))$s)),
  weight_cv = sd(obs$weight[sub]) / mean(obs$weight[sub]),
  fa_gbound = FA_GBOUND, g_oof_bound = G_OOF_BOUND)

out <- list(example = ex_id, label = label, rung = rung, learners = learners,
            B = B, summary = summ, per_split = allres, diagnostics = diagnostics)
saveRDS(out, fn)   # `fn` computed early for the checkpoint guard
cat(sprintf("[task %d] saved %s\n", task, fn))
print(summ[, c("method","b","se","lcl","ucl","df","b_split_sd")], row.names = FALSE)

# ---- DECISION-RULE readout (point-estimate divergence; FLOOR-ROBUST) --------
# Lock the framing on the FA-vs-CF POINT divergence, which cannot move with the
# floor (FA psi is computed from tmle Qstar, not the EIF). FA and FA-h05 share an
# identical b by construction -> their |delta b| is exactly 0 (a built-in check).
b_of <- function(m) summ$b[summ$method == m]
d_fa_cf   <- b_of("Fully-Aware") - b_of("Fully-Aware-CF")          # single-fit vs cross-fit (the headline)
d_fa_h    <- b_of("Fully-Aware") - b_of("Fully-Aware-h05")         # floor effect on the POINT (must be ~0)
d_cfw_cf  <- b_of("Fully-Aware-CF-wOOF") - b_of("Fully-Aware-CF")  # weighting effect within cross-fit
d_cfw_fa  <- b_of("Fully-Aware-CF-wOOF") - b_of("Fully-Aware")     # does weighted-OOF track FA or CF?
cat(sprintf(
  "[DECISION] |b_FA - b_CF|=%.4f (headline divergence)  |b_FA - b_FAh05|=%.4f (floor effect on POINT; ~0)\n           b_CFwOOF=%.4f vs b_CF=%.4f (|d|=%.4f)  vs b_FA=%.4f (|d|=%.4f)\n",
  abs(d_fa_cf), abs(d_fa_h),
  b_of("Fully-Aware-CF-wOOF"), b_of("Fully-Aware-CF"), abs(d_cfw_cf),
  b_of("Fully-Aware"), abs(d_cfw_fa)))

# ---- reproducibility manifest (mirror nhanes_arc.R) -------------------------
manifest <- list(
  run_id = "R05_harmonized_floor",
  example = ex_id, label = label, rung = rung, learners = learners, B = B,
  fa_gbound = FA_GBOUND, g_oof_bound = G_OOF_BOUND,
  n_domain = length(sub), covariates = covs, dummy_cols = wnames,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger","mice"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = arc_git_sha(),
  sysname = Sys.info()[["nodename"]])
saveRDS(manifest, file.path(MAN, sprintf("manifest_%s_%s.rds", ex_id, rung)))
cat(sprintf("[task %d] DEFF(eif,1e-3)=%.2f DEFF(eif,h05)=%.2f | g_FA=[%.3f,%.3f] g_CF(unwt)=[%.3f,%.3f] g_CF(wt)=[%.3f,%.3f] | near-bound(FA@%.2f)=%.1f%%\n",
            task, dd$deff_clust, ddh$deff_clust,
            min(r1$g_fa), max(r1$g_fa), min(r1$g_cf), max(r1$g_cf),
            min(r1$g_cfw), max(r1$g_cfw), gb_fa, 100 * diagnostics$g_fa_near_bound))
