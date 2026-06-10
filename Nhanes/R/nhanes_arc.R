# =====================================================================
# Nhanes/R/nhanes_arc.R  —  NHANES application as a SLURM ARRAY job on UBC ARC
#
# Mirrors codes/run_sim.R, but on the FOUR fixed real analytic datasets instead of
# simulated draws. One array task = one (example x SL-library-rung) cell. Each task
# runs the full FIVE-arm suite at that rung and, using the cores SLURM gave it,
# repeats the analysis over NHANES_B random cross-fitting / Super-Learner splits to
# quantify split-induced variability (a "fuller" diagnostic). Rich per-cell
# diagnostics are saved (DEFF on the influence function, propensity overlap for the
# single-fit vs cross-fitted estimators, design df, min PSUs/stratum, weight CV,
# targeting epsilon). Aggregate afterwards with aggregate_nhanes.R.
#
# The example x rung grid is the APPLIED echo of the simulation's Figure 1: as the
# library climbs L1->L4 (Donsker -> non-Donsker), single-fit Fully-Aware can drift /
# under-cover while Fully-Aware-CF stays calibrated -- now on real data.
#
# Env vars (local-test fallbacks):
#   SLURM_ARRAY_TASK_ID  1-based task index            (default 1)
#   SLURM_CPUS_PER_TASK  cores                         (default detectCores-1)
#   NHANES_B             cross-fit-split repeats        (default 20)
#   SIM_CODE             engine dir (codes/)            (default codes)
#   NH_ANA               analytic dir (imputed rds)     (default Nhanes/analytic)
#   NH_OUT               per-task output dir            (default Nhanes/nhanes_output/intermediate)
#   NH_MANIFEST          manifest dir                   (default Nhanes/nhanes_output/manifest)
# =====================================================================

CODE <- Sys.getenv("SIM_CODE", "codes")
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
source(file.path(CODE, "learners.R"))
suppressMessages({ library(parallel); library(survey); library(tmle); library(SuperLearner) })

# one-way ANOVA ICC (deff_clust needs it; defined here so we need not source dgp.R)
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

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
task  <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
B     <- geti("NHANES_B", 20L)
ANA   <- Sys.getenv("NH_ANA", "Nhanes/analytic")
OUT   <- Sys.getenv("NH_OUT", "Nhanes/nhanes_output/intermediate")
MAN   <- Sys.getenv("NH_MANIFEST", "Nhanes/nhanes_output/manifest")
for (d in c(OUT, MAN)) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- cell grid: example x rung ---------------------------------------------
EXAMPLES_ID <- c("E1", "E2", "E3", "E4")
RUNGS <- names(SL_LADDER)
grid <- expand.grid(example = EXAMPLES_ID, rung = RUNGS,
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
stopifnot(task >= 1L, task <= nrow(grid))
ex_id <- grid$example[task]; rung <- grid$rung[task]; learners <- SL_LADDER[[rung]]

# ---- load the (already imputed) analytic data; encode covariates -----------
obs  <- readRDS(file.path(ANA, paste0(ex_id, "_imputed.rds")))
covs <- attr(obs, "covs"); label <- attr(obs, "example"); sub <- which(obs$inpop)
covdf <- droplevels(obs[sub, covs, drop = FALSE])
mm <- model.matrix(~ ., data = covdf)[, -1, drop = FALSE]
stopifnot(nrow(mm) == length(sub))
wnames <- make.names(colnames(mm), unique = TRUE); colnames(mm) <- wnames
for (j in seq_along(wnames)) { obs[[wnames[j]]] <- NA_real_; obs[[wnames[j]]][sub] <- mm[, j] }
obs$strata <- obs$SDMVSTRA; obs$cluster <- obs$SDMVPSU; obs$weight <- obs$WTMEC_POOLED

cat(sprintf("[task %d] %s (%s) | rung=%s | learners={%s} | n_domain=%d | B=%d | cores=%d\n",
            task, ex_id, label, rung, paste(learners, collapse = ","), length(sub), B, cores))

# ---- repeated cross-fit-split runs across cores ----------------------------
one_run <- function(seed) {
  set.seed(seed)
  r <- run_estimators(obs, learners = learners, V_cf = 5L, inner_cv_folds = 5L,
                      W_cols = wnames, nest = TRUE, inpop = obs$inpop)
  list(results = r$results, drow = r$diagnostics$drow,
       eif_fa = r$diagnostics$eif_fa, g_fa = r$diagnostics$g_fa, g_cf = r$diagnostics$g_cf)
}
seeds <- 20260607L + seq_len(B)
ncl <- max(1L, min(cores, B))
cl <- makeCluster(ncl); on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, {
  CODE <- Sys.getenv("SIM_CODE", "codes")
  source(file.path(CODE, "estimators.R")); source(file.path(CODE, "learners.R"))
  suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(ranger);
                     library(earth); library(gam); library(glmnet) })
  TRUE
}))
clusterExport(cl, c("obs", "learners", "wnames"), envir = environment())
runs <- parLapply(cl, seeds, one_run)

# ---- summarize across the B splits: per-arm mean estimate + split-SD --------
allres <- do.call(rbind, Map(function(r, b) cbind(split = b, r$results), runs, seq_along(runs)))
arms <- unique(allres$method)
summ <- do.call(rbind, lapply(arms, function(m) {
  s <- allres[allres$method == m, ]
  data.frame(example = ex_id, label = label, rung = rung, method = m,
             b = mean(s$b), b_split_sd = sd(s$b), se = mean(s$se), df = s$df[1])
}))
crit <- qt(0.975, pmax(1, summ$df)); summ$lcl <- summ$b - crit * summ$se; summ$ucl <- summ$b + crit * summ$se

# ---- rich diagnostics (from split 1) ---------------------------------------
r1 <- runs[[1]]
ucl_sub <- paste(obs$strata[sub], obs$cluster[sub], sep = "_")   # globally-unique PSU id (SDMVPSU repeats)
dd <- deff_clust(r1$eif_fa, obs$strata[sub], ucl_sub, obs$weight[sub])
gbnd <- 0.05
diagnostics <- list(
  drow = r1$drow,
  deff_clust = dd$deff_clust, icc_eif = dd$icc_eif,
  g_fa = r1$g_fa, g_cf = r1$g_cf,                       # full vectors for overlap plots
  g_fa_q = quantile(r1$g_fa, c(.01,.05,.25,.5,.75,.95,.99)),
  g_cf_q = quantile(r1$g_cf, c(.01,.05,.25,.5,.75,.95,.99)),
  g_fa_near_bound = mean(r1$g_fa < gbnd | r1$g_fa > 1 - gbnd),
  n_domain = length(sub), A_prev = mean(obs$A[sub]), Y_prev = mean(obs$Y[sub]),
  n_psu = length(unique(paste(obs$strata[sub], obs$cluster[sub]))),
  n_strata = length(unique(obs$strata[sub])),
  design_df = summ$df[summ$method == "Fully-Aware-CF"][1],
  min_psu_stratum = min(table(unique(data.frame(s = obs$strata[sub], c = obs$cluster[sub]))$s)),
  weight_cv = sd(obs$weight[sub]) / mean(obs$weight[sub]))

out <- list(example = ex_id, label = label, rung = rung, learners = learners,
            B = B, summary = summ, per_split = allres, diagnostics = diagnostics)
fn <- file.path(OUT, sprintf("nh_%s_%s.rds", ex_id, rung))
saveRDS(out, fn)
cat(sprintf("[task %d] saved %s\n", task, fn)); print(summ[, c("method","b","se","lcl","ucl","df","b_split_sd")], row.names = FALSE)

# ---- reproducibility manifest ----------------------------------------------
manifest <- list(
  example = ex_id, label = label, rung = rung, learners = learners, B = B,
  n_domain = length(sub), covariates = covs, dummy_cols = wnames,
  R_version = R.version.string,
  packages = sapply(c("SuperLearner","tmle","survey","surveyCV","earth","gam","glmnet","ranger","mice"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA),
  sysname = Sys.info()[["nodename"]])
saveRDS(manifest, file.path(MAN, sprintf("manifest_%s_%s.rds", ex_id, rung)))
cat(sprintf("[task %d] DEFF(eif)=%.2f  g_FA=[%.3f,%.3f]  g_CF=[%.3f,%.3f]  near-bound(FA)=%.1f%%\n",
            task, dd$deff_clust, min(r1$g_fa), max(r1$g_fa), min(r1$g_cf), max(r1$g_cf),
            100 * diagnostics$g_fa_near_bound))
