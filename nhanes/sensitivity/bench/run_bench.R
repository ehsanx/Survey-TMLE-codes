# =====================================================================
# nhanes/sensitivity/bench/run_bench.R
#   Compute-time + variance-method benchmark for the NHANES application.
#
# GOAL. (1) Close the manuscript's stated limitation (Discussion: "we did not
# compare [the linearization] to the survey bootstrap or balanced repeated
# replication") -- it names BOTH the survey bootstrap AND BRR, so BOTH must
# appear. (2) Put a number behind "the linearization SE is ... cheap ... as
# opposed to computationally intensive bootstraps".
#
# WHAT. On ONE example (default E3, the smallest-n / fewest-PSU domain), under
# SINGLE imputation and the CURRENT engine, for the primary Fully-Aware-CF
# estimator, time + compare SIX variance routes and record diagnostics:
#   PANEL A -- design-based replicate SEs on the FROZEN fitted EIF (no refit):
#     1. LINEARIZATION (Eq-8 TSL)        -- our method: plug the fitted EIF.
#     2. delete-one-PSU JACKKNIFE (JKn)
#     3. balanced repeated replication (BRR, half-sample Hadamard)
#     4. Fay's BRR (rho = 0.3, the NCHS/NHANES default)
#     5. Rao-Wu / Rao-Wu-Yue survey bootstrap WEIGHTS (subbootstrap, R reps)
#   PANEL B -- full-refit resampling:
#     6. full-refit cluster BOOTSTRAP (B reps; re-draws PSUs within strata,
#        re-cross-fits + refits all SL nuisances).
# Panel A routes all consume the SAME fitted EIF (sub-second, no nuisance refit);
# Panel B is the expensive full-refit baseline. Also records: design df + deff,
# #strata/#PSU; B-convergence (SE at B={100,200,500}) + MC error + bias +
# percentile & BCa CIs from the retained B estimates; wall-clock, CPU time
# (user+sys), and core-seconds (the parallelism-invariant cost) per route; and
# the 5-arm suite wall-clock + machine configuration.
#
# NON-DESTRUCTIVE: SOURCES R/estimators.R + the R06_mi mi_helpers.R read-only;
# writes only under its own output dir.
#
# Env: SIM_CODE, NH_ANA, BENCH_EX (E3), BENCH_B (500), BENCH_SEED (20260607),
#      BENCH_LIB (override for a fast smoke), SLURM_CPUS_PER_TASK, BENCH_RESULTS.
# =====================================================================
t0 <- proc.time()[3]

CODE <- Sys.getenv("SIM_CODE", "R")
.this <- tryCatch({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)]); if (length(f)) normalizePath(f) else NA }, error = function(e) NA)
RUN_DIR <- if (!is.na(.this)) dirname(.this) else file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "sensitivity", "bench")
source(file.path(CODE, "estimators.R"))                              # run_estimators, .se_des, make_cf_folds
source(file.path(dirname(RUN_DIR), "R06_mi", "mi_helpers.R"))        # impute_m, encode_domain
suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(earth); library(glmnet); library(mice); library(parallel) })

EX    <- Sys.getenv("BENCH_EX", "E3")
B     <- { x <- Sys.getenv("BENCH_B"); if (nzchar(x)) as.integer(x) else 500L }
SEED  <- { x <- Sys.getenv("BENCH_SEED"); if (nzchar(x)) as.integer(x) else 20260607L }
cores <- { x <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(x)) as.integer(x) else max(1L, parallel::detectCores() - 1L) }
ANA   <- Sys.getenv("NH_ANA", file.path(Sys.getenv("REPO_ROOT","."), "nhanes", "analytic"))
RES   <- Sys.getenv("BENCH_RESULTS", file.path(Sys.getenv("REPO_ROOT","."), "results", "arc"))
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE, showWarnings = FALSE)
LIB   <- { x <- Sys.getenv("BENCH_LIB"); if (nzchar(x)) trimws(strsplit(x, ",")[[1]]) else c("SL.glm","SL.earth","SL.glmnet") }   # primary rung; override via BENCH_LIB for a fast smoke

# ---- machine configuration --------------------------------------------------
pk  <- sapply(c("SuperLearner","tmle","survey","surveyCV","earth","glmnet","mice"),
              function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA))
cpu <- tryCatch(paste(grep("Model name|^CPU\\(s\\)|Socket|Core", system("lscpu", intern = TRUE, ignore.stderr = TRUE), value = TRUE), collapse = " ; "), error = function(e) NA)
mem <- tryCatch(paste(system("free -h", intern = TRUE, ignore.stderr = TRUE)[1:2], collapse = " ; "), error = function(e) NA)

# ---- load example (single imputation) ---------------------------------------
obs0 <- readRDS(file.path(ANA, paste0(EX, "_analytic.rds")))
comp <- impute_m(obs0, M = 1L, seed = SEED)[[1]]
enc  <- encode_domain(comp)
ob   <- enc$obs; Wc <- enc$wnames
n_dom <- sum(ob$inpop)
cat(sprintf("[bench] EX=%s  n_domain=%d  B=%d  cores=%d  LIB={%s}\n", EX, n_dom, B, cores, paste(LIB, collapse=",")))

# ---- fit the 5 arms once (timed) --------------------------------------------
set.seed(SEED)
st_arms <- system.time(r <- run_estimators(ob, learners = LIB, V_cf = 5L, inner_cv_folds = 5L,
                                          W_cols = Wc, nest = TRUE, inpop = ob$inpop))
t_arms <- st_arms[["elapsed"]]
cf     <- r$results[r$results$method == "Fully-Aware-CF", ]
cf_eif <- r$diagnostics$arms[["Fully-Aware-CF"]]$eif1 - r$diagnostics$arms[["Fully-Aware-CF"]]$eif0   # CF ATE influence function (domain rows)

# ---- shared survey design (full design, subset to the domain) ---------------
full_e <- numeric(nrow(ob)); full_e[ob$inpop] <- cf_eif           # CF EIF on domain rows, 0 elsewhere
dd   <- data.frame(strata = ob$strata, cluster = ob$cluster, weight = ob$weight, e = full_e, dom = ob$inpop)
des  <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = dd, nest = TRUE)
dsub <- subset(des, dom)

# ---- design metadata (df, deff, #strata, #PSU) ------------------------------
df_design <- tryCatch(as.numeric(degf(dsub)), error = function(e) NA_real_)
deff_e    <- tryCatch(as.numeric(attr(svymean(~e, dsub, deff = TRUE), "deff")), error = function(e) NA_real_)
key       <- paste(ob$strata[ob$inpop], ob$cluster[ob$inpop], sep = "_")
n_strata  <- length(unique(ob$strata[ob$inpop])); n_psu <- length(unique(key))
tcrit     <- tryCatch(qt(0.975, df_design), error = function(e) qnorm(0.975))
ci_of     <- function(se) cf$b + c(-1, 1) * tcrit * se           # t_df Wald interval

# ---- helper: time a replicate-weight SE on the FROZEN EIF -------------------
rep_route <- function(type, ...) {
  st <- system.time({
    rep <- as.svrepdesign(dsub, type = type, ...)
    jk  <- withReplicates(rep, function(w, d) sum(w * d$e) / sum(w))   # Hajek mean of the EIF over the domain
    se  <- as.numeric(sqrt(attr(jk, "var")))
  })
  list(se = se, wall = st[["elapsed"]], cpu = st[["user.self"]] + st[["sys.self"]], cores = 1L)
}
safe_route <- function(type, ...) tryCatch(rep_route(type, ...), error = function(e) {
  cat(sprintf("[bench] route %s failed: %s\n", type, conditionMessage(e)))
  list(se = NA_real_, wall = NA_real_, cpu = NA_real_, cores = 1L) })

# ---- PANEL A: design-based replicate SEs on the frozen EIF (no refit) -------
se_lin  <- cf$se                                                  # authoritative CF linearization SE
st_lin  <- system.time(invisible(svymean(~e, dsub)))             # time the design-SE compute
A_lin <- list(se = se_lin, wall = st_lin[["elapsed"]], cpu = st_lin[["user.self"]] + st_lin[["sys.self"]], cores = 1L)
A_jkn <- safe_route("JKn")
A_brr <- safe_route("BRR")
A_fay <- safe_route("Fay", fay.rho = 0.3)
A_rw  <- safe_route("subbootstrap", replicates = 500L)           # Rao-Wu-Yue survey bootstrap weights

# ---- PANEL B: full-refit cluster bootstrap (re-cross-fits) ------------------
# Resample PSUs within strata WITH replacement; relabel each resampled PSU
# uniquely (so nest=TRUE + the CF folds treat duplicates as distinct PSUs);
# re-encode + refit. Returns (estimate, per-replicate elapsed seconds).
boot_one <- function(b) {
  set.seed(SEED + 7919L * b)
  cb  <- comp
  str <- cb$SDMVSTRA; psu <- cb$SDMVPSU
  rows <- integer(0); newpsu <- integer(0); nxt <- 0L
  for (s in unique(str)) {
    ps <- unique(psu[str == s]); mh <- length(ps)
    samp <- if (mh < 2L) ps else sample(ps, mh, replace = TRUE)
    for (p in samp) { ix <- which(str == s & psu == p); nxt <- nxt + 1L
                      rows <- c(rows, ix); newpsu <- c(newpsu, rep.int(nxt, length(ix))) }
  }
  cb <- cb[rows, , drop = FALSE]; cb$SDMVPSU <- newpsu
  attr(cb, "covs") <- attr(comp, "covs"); attr(cb, "example") <- attr(comp, "example")
  st <- system.time(rr <- tryCatch({
    e2 <- encode_domain(cb)
    rx <- run_estimators(e2$obs, learners = LIB, V_cf = 5L, inner_cv_folds = 5L,
                         W_cols = e2$wnames, nest = TRUE, inpop = e2$obs$inpop)
    rx$results$b[rx$results$method == "Fully-Aware-CF"]
  }, error = function(e) NA_real_))
  c(est = rr, secs = st[["elapsed"]])
}
st_boot <- system.time({
  ncl <- max(1L, min(cores, B))
  if (ncl > 1L) {
    cl <- makeCluster(ncl)
    invisible(clusterEvalQ(cl, { source(file.path(Sys.getenv("SIM_CODE","R"), "estimators.R"))
      suppressMessages({ library(survey); library(tmle); library(SuperLearner); library(earth); library(glmnet) }); TRUE }))
    clusterExport(cl, c("comp", "LIB", "SEED", "encode_domain"), envir = environment())
    bm <- parSapply(cl, seq_len(B), boot_one); stopCluster(cl)
  } else bm <- sapply(seq_len(B), boot_one)
})
bvals     <- as.numeric(bm["est", ]); boot_secs <- as.numeric(bm["secs", ])
n_boot_ok <- sum(!is.na(bvals))
se_boot   <- sd(bvals, na.rm = TRUE)
boot_wall <- st_boot[["elapsed"]]
boot_coresec <- sum(boot_secs, na.rm = TRUE)                      # total CPU work across replicates
ncl_boot  <- max(1L, min(cores, B))

# ---- bootstrap distribution diagnostics (FREE re-summaries of bvals) --------
bv     <- bvals[!is.na(bvals)]
seB    <- sapply(c(100L, 200L, 500L), function(k) if (length(bv) >= 2L) sd(bv[seq_len(min(k, length(bv)))]) else NA_real_)
mcse_seboot <- if (n_boot_ok > 1L) se_boot / sqrt(2 * (n_boot_ok - 1)) else NA_real_
boot_bias   <- mean(bv) - cf$b
pct_ci <- tryCatch(as.numeric(quantile(bv, c(.025, .975))), error = function(e) c(NA, NA))
# BCa: z0 from bvals, acceleration from the EIF (the empirical influence values)
z0  <- tryCatch(qnorm(mean(bv < cf$b)), error = function(e) NA_real_)
Lc  <- cf_eif - mean(cf_eif)
a_hat <- tryCatch(sum(Lc^3) / (6 * (sum(Lc^2))^1.5), error = function(e) NA_real_)
bca_q <- function(alpha) { z <- qnorm(alpha); pnorm(z0 + (z0 + z) / (1 - a_hat * (z0 + z))) }
bca_ci <- tryCatch(as.numeric(quantile(bv, c(bca_q(.025), bca_q(.975)))), error = function(e) c(NA, NA))

# ---- assemble S35 (six routes) ----------------------------------------------
mk <- function(name, panel, refit, nuis, route, ci = NULL, coresec = NULL) {
  cc <- if (is.null(coresec)) route$cpu else coresec
  ci <- if (is.null(ci)) ci_of(route$se) else ci
  data.frame(panel = panel, method = name, refit = refit, nuisance_unc = nuis,
             se = route$se, se_ratio = route$se / se_lin,
             ci_lo = ci[1], ci_hi = ci[2],
             wall_s = route$wall, cpu_s = route$cpu, cores = route$cores,
             coresec = cc, stringsAsFactors = FALSE)
}
boot_route <- list(se = se_boot, wall = boot_wall, cpu = boot_coresec, cores = ncl_boot)
S35 <- rbind(
  mk("Linearization (TSL)",            "A", "No",  "No",  A_lin),
  mk("Delete-one-PSU jackknife (JKn)", "A", "No",  "No",  A_jkn),
  mk("Balanced repeated replication",  "A", "No",  "No",  A_brr),
  mk("Fay's BRR (rho=0.3)",            "A", "No",  "No",  A_fay),
  mk("Rao-Wu bootstrap weights",       "A", "No",  "No",  A_rw),
  mk("Full-refit cluster bootstrap",   "B", "Yes", "Yes", boot_route, ci = pct_ci, coresec = boot_coresec)
)
lin_coresec   <- max(A_lin$cpu, 1e-4)
S35$coresec_ratio <- round(S35$coresec / lin_coresec, 1)
S35$se        <- round(S35$se, 5);    S35$se_ratio <- round(S35$se_ratio, 3)
S35$ci_lo     <- round(S35$ci_lo, 4); S35$ci_hi <- round(S35$ci_hi, 4)
S35$wall_s    <- round(S35$wall_s, 3); S35$cpu_s <- round(S35$cpu_s, 3); S35$coresec <- round(S35$coresec, 2)

diag <- data.frame(
  example = EX, n_domain = n_dom, n_strata = n_strata, n_psu = n_psu, df_design = df_design,
  deff = round(deff_e, 3), estimate_cf = cf$b, B = B, n_boot_ok = n_boot_ok,
  se_lin = round(se_lin, 5), se_boot = round(se_boot, 5),
  seB100 = round(seB[1], 5), seB200 = round(seB[2], 5), seB500 = round(seB[3], 5),
  mcse_seboot = round(mcse_seboot, 6), boot_bias = round(boot_bias, 6),
  pct_lo = round(pct_ci[1], 4), pct_hi = round(pct_ci[2], 4),
  bca_lo = round(bca_ci[1], 4), bca_hi = round(bca_ci[2], 4),
  stringsAsFactors = FALSE)

# ---- write outputs ----------------------------------------------------------
write.csv(S35,  file.path(RES, sprintf("bench_variance_%s.csv", EX)), row.names = FALSE)
write.csv(diag, file.path(RES, sprintf("bench_diagnostics_%s.csv", EX)), row.names = FALSE)
writeLines(c("== machine configuration ==",
             paste("R:", R.version.string), paste("node:", Sys.info()[["nodename"]]),
             paste("os:", Sys.info()[["sysname"]], Sys.info()[["release"]]),
             paste("cores_used:", cores), paste("cpu:", cpu), paste("mem:", mem),
             paste("RNGkind:", paste(RNGkind(), collapse = " ")), paste("master_seed:", SEED),
             sprintf("E3_design: strata=%d PSU=%d df=%s deff=%.3f", n_strata, n_psu, df_design, deff_e),
             sprintf("bootstrap: B=%d reps, PSUs resampled within strata w/ replacement, re-cross-fits + refits all nuisances; parallelized over %d cores", B, ncl_boot),
             sprintf("SL_library(%d): %s", length(LIB), paste(LIB, collapse = ", ")),
             "== package versions ==", paste(names(pk), pk, sep = ": ")),
           file.path(RES, "bench_machine_config.txt"))
saveRDS(list(S35 = S35, diag = diag, five_arm_seconds = t_arms, five_arm = r$results,
             cf_estimate = cf$b, cf_eif_n = length(cf_eif), boot_vals = bvals, boot_secs = boot_secs,
             machine = readLines(file.path(RES, "bench_machine_config.txt"))),
        file.path(RES, sprintf("bench_%s.rds", EX)))

cat("\n==== VARIANCE-METHOD BENCHMARK (", EX, ", single imputation) ====\n", sep = "")
print(S35[, c("panel","method","se","se_ratio","ci_lo","ci_hi","wall_s","coresec","coresec_ratio")], row.names = FALSE)
cat(sprintf("\nE3 design: strata=%d, PSU=%d, df=%g, deff=%.3f | CF estimate %.4f | 5-arm suite %.1fs\n",
            n_strata, n_psu, df_design, deff_e, cf$b, t_arms))
cat(sprintf("B-convergence SE@{100,200,500} = {%.5f, %.5f, %.5f}; MC-SE(boot SE)=%.5g; bias=%.5g\n",
            seB[1], seB[2], seB[3], mcse_seboot, boot_bias))
cat(sprintf("percentile CI=(%.4f, %.4f); BCa CI=(%.4f, %.4f)\n", pct_ci[1], pct_ci[2], bca_ci[1], bca_ci[2]))
cat(sprintf("[bench] done in %.1f min.\n", (proc.time()[3] - t0) / 60))
