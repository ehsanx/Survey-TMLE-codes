# =====================================================================
# R13_null_typeI / aggregate.R
#   Combine the per-task RDS from run.R into ONE summary CSV mirroring
#   aggregate_sim.R's column conventions EXACTLY, with the run's knob
#   column `te` leading and two new columns: reject_rate (test of
#   H0: psi = 0) + its MCSE. Writes results/arc/R13_null_typeI_summary.csv
#   (never the locked results/sim_full_summary.csv).
#
# Conventions copied verbatim from aggregate_sim.R so numbers are
# directly comparable to the headline table:
#   crit        = qt(0.975, max(1, df))                  (t reference)
#   coverage    = mean(|b - Psi| <= crit*se)             (true-Psi coverage)
#   bias        = mean(b) - Psi ; emp_sd = sd(b)
#   se_ratio    = mean(se)/sd(b) ; mcse_cov = sqrt(cov*(1-cov)/n_reps)
# New (this run):
#   reject_rate = mean(|b|/se > crit)                    (test of H0: psi=0)
#   mcse_reject = sqrt(rr*(1-rr)/n_reps)                 (empirical MCSE)
#
# Decision views printed:
#   (1) type-I table at te=0 — flag any method/rung outside 0.05 +/- 2*MCSE,
#       where the flag band uses the NOMINAL-rate MCSE sqrt(.05*.95/n)
#       (the empirical MCSE degenerates to 0 when rr hits 0 or 1);
#       the Fully-Aware-CF rows are the audit target. A FLAG must be read
#       JOINTLY with bias + se_ratio (printed guidance below): CF carries the
#       known nuisance-misspecification bias at L1-L3 under the locked complex
#       DGP, and its conservative SEs can push the rate below the band.
#   (2) power table at te=0.3 (FA/CF/CV rows are the interpretable power
#       numbers; PA/NA "power" rides on invalid SEs).
#
# Usage:
#   Rscript codes/arc_runs/R13_null_typeI/aggregate.R
#   (reads R13_OUT; skips SMOKE_* files and the manifest/ subdir;
#    writes results/arc/R13_null_typeI_summary.csv +
#    R13_OUT/R13_null_typeI_combined.rds)
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from the repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

RUN_ID <- "R13_null_typeI"
OUT <- Sys.getenv("R13_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RESULTS)) dir.create(ARC_RESULTS, recursive = TRUE, showWarnings = FALSE)

# ^r13_ skips SMOKE_* outputs; manifests live in the manifest/ subdir (not listed)
files <- list.files(OUT, pattern = "^r13_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no r13_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
# o$results already carries the per-rep knob column `te` (tagged in run.R)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))

# null-truth audit across te=0 tasks (recorded per task in run.R; warn-not-stop)
nulltab <- unique(do.call(rbind, lapply(objs, function(o)
  data.frame(te = o$te, scenario = o$scenario, Psi = o$Psi,
             se_mc = o$truth$se_mc, null_ok = o$null_ok))))
cat("\nnull-truth audit (per population; null_ok = |psi| <= 4*se_mc at te=0):\n")
print(nulltab[order(nulltab$te, nulltab$scenario), ], row.names = FALSE)
if (any(nulltab$te == 0 & !vapply(nulltab$null_ok, isTRUE, logical(1))))
  cat("WARNING: a te=0 population failed the null truth check (see table above)\n")

z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows, list(rows$te, rows$scenario, rows$rung, rows$method),
  function(d) {
    if (is.null(d) || !nrow(d)) return(NULL)
    Psi <- d$Psi[1]; crit <- z_or_t(d$df)
    cov <- mean(abs(d$b - Psi) <= crit * d$se)
    rej <- mean(abs(d$b) / d$se > crit)          # test of H0: psi = 0
    data.frame(te = d$te[1], scenario = d$scenario[1], rung = d$rung[1],
               method = d$method[1], n_reps = nrow(d), Psi = Psi,
               bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
               se_ratio = mean(d$se) / sd(d$b),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
               reject_rate = rej, mcse_reject = sqrt(rej * (1 - rej) / nrow(d)),
               deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif))
  }))
rownames(agg) <- NULL
ord_m   <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
rung_lv <- c("L1_param", "L2_smooth", "L3_adaptive", "L4_aggressive")
agg <- agg[order(agg$te, agg$scenario, match(agg$rung, rung_lv),
                 match(agg$method, ord_m)), ]

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
cat("\nfull summary (te, scenario, rung, method):\n")
print(agg_print, row.names = FALSE)

# ---- DECISION VIEW 1: type-I error at the null (te = 0; target 0.05) --------
# Flag band: 0.05 +/- 2 * sqrt(0.05*0.95/n_reps) — the MCSE of the rejection
# rate UNDER nominal 0.05 (the empirical mcse_reject is reported in the CSV but
# degenerates to 0 when rr = 0 or 1, which would flag everything/nothing).
cat("\n=== DECISION VIEW 1: type-I error at te=0 (test of H0: psi=0; target 0.05) ===\n")
t1 <- agg[agg$te == 0, c("scenario", "rung", "method", "n_reps",
                         "reject_rate", "mcse_reject")]
mcse_nom <- sqrt(0.05 * 0.95 / t1$n_reps)
t1$lo_2mcse <- 0.05 - 2 * mcse_nom
t1$hi_2mcse <- 0.05 + 2 * mcse_nom
t1$flag  <- ifelse(t1$reject_rate < t1$lo_2mcse | t1$reject_rate > t1$hi_2mcse,
                   "FLAG", "ok")
t1$audit <- ifelse(t1$method == "Fully-Aware-CF", "<= AUDIT TARGET", "")
t1p <- t1; nn <- sapply(t1p, is.numeric); t1p[nn] <- round(t1p[nn], 4)
print(t1p, row.names = FALSE)
cf <- t1[t1$method == "Fully-Aware-CF", ]
cat(sprintf("\nFully-Aware-CF type-I range: [%.4f, %.4f] across %d cells; %d flagged.\n",
            min(cf$reject_rate), max(cf$reject_rate), nrow(cf),
            sum(cf$flag == "FLAG")))
cat(sprintf("All methods: %d of %d (te=0) rows flagged outside 0.05 +/- 2*MCSE.\n",
            sum(t1$flag == "FLAG"), nrow(t1)))
cat("Expected: CF ~0.05 everywhere; PA/NA over-reject (invalid SEs);\n",
    "single-fit FA over-rejects at L4 (the Donsker failure).\n", sep = "")
cat("HOW TO READ A CF FLAG (jointly with the bias + se_ratio columns):\n",
    "under the locked model_type='complex' DGP the CF arm carries the known\n",
    "nuisance-misspecification bias at L1-L3 (~ +0.018 at the null, matching the\n",
    "headline bias rows), so the te=0 rejection of H0: psi=0 conflates bias with\n",
    "SE calibration. CF rows can leave the band UPWARD from that bias (L1-L3) or\n",
    "DOWNWARD from the conservative CF SEs (se_ratio ~1.4 => effective critical\n",
    "value ~2.8 empirical-SD units). A FLAG is not by itself a type-I failure:\n",
    "R01 simple-control is the matched evidence separating misspecification from\n",
    "variance miscalibration; the cleanest audit-target cell is CF at L4 (bias\n",
    "~0.002-0.006).\n", sep = "")

# ---- DECISION VIEW 2: power at te = 0.3 --------------------------------------
cat("\n=== DECISION VIEW 2: power at te=0.3 (reject H0: psi=0) ===\n")
t2 <- agg[agg$te == 0.3, c("scenario", "rung", "method", "n_reps", "Psi",
                           "reject_rate", "mcse_reject", "bias", "coverage")]
t2p <- t2; nn <- sapply(t2p, is.numeric); t2p[nn] <- round(t2p[nn], 4)
print(t2p, row.names = FALSE)
cat("NOTE: only the design-valid arms (FA / CF / CV) give honest power;\n",
    "PA/NA rejection rates ride on invalid SEs and are shown for completeness.\n", sep = "")

# ---- outputs ------------------------------------------------------------------
out_csv <- file.path(ARC_RESULTS, sprintf("%s_summary.csv", RUN_ID))
write.csv(agg_print, out_csv, row.names = FALSE)
saveRDS(list(per_rep = rows, summary = agg, null_audit = nulltab,
             typeI_view = t1, power_view = t2),
        file.path(OUT, sprintf("%s_combined.rds", RUN_ID)))
cat("\nSaved:\n -", out_csv,
    "\n -", file.path(OUT, sprintf("%s_combined.rds", RUN_ID)), "\n")
