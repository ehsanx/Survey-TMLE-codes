# =====================================================================
# aggregate.R  —  R15_aipw_benchmark: combine per-task RDS into the summary
# (spec item A7, Writing/comments/phase4-arc-sim-specs.md)
#
# Usage:  Rscript codes/arc_runs/R15_aipw_benchmark/aggregate.R
#         (env: R15_OUT to point at the per-task RDS dir; SIM_CODE optional)
#
# Reads  $R15_OUT/r15_<scenario>_<rung>_chunk###.rds  (skips SMOKE_* files;
# manifests live in manifest/ and are not matched). Coverage is computed
# TWICE per arm -- once with the JKn replicate SE, once with the Eq-8
# linearized SE -- both against a t(df) reference, using aggregate_sim.R's
# EXACT formulas: crit = qt(.975, pmax(1, df)); coverage =
# mean(|b - Psi| <= crit*se); bias = mean(b) - Psi; emp_sd = sd(b);
# se_ratio = mean(se)/emp_sd; mcse_cov = sqrt(cov*(1-cov)/n_reps).
#
# Writes results/arc/R15_aipw_benchmark_summary.csv  (long: se_type leading)
#        $R15_OUT/R15_aipw_benchmark_combined.rds    (per-rep + summary + diag)
# Prints the decision views:
#   (1) AIPW arms beside the LOCKED five paper arms (results/sim_full_summary.csv)
#       per (scenario, rung) -- bias / emp_sd / mean_se / se_ratio / coverage;
#   (2) the L4 two-fold-splitting SE-cost comparison (locked CF SE/SD ~ 1.44
#       at L4 vs the AIPW arms);
#   (3) JKn vs linearization SE agreement per cell.
# =====================================================================

# locate config.R robustly: SIM_CODE env, else repo-relative, else Windows default
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from the repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R15_OUT", file.path(DATA_ROOT, "arc_runs", "R15_aipw_benchmark"))
ARC_RESULTS <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RESULTS)) dir.create(ARC_RESULTS, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^r15_.*_chunk[0-9]{3}\\.rds$", full.names = TRUE)
if (!length(files)) stop("no r15_*_chunk*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))
diags <- do.call(rbind, lapply(objs, function(o)
  if (!is.null(o$diagnostics)) cbind(scenario = o$scenario, rung = o$rung, o$diagnostics)))

# ---- integrity guards (chunk-aliasing / partial-failure protection) ---------
# (1) HARD STOP on any duplicated (scenario, rung, method, rep): duplicates can
#     only arise from chunk files written under DIFFERENT SIM_CHUNK values
#     aliasing each other's chunk numbers (see NOTES.md "Walltime recovery") --
#     rbind would silently double-count those reps in every summary.
stopifnot(!anyDuplicated(rows[, c("scenario", "rung", "method", "rep")]))
# (2) LOUD per-(scenario, rung) completeness check: WARN (not stop) when the
#     distinct rep count != the expected N. Expected N = R15_EXPECTED_REPS env
#     (default 1000), raised to max(rep) seen so over-runs also flag.
exp_reps <- as.integer(Sys.getenv("R15_EXPECTED_REPS", "1000"))
exp_reps <- max(exp_reps, max(rows$rep))
n_incomplete <- 0L
for (cell in split(rows, list(rows$scenario, rows$rung), drop = TRUE)) {
  n_dist <- length(unique(cell$rep))
  if (n_dist != exp_reps) {
    n_incomplete <- n_incomplete + 1L
    cat(sprintf("**** WARN [R15 completeness] %s x %s: %d distinct reps != expected %d -- missing/failed chunks or reps; summaries below are based on the incomplete cell ****\n",
                cell$scenario[1], cell$rung[1], n_dist, exp_reps))
  }
}
if (n_incomplete == 0L) {
  cat(sprintf("[R15 completeness] OK: all %d (scenario x rung) cells have %d distinct reps\n",
              length(unique(paste(rows$scenario, rows$rung))), exp_reps))
} else {
  cat(sprintf("**** WARN [R15 completeness] %d incomplete cell(s) -- see lines above (set R15_EXPECTED_REPS to silence if intentional) ****\n",
              n_incomplete))
}

# ---- summary: per (scenario, rung, method) x se_type (jkn, lin) -------------
# bias_med is reported alongside the conventional mean bias: AIPW has no
# targeting step, so occasional extreme nuisance fits pass straight into the
# pseudo-outcome and can fatten mean bias / emp_sd (collapse shares tracked
# below; expected ~0 after the helper's mean-1 weight normalization).
agg_one <- function(d, se, se_type) {
  Psi  <- d$Psi[1]
  crit <- qt(0.975, pmax(1, d$df))
  cov  <- mean(abs(d$b - Psi) <= crit * se)
  data.frame(se_type = se_type,
             scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
             n_reps = nrow(d), Psi = Psi,
             bias = mean(d$b) - Psi, bias_med = median(d$b) - Psi,
             emp_sd = sd(d$b), mean_se = mean(se),
             se_ratio = mean(se) / sd(d$b),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
             deff_clust = mean(d$deff), icc_D = mean(d$icc_D))
}
agg <- do.call(rbind, by(rows, list(rows$scenario, rows$rung, rows$method), function(d) {
  if (is.null(d)) return(NULL)
  rbind(agg_one(d, d$se_jkn, "jkn"), agg_one(d, d$se_lin, "lin"))
}))
rownames(agg) <- NULL
ord_m <- c("AIPW-SF", "AIPW-CF")
agg <- agg[order(agg$scenario, agg$rung, match(agg$method, ord_m), agg$se_type), ]

num <- sapply(agg, is.numeric); agg_print <- agg; agg_print[num] <- round(agg_print[num], 4)
cat("\n==== R15 AIPW benchmark summary (coverage TWICE: JKn + Eq-8 lin) ====\n")
print(agg_print, row.names = FALSE)

# ---- nuisance-collapse diagnostic per cell: counts reps whose raw g is
# ENTIRELY outside the floor (genuine weighted quasi-separation; expected 0
# after the helper's mean-1 weight normalization -- investigate otherwise) ----
diag_sum <- NULL
if (!is.null(diags)) {
  diag_sum <- do.call(rbind, by(diags, list(diags$scenario, diags$rung), function(d) {
    if (is.null(d)) return(NULL)
    data.frame(scenario = d$scenario[1], rung = d$rung[1], n_reps = nrow(d),
               n_gsf_collapsed = sum(d$g_sf_outside_floor >= 1),
               n_gcf_collapsed = sum(d$g_cf_outside_floor >= 1),
               g_sf_outside_floor = round(mean(d$g_sf_outside_floor), 4),
               g_cf_outside_floor = round(mean(d$g_cf_outside_floor), 4),
               cf_V_eff = d$cf_V_eff[1], deff = round(mean(d$deff), 4))
  }))
  rownames(diag_sum) <- NULL
  diag_sum <- diag_sum[order(diag_sum$scenario, diag_sum$rung), ]
  cat("\n==== nuisance-collapse diagnostics (per cell) ====\n")
  print(diag_sum, row.names = FALSE)
}

# ---- save -----------------------------------------------------------------
saveRDS(list(per_rep = rows, summary = agg, per_rep_diag = diags,
             diag_summary = diag_sum),
        file.path(OUT, "R15_aipw_benchmark_combined.rds"))
write.csv(agg_print, file.path(ARC_RESULTS, "R15_aipw_benchmark_summary.csv"),
          row.names = FALSE)

# ============================================================================
# DECISION VIEW 1: AIPW vs the LOCKED five paper arms, per (scenario, rung)
# ============================================================================
locked_fn <- file.path(RESULTS_DIR, "sim_full_summary.csv")
if (file.exists(locked_fn)) {
  lk <- read.csv(locked_fn, stringsAsFactors = FALSE)
  lk_v <- data.frame(scenario = lk$scenario, rung = lk$rung,
                     method = lk$method, se_type = "lin(Eq8)",
                     bias = lk$bias, emp_sd = lk$emp_sd, mean_se = lk$mean_se,
                     se_ratio = lk$se_ratio, coverage = lk$coverage,
                     stringsAsFactors = FALSE)
  ai_v <- data.frame(scenario = agg$scenario, rung = agg$rung,
                     method = agg$method, se_type = agg$se_type,
                     bias = round(agg$bias, 4), emp_sd = round(agg$emp_sd, 4),
                     mean_se = round(agg$mean_se, 4),
                     se_ratio = round(agg$se_ratio, 4),
                     coverage = round(agg$coverage, 4), stringsAsFactors = FALSE)
  ord_all <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV",
               "Partially-Aware", "Non-Aware", "AIPW-SF", "AIPW-CF")
  cat("\n==== DECISION VIEW 1: AIPW beside the locked five arms ====\n")
  for (sc in unique(ai_v$scenario)) for (rg in unique(ai_v$rung[ai_v$scenario == sc])) {
    v <- rbind(lk_v[lk_v$scenario == sc & lk_v$rung == rg, ],
               ai_v[ai_v$scenario == sc & ai_v$rung == rg, ])
    v <- v[order(match(v$method, ord_all)), ]
    cat(sprintf("\n-- %s x %s --\n", sc, rg))
    print(v[, c("method", "se_type", "bias", "emp_sd", "mean_se", "se_ratio", "coverage")],
          row.names = FALSE)
  }

  # ==========================================================================
  # DECISION VIEW 2: the L4 fold-splitting SE cost ("compared to what?")
  # Locked Fully-Aware-CF pays SE/SD ~ 1.44 at L4 (the cost of PSU-level
  # cross-fitting + pooled targeting). Does the AIPW competitor pay it too?
  # ==========================================================================
  cat("\n==== DECISION VIEW 2: L4 SE/SD (fold-splitting SE cost) ====\n")
  l4 <- rbind(
    lk_v[lk_v$rung == "L4_aggressive" & lk_v$method %in% c("Fully-Aware", "Fully-Aware-CF"), ],
    ai_v[ai_v$rung == "L4_aggressive", ])
  l4 <- l4[order(l4$scenario, match(l4$method, ord_all), l4$se_type), ]
  print(l4[, c("scenario", "method", "se_type", "bias", "se_ratio", "coverage")],
        row.names = FALSE)
  cat("(locked Fully-Aware-CF se_ratio at L4: standard 1.4389, R1 1.4403)\n")
} else {
  cat("\n[note] locked", locked_fn, "not found -- skipping decision views 1-2\n")
}

# ============================================================================
# DECISION VIEW 3: JKn vs Eq-8 linearization agreement (per cell x method)
# ============================================================================
cat("\n==== DECISION VIEW 3: mean se_jkn / mean se_lin per cell ====\n")
jk <- agg[agg$se_type == "jkn", c("scenario", "rung", "method", "mean_se")]
ln <- agg[agg$se_type == "lin", c("scenario", "rung", "method", "mean_se")]
names(jk)[4] <- "mean_se_jkn"; names(ln)[4] <- "mean_se_lin"
v3 <- merge(jk, ln, by = c("scenario", "rung", "method"))
v3$jkn_over_lin <- round(v3$mean_se_jkn / v3$mean_se_lin, 4)
v3$mean_se_jkn <- round(v3$mean_se_jkn, 4); v3$mean_se_lin <- round(v3$mean_se_lin, 4)
v3 <- v3[order(v3$scenario, v3$rung, match(v3$method, ord_m)), ]
print(v3, row.names = FALSE)

cat("\nSaved:\n -", file.path(ARC_RESULTS, "R15_aipw_benchmark_summary.csv"),
    "\n -", file.path(OUT, "R15_aipw_benchmark_combined.rds"), "\n")
