# =====================================================================
# aggregate.R  —  R21_deployable_cvcf
# Pools the per-task RDS into a per-(scenario, method) summary at the
# DEPLOYABLE non-Donsker library L6, and prints the headline 3x2 DECISION
# VIEW: {single-fit, cluster-CV, CF} x {TMLE, AIPW} with BOTH coverage AND
# se_ratio. Uses aggregate_sim.R's EXACT formulas (qt(.975, df)).
#
# Coverage uses the unified `se` column (TMLE: tmle EIF SE; AIPW: the Eq-8
# linearized SE), so the TMLE-vs-AIPW cells are apples-to-apples (both
# EIF/linearization-based). The se_ratio column is MANDATORY: the honest claim
# is ONE-SIDED -- "single-fit/cluster-CV under-cover; CF does NOT under-cover"
# -- and CF is CONSERVATIVE (se_ratio>1), not exactly nominal, so se_ratio must
# be shown to avoid over-claiming CF as "nominal".
#
# Env: SIM_CODE, REPO_ROOT, DATA_ROOT, R21_OUT (per-task RDS dir),
#      RESULTS_DIR (summary CSV dir; via config.R = REPO_ROOT/results or the
#      ARC bundle's output/results_arc).
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "R")
if (!file.exists(file.path(CODE, "config.R"))) CODE <- "R"
source(file.path(CODE, "config.R"))

RUN_ID <- "R21_deployable_cvcf"
OUT  <- Sys.getenv("R21_OUT", file.path(DATA_ROOT, "arc_runs", RUN_ID))
ARCR <- file.path(RESULTS_DIR, "arc"); if (!dir.exists(ARCR)) dir.create(ARCR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(OUT, pattern = "^r21(hal)?_.*\\.rds$", full.names = TRUE)  # L6 (r21_) + L7_hal (r21hal_); skips SMOKE_, manifest/, *_combined
files <- files[!grepl("_combined\\.rds$", files)]
if (!length(files)) stop("no r21_*.rds in ", OUT, " -- run the array first")

per_rep <- do.call(rbind, lapply(files, function(f) {
  o <- readRDS(f)
  cbind(scenario = o$scenario, Psi = o$Psi, o$results)   # o$results carries rung (L6_deployable / L7_hal) + base_m
}))

# ---- integrity guard: no duplicated (scenario, base_m, method, rep) ----------
key <- per_rep[, c("scenario", "base_m", "rung", "method", "rep")]
stopifnot(!anyDuplicated(key))

# exact aggregate_sim.R math, split by (scenario, base_m, method)
bm_grp <- ifelse(is.na(per_rep$base_m), "_default", as.character(per_rep$base_m))  # NA-safe group key for the default (non-sweep) run
summ <- do.call(rbind, by(per_rep, list(per_rep$scenario, bm_grp, per_rep$rung, per_rep$method),
                          function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  Psi <- d$Psi[1]; df <- pmax(1, d$df); crit <- qt(0.975, df)
  cov <- mean(abs(d$b - Psi) <= crit * d$se); n <- nrow(d)
  data.frame(scenario = d$scenario[1], rung = d$rung[1],
             base_m = d$base_m[1], estimator = d$estimator[1], method = d$method[1],
             n_reps = n, Psi = round(Psi, 4), bias = round(mean(d$b) - Psi, 4),
             emp_sd = round(sd(d$b), 4), mean_se = round(mean(d$se), 4),
             se_ratio = round(mean(d$se) / sd(d$b), 4), coverage = round(cov, 4),
             mcse_cov = round(sqrt(cov * (1 - cov) / n), 4),
             deff_clust = round(mean(d$deff), 4), icc_eif = round(mean(d$icc_eif), 4),
             stringsAsFactors = FALSE)
}))
rownames(summ) <- NULL
ord_m <- c("Fully-Aware", "Fully-Aware-CV", "Fully-Aware-CF",
           "AIPW-SF", "AIPW-CV", "AIPW-CF")
summ <- summ[order(summ$scenario, summ$base_m, summ$rung, match(summ$method, ord_m)), ]

csv <- file.path(ARCR, paste0(RUN_ID, "_summary.csv"))
write.csv(summ, csv, row.names = FALSE)
saveRDS(list(summary = summ, per_rep = per_rep), file.path(OUT, paste0(RUN_ID, "_combined.rds")))
cat("wrote", csv, "\n\n=== FULL SUMMARY (L6 = glm + earth + glmnet + tuned-xgb + deep RF) ===\n")
print(summ, row.names = FALSE)

# ---- HEADLINE: 3x2 DECISION VIEW -------------------------------------------
# {single-fit, cluster-CV, CF} x {TMLE, AIPW}, coverage AND se_ratio.
map3 <- list(
  TMLE = c("single-fit" = "Fully-Aware", "cluster-CV" = "Fully-Aware-CV",
           "CF"          = "Fully-Aware-CF"),
  AIPW = c("single-fit" = "AIPW-SF",     "cluster-CV" = "AIPW-CV",
           "CF"          = "AIPW-CF"))
designs <- c("single-fit", "cluster-CV", "CF")

cat("\n=== 3x2 DECISION VIEW: {single-fit, cluster-CV, CF} x {TMLE, AIPW} at L6 ===\n")
cat("Question: at a DEPLOYABLE non-Donsker library, is cluster-CV != CF\n",
    "ESTIMATOR-AGNOSTIC? i.e. do BOTH TMLE-CV and AIPW-CV under-cover where the\n",
    "single-fit arms do, while only CF holds? (CF is conservative, se_ratio>1.)\n", sep = "")

cells_sd <- unique(summ[, c("scenario", "base_m", "rung")])
for (k in seq_len(nrow(cells_sd))) {
  sc <- cells_sd$scenario[k]; bm <- cells_sd$base_m[k]; rg <- cells_sd$rung[k]
  blab <- if (is.na(bm)) "" else sprintf(" (base_m=%d)", bm)
  sub <- summ[summ$scenario == sc & summ$rung == rg & ((is.na(summ$base_m) & is.na(bm)) | summ$base_m %in% bm), ]
  get <- function(method, col) { r <- sub[sub$method == method, col]; if (length(r)) r[1] else NA }
  cat(sprintf("\n-- %s | %s%s --\n", rg, sc, blab))
  cat(sprintf("   %-12s | %-26s | %-26s\n", "design", "TMLE (cov / se_ratio)", "AIPW (cov / se_ratio)"))
  cat(sprintf("   %-12s-+-%-26s-+-%-26s\n", strrep("-", 12), strrep("-", 26), strrep("-", 26)))
  for (dz in designs) {
    tm <- map3$TMLE[[dz]]; ai <- map3$AIPW[[dz]]
    tcov <- get(tm, "coverage"); tsr <- get(tm, "se_ratio")
    acov <- get(ai, "coverage"); asr <- get(ai, "se_ratio")
    cat(sprintf("   %-12s | cov=%-6s se_ratio=%-7s | cov=%-6s se_ratio=%-7s\n",
                dz,
                ifelse(is.na(tcov), "NA", sprintf("%.3f", tcov)),
                ifelse(is.na(tsr),  "NA", sprintf("%.3f", tsr)),
                ifelse(is.na(acov), "NA", sprintf("%.3f", acov)),
                ifelse(is.na(asr),  "NA", sprintf("%.3f", asr))))
  }
  # per-cell verdict (one-sided): CF holds (>=0.93) AND both CV arms under-cover (<0.90)?
  tcf <- get("Fully-Aware-CF", "coverage"); acf <- get("AIPW-CF", "coverage")
  tcv <- get("Fully-Aware-CV", "coverage"); acv <- get("AIPW-CV", "coverage")
  if (all(is.finite(c(tcf, acf, tcv, acv)))) {
    verdict <-
      if (tcf >= 0.93 && acf >= 0.93 && tcv < 0.90 && acv < 0.90)
        "BOTH CF arms hold; BOTH cluster-CV arms under-cover -> cluster-CV != CF is ESTIMATOR-AGNOSTIC on a deployable class (the C3 cell; strengthens A1)."
      else if (tcf >= 0.93 && acf >= 0.93 && tcv >= 0.93 && acv >= 0.93)
        "BOTH CF and BOTH CV arms hold -> internal CV's down-weighting recovered coverage on L6 (HEADLINE FLIPS; reportable: CV retreats to smooth members, CF keeps the flexible learner with a guarantee)."
      else if (tcf >= 0.93 && acf >= 0.93 && acv >= 0.93 && tcv < 0.90)
        "CF holds (both); AIPW-CV nominal but TMLE-CV under-covers -> CV pathology LOCALIZES to TMLE; A1's framing WEAKENS (the gating outcome flagged in the spec)."
      else "mixed pattern -> inspect bias/se_ratio jointly per arm before quoting."
    cat(sprintf("   => %s\n", verdict))
  }
}

# ---- locked-target context (graceful skip) ----------------------------------
cat("\n=== LOCKED TARGETS for reference ===\n")
cat(" R20/L5 TMLE: standard CF 0.948 (se_ratio 1.066) / single-fit 0.905 / cluster-CV 0.850\n")
cat(" R15/L4 AIPW: AIPW-SF 0.621 (se_ratio ~0.43) / AIPW-CF 0.96\n")
cat(" Expect at L6: TMLE-SF ~0.90, TMLE-CF ~0.93-0.95, AIPW-SF under-covers, AIPW-CF ~0.96;\n")
cat("   HYPOTHESIS: BOTH cluster-CV arms ALSO under-cover (estimator-agnostic).\n")

locked <- file.path(RESULTS_DIR, "sim_full_summary.csv")
if (file.exists(locked)) {
  L <- read.csv(locked, stringsAsFactors = FALSE)
  ctx <- L[L$method %in% c("Fully-Aware-CV", "Fully-Aware-CF", "Fully-Aware") &
           L$rung %in% c("L3_adaptive", "L4_aggressive"),
           c("scenario", "rung", "method", "coverage", "se_ratio")]
  if (nrow(ctx)) {
    cat("\n=== LADDER CONTEXT (locked headline run) ===\n")
    print(ctx[order(ctx$scenario, ctx$rung, ctx$method), ], row.names = FALSE)
  }
} else {
  cat("\n[note] locked results/sim_full_summary.csv not found under RESULTS_DIR; ladder-context view skipped.\n")
}

cat("\nSaved:\n -", csv, "\n -", file.path(OUT, paste0(RUN_ID, "_combined.rds")), "\n")
