# =====================================================================
# aggregate.R  —  combine R18_cvu per-task RDS into the CV-u summary + the
#                 DECISION VIEW: CVu coverage vs the locked CV-w / CF / FA rows
# Usage:  Rscript codes/arc_runs/R18_cvu/aggregate.R   (after the array job)
#
# Mirrors codes/aggregate_sim.R conventions EXACTLY (same coverage / se_ratio /
# MCSE math, same column names) so CVu rows drop next to the locked table.
# Reads the LOCKED results/sim_full_summary.csv (read-only) for the comparator
# rows; writes results/arc/R18_cvu_summary.csv -- never clobbers locked files.
# Skips SMOKE_* files automatically (pattern anchored at ^r18_) and the
# manifest/ subdir (list.files is non-recursive).
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R18_OUT", file.path(DATA_ROOT, "arc_runs", "R18_cvu"))
files <- list.files(OUT, pattern = "^r18_.*\\.rds$", full.names = TRUE)  # excludes SMOKE_* + manifest/
if (!length(files)) stop("no r18_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))
diags <- do.call(rbind, lapply(objs, function(o)
  if (!is.null(o$diagnostics)) cbind(scenario = o$scenario, rung = o$rung, o$diagnostics)))

# ---- summary: EXACT aggregate_sim.R formulas ---------------------------------
z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows, list(rows$scenario, rows$rung, rows$method), function(d) {
  if (is.null(d)) return(NULL)
  Psi <- d$Psi[1]; crit <- z_or_t(d$df)
  cov <- mean(abs(d$b - Psi) <= crit * d$se)
  data.frame(cv_nuis_wts = "unweighted",          # the run's knob (vs locked CV-w "weighted")
             scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
             n_reps = nrow(d), Psi = Psi,
             bias = mean(d$b) - Psi, emp_sd = sd(d$b), mean_se = mean(d$se),
             se_ratio = mean(d$se) / sd(d$b),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / nrow(d)),
             deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif))
}))
rownames(agg) <- NULL
rung_lv <- c("L2_smooth", "L3_adaptive")
agg <- agg[order(agg$scenario, match(agg$rung, rung_lv)), ]

num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
cat("\n==== R18_cvu summary (new Fully-Aware-CVu rows) ====\n")
print(ap, row.names = FALSE)

dir.create(file.path(RESULTS_DIR, "arc"), recursive = TRUE, showWarnings = FALSE)
csv <- file.path(RESULTS_DIR, "arc", "R18_cvu_summary.csv")
write.csv(ap, csv, row.names = FALSE)
saveRDS(list(per_rep = rows, summary = agg, per_rep_diag = diags),
        file.path(OUT, "R18_cvu_combined.rds"))
cat("\nSaved:\n -", csv, "\n -", file.path(OUT, "R18_cvu_combined.rds"), "\n")

# ---- DECISION VIEW: CVu next to the LOCKED CV-w / CF / single-fit FA rows ----
# Same (scenario, rung) cells; samples were seed-identical, so differences are
# attributable to the arm definition, not the draws.
locked_csv <- file.path(RESULTS_DIR, "sim_full_summary.csv")
if (!file.exists(locked_csv)) {
  cat("\n[WARN] locked", locked_csv, "not found -- decision view skipped.\n")
} else {
  cols <- c("scenario","rung","method","n_reps","Psi","bias","emp_sd",
            "mean_se","se_ratio","coverage","mcse_cov")
  lk <- read.csv(locked_csv, stringsAsFactors = FALSE)
  lk <- lk[lk$rung %in% rung_lv &
           lk$method %in% c("Fully-Aware-CV", "Fully-Aware-CF", "Fully-Aware"), cols]
  lk$src <- "locked"
  nw <- agg[, cols]; nw$src <- "R18_new"
  cmp <- rbind(nw, lk)
  ord_m <- c("Fully-Aware-CVu", "Fully-Aware-CV", "Fully-Aware-CF", "Fully-Aware")
  cmp <- cmp[order(cmp$scenario, match(cmp$rung, rung_lv), match(cmp$method, ord_m)), ]
  numc <- sapply(cmp, is.numeric); cp <- cmp; cp[numc] <- round(cp[numc], 4)
  cat("\n==== DECISION VIEW (long): CVu vs locked CV-w / CF / FA, per scenario x rung ====\n")
  print(cp, row.names = FALSE)

  # wide coverage view: the headline decision numbers
  get_cov <- function(s, r, m) {
    v <- cmp$coverage[cmp$scenario == s & cmp$rung == r & cmp$method == m]
    if (length(v)) v[1] else NA_real_
  }
  cellz <- unique(cmp[, c("scenario", "rung")])
  wide <- do.call(rbind, lapply(seq_len(nrow(cellz)), function(i) {
    s <- cellz$scenario[i]; r <- cellz$rung[i]
    data.frame(scenario = s, rung = r,
               cov_CVu = get_cov(s, r, "Fully-Aware-CVu"),
               cov_CVw = get_cov(s, r, "Fully-Aware-CV"),
               cov_CF  = get_cov(s, r, "Fully-Aware-CF"),
               cov_FA  = get_cov(s, r, "Fully-Aware"))
  }))
  wide$CVu_minus_CVw <- wide$cov_CVu - wide$cov_CVw
  wide$CVu_minus_CF  <- wide$cov_CVu - wide$cov_CF
  wide <- wide[order(wide$scenario, match(wide$rung, rung_lv)), ]
  wp <- wide; wn <- sapply(wp, is.numeric); wp[wn] <- round(wp[wn], 3)
  cat("\n==== DECISION VIEW (wide): coverage, CVu vs CV-w vs CF at L2-L3 ====\n")
  print(wp, row.names = FALSE)
  cat("\nINTERPRETATION (see NOTES.md):\n",
      " - EXPECTED: at L3 CVu tracks CV-w (under-covers, ~0.85-0.88) and sits BELOW\n",
      "   CF (~0.94-0.95) -> internal CV != cross-fitting holds for the UNWEIGHTED\n",
      "   object the application reports. At L2 all three should be ~nominal.\n",
      " - SURPRISE: if CVu covers ~nominal at L3 (>= CF), the weighted nuisance fit,\n",
      "   not the missing sample-split, drove the locked CV under-coverage -- STOP\n",
      "   and report before using A4 in the paper.\n", sep = "")
}
