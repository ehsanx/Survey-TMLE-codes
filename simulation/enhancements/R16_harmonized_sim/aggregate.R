# =====================================================================
# aggregate.R  —  combine R16_harmonized_sim per-task RDS into a summary CSV
# Usage:  Rscript aggregate.R     (after the full --array=1-80 job completes)
#
# Mirrors codes/aggregate_sim.R conventions EXACTLY (crit = qt(.975, max(1,df));
# coverage = mean(|b-Psi| <= crit*se); bias = mean(b)-Psi; emp_sd = sd(b);
# se_ratio = mean(se)/emp_sd; mcse_cov = sqrt(cov*(1-cov)/n_reps)) so the result
# drops straight into the figures. Writes results/arc/R16_harmonized_sim_summary.csv
# + R16_harmonized_sim_combined.rds under R16_OUT -- it does NOT clobber the
# locked results/sim_full_summary.csv.
#
# DECISION VIEW: reads the LOCKED results/sim_full_summary.csv and prints
# harmonized-vs-locked coverage and se_ratio side by side per (scenario, rung,
# base arm). Headline question: does the L4 single-fit undercoverage SURVIVE
# harmonization (expected YES per R03_isolation_2x2)?
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R16_OUT", file.path(DATA_ROOT, "arc_runs", "R16_harmonized_sim"))
# full-run task files only: skip SMOKE_* and the manifest/ subdir (non-recursive
# listing already excludes manifest/; the ^r16_ anchor already excludes SMOKE_,
# the explicit filter is belt-and-braces).
files <- list.files(OUT, pattern = "^r16_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("^SMOKE_", basename(files))]
if (!length(files)) stop("no full r16_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)

# ---- mixed-knob guard --------------------------------------------------------
# Per-task FILENAMES do not encode g_floor/q_lo, and the per-chunk checkpoint
# skips already-written files -- so re-submitting with changed R16_GFLOOR/R16_QLO
# into the SAME R16_OUT would silently pool mixed-floor chunks. Refuse to
# aggregate unless the knobs are constant across ALL loaded task files.
knob_tbl <- data.frame(file    = basename(files),
                       g_floor = vapply(objs, function(o) o$g_floor, numeric(1)),
                       q_lo    = vapply(objs, function(o) o$q_lo,    numeric(1)),
                       stringsAsFactors = FALSE)
uk <- unique(knob_tbl[, c("g_floor", "q_lo")])
if (nrow(uk) > 1L) {
  cat("\nERROR: mixed truncation knobs across the loaded task files:\n")
  for (k in seq_len(nrow(uk))) {
    sel <- knob_tbl$g_floor == uk$g_floor[k] & knob_tbl$q_lo == uk$q_lo[k]
    ex  <- knob_tbl$file[sel]
    cat(sprintf("  g_floor=%g q_lo=%g : %d file(s) -- %s%s\n",
                uk$g_floor[k], uk$q_lo[k], sum(sel),
                paste(head(ex, 5), collapse = ", "),
                if (length(ex) > 5) sprintf(" (+%d more)", length(ex) - 5) else ""))
  }
  stop("g_floor/q_lo are not constant across the loaded chunks. The checkpoint ",
       "cannot distinguish knob values, so changed R16_GFLOOR/R16_QLO chunks need ",
       "a FRESH R16_OUT (or delete the old chunks) before aggregating.")
}

rows <- do.call(rbind, lapply(objs, function(o)
  cbind(g_floor = o$g_floor, q_lo = o$q_lo,
        scenario = o$scenario, rung = o$rung, Psi = o$Psi, o$results)))

z_or_t <- function(df) qt(0.975, pmax(1, df))
agg <- do.call(rbind, by(rows, list(rows$scenario, rows$rung, rows$method), function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  Psi <- d$Psi[1]
  ok  <- is.finite(d$b) & is.finite(d$se)     # drop non-convergent (diverged) reps
  nd  <- sum(!ok); nrok <- sum(ok)
  bo  <- d$b[ok]; so <- d$se[ok]; crit <- z_or_t(d$df[ok])
  cov <- mean(abs(bo - Psi) <= crit * so)
  data.frame(g_floor = d$g_floor[1], q_lo = d$q_lo[1],       # knob columns LEADING
             scenario = d$scenario[1], rung = d$rung[1], method = d$method[1],
             n_reps = nrok, n_diverged = nd, Psi = Psi,
             bias = mean(bo) - Psi, emp_sd = sd(bo), mean_se = mean(so),
             se_ratio = mean(so) / sd(bo),
             coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / max(1, nrok)),
             # deff/icc are CELL-level diagnostics (per-rep values from the FA-h
             # EIF, duplicated across every method row of a rep): average over
             # the FULL rep set -- the same set for every method row of a cell --
             # not the per-method ok-subset (matches aggregate_sim.R / R03).
             deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif))
}))
rownames(agg) <- NULL
ord_m <- c("Fully-Aware-h", "Fully-Aware-CF-h", "Fully-Aware-CV-h",
           "Partially-Aware-h", "Non-Aware-h")
agg <- agg[order(agg$scenario, agg$rung, match(agg$method, ord_m)), ]  # L1_.. ladder order

num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
cat("\n==== R16 harmonized-truncation summary (all five arms, one rule) ====\n")
print(ap, row.names = FALSE)

# ---- errored-rep accounting (WARN, not stop) ---------------------------------
# Hard-errored reps (run.R one_rep tryCatch -> NULL) never reach the per-task
# rows, so they are visible only as a shortfall: n_reps (finite) + n_diverged
# (guard-flagged NAs) should equal the expected rep count per cell.
exp_reps <- as.integer(Sys.getenv("R16_EXPECTED_REPS", "1000"))
short <- agg[agg$n_reps + agg$n_diverged < exp_reps, , drop = FALSE]
if (nrow(short)) {
  cat(sprintf("\nWARN: %d summary row(s) short of the expected %d reps/cell (errored or not-yet-run reps; set R16_EXPECTED_REPS if intentional):\n",
              nrow(short), exp_reps))
  print(short[, c("scenario", "rung", "method", "n_reps", "n_diverged")],
        row.names = FALSE)
}

dir.create(file.path(RESULTS_DIR, "arc"), recursive = TRUE, showWarnings = FALSE)
csv <- file.path(RESULTS_DIR, "arc", "R16_harmonized_sim_summary.csv")
write.csv(ap, csv, row.names = FALSE)
combined <- file.path(OUT, "R16_harmonized_sim_combined.rds")
saveRDS(list(per_rep = rows, summary = agg), combined)
cat("\nSaved:\n -", csv, "\n -", combined, "\n")

# ---- DECISION VIEW: harmonized (h) vs LOCKED headline, side by side ----------
locked_csv <- file.path(RESULTS_DIR, "sim_full_summary.csv")
cat("\n==== DECISION VIEW: harmonized vs locked, per (scenario, rung, base arm) ====\n")
if (!file.exists(locked_csv)) {
  cat("locked summary", locked_csv, "not found -- skipping the side-by-side view.\n",
      "(Run this aggregate where results/sim_full_summary.csv is present.)\n")
} else {
  lk <- read.csv(locked_csv, stringsAsFactors = FALSE)
  base_map <- c("Fully-Aware-h"     = "Fully-Aware",
                "Partially-Aware-h" = "Partially-Aware",
                "Non-Aware-h"       = "Non-Aware",
                "Fully-Aware-CV-h"  = "Fully-Aware-CV",
                "Fully-Aware-CF-h"  = "Fully-Aware-CF")
  cmp <- agg[, c("scenario", "rung", "method", "coverage", "se_ratio")]
  cmp$base <- unname(base_map[cmp$method])
  m <- merge(cmp, lk[, c("scenario", "rung", "method", "coverage", "se_ratio")],
             by.x = c("scenario", "rung", "base"),
             by.y = c("scenario", "rung", "method"),
             suffixes = c("_h", "_locked"), all.x = TRUE)
  m$d_cov <- m$coverage_h - m$coverage_locked
  m <- m[order(m$scenario, m$rung, match(m$method, ord_m)),
         c("scenario", "rung", "method", "coverage_h", "coverage_locked", "d_cov",
           "se_ratio_h", "se_ratio_locked")]
  mn <- sapply(m, is.numeric); mp <- m; mp[mn] <- round(mp[mn], 4)
  print(mp, row.names = FALSE)

  # ---- the headline question -------------------------------------------------
  cat("\n==== HEADLINE: does the L4 single-fit undercoverage SURVIVE harmonization? ====\n")
  fa <- m[m$rung == "L4_aggressive" & m$method == "Fully-Aware-h", ]
  cf <- m[m$rung == "L4_aggressive" & m$method == "Fully-Aware-CF-h", ]
  for (i in seq_len(nrow(fa))) {
    j <- match(fa$scenario[i], cf$scenario)
    cat(sprintf("  %-8s  FA-h L4 coverage = %.3f (locked FA %.3f)  |  CF-h = %.3f (locked CF %.3f)\n",
                fa$scenario[i], fa$coverage_h[i], fa$coverage_locked[i],
                if (is.na(j)) NA_real_ else cf$coverage_h[j],
                if (is.na(j)) NA_real_ else cf$coverage_locked[j]))
  }
  # Verdict per the NOTES decision rule: CONFIRMS needs FA-h < 0.90 in BOTH
  # scenarios AND the CF-h anchor ~nominal (>= 0.93) in those same L4 cells.
  scen_need <- c("standard", "R1")
  both_present <- all(scen_need %in% fa$scenario) && all(scen_need %in% cf$scenario)
  if (!both_present) {
    cat("  (no L4_aggressive Fully-Aware-h / Fully-Aware-CF-h rows for BOTH scenarios\n",
        "   aggregated yet -- run incomplete; no verdict. Re-run after tasks 61-80\n",
        "   finish.)\n", sep = "")
  } else {
    cf_ok <- all(is.finite(cf$coverage_h)) && all(cf$coverage_h >= 0.93)
    cat(sprintf("  CF-h anchor check (coverage_h >= 0.93 in both L4 cells): %s  [%s]\n",
                if (cf_ok) "PASS" else "FAIL",
                paste(sprintf("%s %.3f", cf$scenario, cf$coverage_h), collapse = ", ")))
    fa_collapsed <- all(is.finite(fa$coverage_h)) && all(fa$coverage_h < 0.90)
    if (fa_collapsed && cf_ok) {
      cat("  >>> CONFIRMS (expected per R03): the L4 single-fit undercoverage SURVIVES\n",
          "      the harmonized truncation while CF-h holds ~nominal -> the Figure 1 /\n",
          "      Table 1 contrast is NOT a floor artifact; cross-fitting remains the\n",
          "      active ingredient. The harmonized table can be reported as the\n",
          "      robustness companion (or adopted as primary per A9a's 'Targets'\n",
          "      note).\n", sep = "")
    } else if (!fa_collapsed) {
      cat("  >>> STOP-AND-REPORT: harmonized truncation rescued single-fit coverage at\n",
          "      L4 (>= 0.90 in some scenario) -> the locked headline contrast is\n",
          "      (partly) a truncation artifact. Do NOT adopt either table as primary\n",
          "      before review; reconcile with R03's isolation result first.\n", sep = "")
    } else {
      cat("  >>> STOP-AND-REPORT: FA-h collapsed as expected BUT the CF-h anchor failed\n",
          "      (coverage < 0.93 in some L4 cell; locked CF: 0.985/0.992) -> suspect a\n",
          "      bug in the harmonized CF path (it is locked-CF-identical up to the EIF\n",
          "      gbound, a no-op at the default floor). Do NOT report either table\n",
          "      before reconciling.\n", sep = "")
    }
  }
}
