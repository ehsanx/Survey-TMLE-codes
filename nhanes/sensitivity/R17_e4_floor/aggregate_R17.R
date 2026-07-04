# =====================================================================
# Nhanes/arc_runs/R17_e4_floor/aggregate_R17.R
#   Combine R17 per-cell nh17_*.rds  ->  results/arc/R17_floor_sensitivity.csv
#                                        results/arc/R17_floor_share.csv
#                                        <R17_OUT>/R17_combined.rds
# Mirrors aggregate_R05.R conventions; writes ONLY to the NON-destructive
# results/arc/ location (never the locked nhanes_output/results).
# Usage:  Rscript Nhanes/arc_runs/R17_e4_floor/aggregate_R17.R
#
# NOTE: this is a REAL-DATA run (no simulation truth), so there is no
# coverage/bias column; the locked NHANES table convention applies
# (b = mean over the B splits, b_split_sd, se = mean design SE, t-based
# lcl/ucl with crit = qt(.975, max(1, df))), plus the A9b floor/share knob
# columns leading.
# =====================================================================

REPO <- Sys.getenv("REPO_ROOT", ".")
OUT  <- Sys.getenv("R17_OUT", file.path(REPO, "nhanes", "nhanes_output", "arc_runs", "R17"))
RES  <- Sys.getenv("R17_RESULTS", file.path(REPO, "results", "arc"))
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# per-cell files only: the ^nh17_ anchor skips SMOKE_* outputs; manifest/ is a
# subdirectory and list.files() is non-recursive, so it is skipped too.
files <- list.files(OUT, pattern = "^nh17_.*\\.rds$", full.names = TRUE)
files <- files[!startsWith(basename(files), "SMOKE_")]   # belt-and-suspenders
if (!length(files)) stop("no nh17_*.rds found in ", OUT, " -- run the array first")
cat("aggregating", length(files), "R17 cell files from", OUT, "\n")
objs <- lapply(files, readRDS)

## (1) full summary: every cell, knob columns (floor) leading ------------------
summary_tbl <- do.call(rbind, lapply(objs, `[[`, "summary"))
summary_tbl <- summary_tbl[, c("floor", "example", "label", "method", "B",
                               "b", "b_split_sd", "se", "df", "lcl", "ucl",
                               "share_clip_w", "share_clip_unw",
                               "mass05_w", "expmass05_w", "g_raw_min", "g_raw_max")]
summary_tbl <- summary_tbl[order(summary_tbl$example, -summary_tbl$floor), ]
rownames(summary_tbl) <- NULL

rnd <- function(d, k = 4) { n <- sapply(d, is.numeric); d[n] <- round(d[n], k); d }

## (2) R17_floor_sensitivity.csv: the E4 rows across floors --------------------
# floor, b (mean over splits), b_split_sd, se, df, t-based lcl/ucl, mean shares.
sens_tbl <- summary_tbl[summary_tbl$example == "E4",
                        c("floor", "example", "label", "method", "B",
                          "b", "b_split_sd", "se", "df", "lcl", "ucl",
                          "share_clip_w", "mass05_w", "expmass05_w")]
sens_tbl <- sens_tbl[order(-sens_tbl$floor), ]
write.csv(rnd(sens_tbl), file.path(RES, "R17_floor_sensitivity.csv"), row.names = FALSE)

## (3) R17_floor_share.csv: E2/E3/E4 at the locked 0.05 floor (the N6 table) ---
share_tbl <- do.call(rbind, lapply(objs, function(o) {
  if (!isTRUE(all.equal(o$floor, 0.05))) return(NULL)
  s <- o$summary; dg <- o$diagnostics
  data.frame(example = o$example, label = o$label, floor = o$floor,
             n_domain = dg$n_domain, A_prev = dg$A_prev,
             share_clip_w = s$share_clip_w, share_clip_unw = s$share_clip_unw,
             mass05_w = s$mass05_w, expmass05_w = s$expmass05_w,
             g_raw_min = s$g_raw_min, g_raw_max = s$g_raw_max,
             b = s$b, se = s$se, lcl = s$lcl, ucl = s$ucl)
}))
# guard: a partial run with only the f025/f010 cells yields share_tbl == NULL
# (every lapply element returns NULL) -- skip the table instead of crashing.
share_csv <- file.path(RES, "R17_floor_share.csv")
if (is.null(share_tbl)) {
  cat("no 0.05-floor cells aggregated yet -- share table skipped",
      "(only f025/f010 cells present; re-run after tasks 1-3 finish).\n")
} else {
  share_tbl <- share_tbl[order(share_tbl$example), ]
  write.csv(rnd(share_tbl), share_csv, row.names = FALSE)
}

## (4) combined object (per-split rows + split-1 g_raw vectors for figures) ----
saveRDS(list(summary = summary_tbl, floor_sensitivity = sens_tbl, floor_share = share_tbl,
             per_split = do.call(rbind, lapply(objs, function(o)
               cbind(example = o$example, o$per_split))),
             cells = objs),
        file.path(OUT, "R17_combined.rds"))

## ---- DECISION VIEWS ----------------------------------------------------------
cat("\n========== R17 DECISION VIEW 1: E4 point-estimate stability across floors ==========\n")
print(rnd(sens_tbl), row.names = FALSE)
if (nrow(sens_tbl) >= 2) {
  bb <- sens_tbl$b
  max_db <- max(abs(outer(bb, bb, "-")))
  cat(sprintf("\n  max |b(f_i) - b(f_j)| across floors = %.4f\n", max_db))
  cat(sprintf("  vs mean design se                   = %.4f  (ratio %.2f)\n",
              mean(sens_tbl$se), max_db / mean(sens_tbl$se)))
  cat(sprintf("  vs max b_split_sd                   = %.4f  (ratio %.2f)\n",
              max(sens_tbl$b_split_sd), max_db / max(sens_tbl$b_split_sd)))
  cat("  READ: the E4 estimate is FLOOR-STABLE if max|db| is small relative to the\n",
      "  design se (ratio << 1); the same-split raw OOF ghat is identical across the\n",
      "  floor cells, so max|db| is the pure floor effect (no split noise).\n", sep = "")
}

cat("\n========== R17 DECISION VIEW 2: share-at-floor across examples (N6 table, floor=0.05) ==========\n")
if (is.null(share_tbl)) {
  cat("  (skipped -- no 0.05-floor cells aggregated yet; only f025/f010 present.)\n")
} else {
  print(rnd(share_tbl), row.names = FALSE)
  cat("\n  READ: share_clip_w / expmass05_w quantify WHERE the truncated (overlap-\n",
      "  restricted) estimand reframe is needed: expect E4 (rare exposure) >> E2/E3.\n", sep = "")
}

cat("\nSaved:\n -", file.path(RES, "R17_floor_sensitivity.csv"), "\n")
if (!is.null(share_tbl)) cat(" -", share_csv, "\n")
cat(" -", file.path(OUT, "R17_combined.rds"), "\n")
