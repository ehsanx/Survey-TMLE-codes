# =====================================================================
# aggregate.R  —  combine R14_dr_factorial per-task RDS into the DR-factorial
# summary (spec item A6, Writing/comments/phase4-arc-sim-specs.md).
# Usage:  Rscript codes/arc_runs/R14_dr_factorial/aggregate.R
#   reads R14_OUT (default DATA_ROOT/arc_runs/R14_dr_factorial)
#   writes results/arc/R14_dr_factorial_summary.csv
#        + R14_dr_factorial_combined.rds under R14_OUT
#
# Summary-stat formulas match codes/aggregate_sim.R EXACTLY:
#   crit = qt(.975, max(1, df)); coverage = mean(|b - Psi| <= crit*se);
#   bias = mean(b) - Psi; emp_sd = sd(b); se_ratio = mean(se)/emp_sd;
#   mcse_cov = sqrt(cov*(1-cov)/n_reps)
# Diverged reps (b/se = NA, the R03 guard) are dropped per group and counted
# in n_diverged. New knob columns (q_spec, g_spec, sampling) lead, per the
# project convention.
#
# Decision views printed at the end:
#   VIEW 1  the 2x2 DR table at L1_param x noninfo: both-correct AND each
#           one-correct cell should be ~unbiased (double robustness),
#           both-wrong biased.
#   VIEW 2  the THREAT view: g-wrong / both-wrong (g_spec = W) cells -- the
#           CF-u minus CF-w bias gap under informative vs non-informative
#           sampling (the Gemini mechanism). PRIMARY: mean of the per-rep
#           PAIRED differences b(CF-u) - b(CF-w) over reps where BOTH CF
#           arms converged (paired MCSE = sd(diff)/sqrt(n_pairs)); the
#           unpaired cell-mean gap is kept as a secondary line. FLAG if the
#           paired gap is materially nonzero under info but ~0 under noninfo.
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R14_OUT", file.path(DATA_ROOT, "arc_runs", "R14_dr_factorial"))
# ^r14_ excludes both the SMOKE_-prefixed smoke outputs and manifest_*;
# list.files() is non-recursive so the manifest/ subdir is never touched.
files <- list.files(OUT, pattern = "^r14_.*\\.rds$", full.names = TRUE)
files <- files[!startsWith(basename(files), "SMOKE_")]
if (!length(files)) stop("no full r14_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

objs <- lapply(files, readRDS)
rows <- do.call(rbind, lapply(objs, function(o)
  cbind(q_spec = o$q_spec, g_spec = o$g_spec, sampling = o$sampling,
        rung = o$rung, Psi = o$Psi,
        o$results[, setdiff(names(o$results), c("q_spec","g_spec","sampling","rung"))])))

z_or_t <- function(df) qt(0.975, pmax(1, df))      # matches aggregate_sim.R
agg <- do.call(rbind, by(rows,
  list(rows$q_spec, rows$g_spec, rows$sampling, rows$rung, rows$method),
  function(d) {
    if (is.null(d)) return(NULL)
    Psi <- d$Psi[1]
    ok  <- is.finite(d$b) & is.finite(d$se)        # drop non-convergent (diverged) reps
    nd  <- sum(!ok); nrok <- sum(ok)
    bo  <- d$b[ok]; so <- d$se[ok]; crit <- z_or_t(d$df[ok])
    cov <- mean(abs(bo - Psi) <= crit * so)
    data.frame(q_spec = d$q_spec[1], g_spec = d$g_spec[1], sampling = d$sampling[1],
               rung = d$rung[1], method = d$method[1],
               n_reps = nrok, n_diverged = nd, Psi = Psi,
               bias = mean(bo) - Psi, emp_sd = sd(bo), mean_se = mean(so),
               se_ratio = mean(so) / sd(bo),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / max(1, nrok)),
               # deff_clust/icc_eif: means of the FA-w-EIF diagnostics over ALL
               # reps INCLUDING diverged ones (the EIF exists even when b/se=NA).
               # In informative cells (weighted fits diverge at the R03 rates)
               # do NOT quote these as clean design-effect estimates --
               # diagnostics only (see NOTES.md caveats).
               deff_clust = mean(d$deff), icc_eif = mean(d$icc_eif),
               w_cv = mean(d$w_cv), df_design = mean(d$df_design),
               stringsAsFactors = FALSE)
  }))
rownames(agg) <- NULL
ord_m <- c("FA-w", "CF-u", "CF-w")
agg <- agg[order(agg$rung, agg$sampling, agg$q_spec, agg$g_spec,
                 match(agg$method, ord_m)), ]

num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 4)
print(ap, row.names = FALSE)

ARC_RES <- file.path(RESULTS_DIR, "arc")
if (!dir.exists(ARC_RES)) dir.create(ARC_RES, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(ARC_RES, "R14_dr_factorial_summary.csv")
write.csv(ap, csv, row.names = FALSE)
rds <- file.path(OUT, "R14_dr_factorial_combined.rds")
saveRDS(list(per_rep = rows, summary = agg), rds)
cat(sprintf("\nSaved:\n - %s\n - %s\n", csv, rds))

# helper: bias MCSE = emp_sd / sqrt(n_reps)
mcse_b <- function(r) r$emp_sd / sqrt(pmax(1, r$n_reps))

# =====================================================================
# VIEW 1 — double-robustness 2x2 at L1_param x noninfo
# Expectation: (C,C), (C,W), (W,C) ~unbiased (|bias| <~ 3*mcse_bias);
# (W,W) clearly biased. (At L1 the wrong spec is sharply Kang-Schafer-
# misspecified; the same table at L2_smooth is printed for the partial-
# recovery contrast.)
# =====================================================================
for (rg in c("L1_param", "L2_smooth")) {
  cat(sprintf("\n==== VIEW 1: DR 2x2 table  [rung=%s, sampling=noninfo] ====\n", rg))
  v1 <- agg[agg$rung == rg & agg$sampling == "noninfo", ]
  if (!nrow(v1)) { cat("  (no rows yet)\n"); next }
  for (mth in ord_m) {
    cat(sprintf("  method %-5s (bias [mcse_bias] / coverage):\n", mth))
    for (q in c("C", "W")) for (g in c("C", "W")) {
      r <- v1[v1$method == mth & v1$q_spec == q & v1$g_spec == g, ]
      if (!nrow(r)) next
      cat(sprintf("    Q=%s g=%s : bias=%+.4f [%.4f]  cov=%.3f  (n=%d, div=%d)\n",
                  q, g, r$bias, mcse_b(r), r$coverage, r$n_reps, r$n_diverged))
    }
    # automated DR readout for this method (only decisive at L1)
    dr  <- v1[v1$method == mth & !(v1$q_spec == "W" & v1$g_spec == "W"), ]
    ww  <- v1[v1$method == mth &   v1$q_spec == "W" & v1$g_spec == "W", ]
    if (nrow(dr) == 3 && nrow(ww) == 1) {
      dr_ok <- all(abs(dr$bias) <= 3 * mcse_b(dr))
      ww_bi <- abs(ww$bias) > 3 * mcse_b(ww)
      cat(sprintf("    -> DR cells ~unbiased: %s | both-wrong biased: %s  ==> %s\n",
                  dr_ok, ww_bi,
                  if (dr_ok && ww_bi) "DR PATTERN CONFIRMED"
                  else if (rg == "L1_param") "DR PATTERN NOT CONFIRMED -- inspect before using"
                  else "see L1 for the decisive read (L2 smooth library partially recovers W)"))
    }
  }
}

# =====================================================================
# VIEW 2 — informative-sampling unweighted-OOF threat (the Gemini mechanism)
# g-wrong or both-wrong (g_spec = W) cells, info vs noninfo sampling.
# PRIMARY (decision) gap: mean of the PER-REP PAIRED differences
#   diff_i = b_i(CF-u) - b_i(CF-w),
# joined on rep within each cell and restricted to reps where BOTH CF arms
# converged; paired MCSE = sd(diff)/sqrt(n_pairs). The arms share the sample
# and the folds, so pairing is exact -- and conditioning BOTH arms on the
# SAME rep subset removes the selection asymmetry of unpaired cell means:
# CF-w diverges on ~17% of informative L1 draws (locked R03 analogue) while
# CF-u converges on ~all, and divergence is concentrated in exactly the
# informative cells the FLAG keys on, so E[b_CFw | CF-w converged] could
# manufacture or mask an "informative-only gap" through selection rather
# than the de-weighting mechanism.
# SECONDARY (reference only): the unpaired cell-mean gap bias(CF-u) -
# bias(CF-w), each arm over its OWN convergent subset, with an MCSE that
# treats the arms as independent.
# FLAG when |paired gap| > 2*paired MCSE under info but <= under noninfo.
# =====================================================================
cat("\n==== VIEW 2: THREAT — CF-u vs CF-w bias gap in g-wrong cells ====\n")
cat("  primary   = paired per-rep diff b(CF-u)-b(CF-w), reps where BOTH CF arms converged\n")
cat("  secondary = unpaired per-arm cell means (each arm's own convergent subset)\n")
threat_flagged <- FALSE
.cf_pairs <- function(d) {                  # d = per-rep rows of ONE cell
  conv <- is.finite(d$b) & is.finite(d$se)  # same convergence rule as `agg`
  u <- d[d$method == "CF-u" & conv, c("rep", "b")]
  w <- d[d$method == "CF-w" & conv, c("rep", "b")]
  ndw <- sum(d$method == "CF-w") - nrow(w)  # CF-w diverged count in this cell
  if (!nrow(u) || !nrow(w))
    return(list(n_pairs = 0L, gap = NA_real_, mcse = NA_real_, nd_w = ndw))
  m   <- merge(u, w, by = "rep")            # join on rep within the cell
  dif <- m$b.x - m$b.y                      # b(CF-u) - b(CF-w), paired
  list(n_pairs = nrow(m), gap = mean(dif),
       mcse = if (nrow(m) >= 2L) sd(dif) / sqrt(nrow(m)) else NA_real_,
       nd_w = ndw)
}
for (rg in c("L1_param", "L2_smooth")) for (q in c("C", "W")) {
  pp <- list()
  ug <- mg <- setNames(rep(NA_real_, 2), c("info", "noninfo"))
  any_rows <- FALSE
  for (sm in c("info", "noninfo")) {
    cell <- rows[rows$rung == rg & rows$q_spec == q & rows$g_spec == "W" &
                 rows$sampling == sm & rows$method %in% c("CF-u", "CF-w"), ]
    if (nrow(cell)) any_rows <- TRUE
    pp[[sm]] <- .cf_pairs(cell)
    # secondary: unpaired cell-mean gap (kept from the original view)
    u <- agg[agg$rung == rg & agg$q_spec == q & agg$g_spec == "W" &
             agg$sampling == sm & agg$method == "CF-u", ]
    w <- agg[agg$rung == rg & agg$q_spec == q & agg$g_spec == "W" &
             agg$sampling == sm & agg$method == "CF-w", ]
    if (nrow(u) == 1 && nrow(w) == 1) {
      ug[sm] <- u$bias - w$bias
      mg[sm] <- sqrt(mcse_b(u)^2 + mcse_b(w)^2)
    }
  }
  if (!any_rows) next                       # cell not run yet
  decisive <- vapply(pp, function(p) p$n_pairs >= 2L && is.finite(p$mcse),
                     logical(1))
  flag <- all(decisive) &&
          abs(pp$info$gap)    >  2 * pp$info$mcse &&
          abs(pp$noninfo$gap) <= 2 * pp$noninfo$mcse
  if (flag) threat_flagged <- TRUE
  fmt_p <- function(p) if (p$n_pairs >= 2L && is.finite(p$mcse))
      sprintf("%+.4f [mcse %.4f, n_pairs=%d, CF-w div=%d]",
              p$gap, p$mcse, p$n_pairs, p$nd_w)
    else sprintf("NA [n_pairs=%d, CF-w div=%d]", p$n_pairs, p$nd_w)
  cat(sprintf("  [%s, Q=%s, g=W] PAIRED gap(CF-u - CF-w):  info %s   noninfo %s   %s\n",
              rg, q, fmt_p(pp$info), fmt_p(pp$noninfo),
              if (flag) ">>> FLAG: unweighted OOF loses DR under informative sampling"
              else if (all(decisive)) "no informative-only gap"
              else "insufficient pairs -- no decision"))
  fmt_u <- function(sm) if (is.finite(ug[sm]))
      sprintf("%+.4f [mcse %.4f]", ug[sm], mg[sm]) else "NA"
  cat(sprintf("      unpaired secondary:                   info %s   noninfo %s\n",
              fmt_u("info"), fmt_u("noninfo")))
}
cat(if (threat_flagged)
      "\n>>> THREAT CONFIRMED in >=1 cell: report as the limitation (unweighted misspecified\n    library converges to the SAMPLE projection under informative sampling; weights restore it).\n"
    else
      "\n>>> No informative-only CF-u/CF-w gap detected: the unweighted-OOF default survives\n    this DR factorial; report as a robustness result.\n")
