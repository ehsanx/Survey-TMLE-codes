# =====================================================================
# estimators_r17.R  —  R17 helper: primary-arm cross-fit at a TUNABLE floor
#                      + RAW-OOF-propensity floor/overlap records (spec A9b)
#
# PURPOSE (spec item A9b, Writing/comments/phase4-arc-sim-specs.md):
# For E4 (6% exposure; single-fit min ghat -> ~0.003 post-N1) the locked 0.05
# OOF propensity floor BINDS. This helper re-runs ONLY the primary
# estimator (Fully-Aware-CF: PSU-within-stratum cross-fit, UNWEIGHTED per-fold
# nuisances, weighted POOLED tmle targeting on the domain, domain design SE)
# with the OOF propensity floored at a CELL-SPECIFIC g_floor and the EIF
# reconstruction gbound = the SAME g_floor, and records, per split, the
# share-of-units-at-the-floor quantities computed from the RAW (PRE-clipping)
# OOF ghat over the domain rows with the domain weights w:
#
#   share_clip_w   = sum(w * (g_raw < g_floor)) / sum(w)     weighted share clipped
#   share_clip_unw = mean(g_raw < g_floor)                   unweighted share clipped
#   mass05_w       = sum(w * (g_raw < 0.05)) / sum(w)        FIXED 0.05 reference mass
#                                                            (same regardless of cell floor)
#   expmass05_w    = sum(w[A==1 & g_raw<0.05]) / sum(w[A==1])  weighted share of the
#                                                            EXPOSED below 0.05
#   g_raw min/max                                            raw OOF overlap range
#
# These feed (i) the E4 floor-sensitivity table (floor in {0.05, 0.025, 0.01})
# and (ii) the N6 share-at-floor table across E2/E3/E4 at the locked 0.05 floor
# (the estimand-reframe: truncated / overlap-restricted target where the floor
# binds).
#
# REUSE (canonical building blocks from codes/estimators.R, sourced read-only):
#   make_cf_folds(strata, cluster, V)  — PSU-within-stratum cross-fit folds
#   .sl(Y, X, weights, learners)       — per-fold SuperLearner predictor closure
#   .eif_from_tmle(fit, Y, A, w, gbound) — EIF reconstruction, TUNABLE floor
#   .se_des(eif, strata, cluster, weight, clustered, nest, inpop)
#                                      — design SE (Eq 8), FULL-design-subset
#                                        domain path (inpop), nest=TRUE for NHANES
# The domain handling (sub <- which(inpop), full-design subset for the variance)
# is byte-identical to run_estimators() / run_estimators_iso() so the design df
# and domain variance match the locked pipeline exactly.
#
# FLOOR MECHANICS (why the cell floor really binds in targeting):
#   * the OOF propensity is clipped to [g_floor, 1-g_floor] BEFORE tmle(), exactly
#     as the locked CF arm clips at [0.05, 0.95];
#   * tmle 2.1.1's internal default gbound (gbound = NULL -> 5/sqrt(n)/log(n);
#     ~0.0051 for the E3/E4 domains, ~0.0028 for E2; ATE bounds = [lb, 1]) is
#     BELOW every cell floor here (0.0051 < 0.01 < 0.025 < 0.05), so .bound is a
#     no-op on the pre-clipped g and the supplied g1W enters targeting unchanged.
#     This is a SOURCE-INSPECTION argument, not a runtime-verifiable one: tmle
#     returns the user-supplied g1W unmodified in cf$g$g1W (its internal bounding
#     applies to an unexposed g1W.total), so the SMOKE-gate g_used_min >= floor
#     line is only a cheap sanity invariant on OUR clipping, never a passthrough
#     verification;
#   * the EIF clever covariate is re-bounded at gbound = g_floor via the exposed
#     .eif_from_tmle — a no-op on the already-clipped g, kept for explicitness.
#   The OOF Q is clipped at [1e-3, 1-1e-3] EXACTLY as the locked CF arm (the A9b
#   floor knob is the PROPENSITY floor only; Q-truncation is not varied).
#
# RNG NOTE (exact floor attribution): the caller set.seed()s per split BEFORE
# calling this function, and make_cf_folds is the FIRST RNG consumer here, so for
# a given split seed the folds, the per-fold SL fits, and hence the RAW OOF ghat
# are IDENTICAL across the three E4 floor cells; the floor enters only via the
# post-hoc clipping + targeting. Differences in b across floors are therefore
# attributable to the floor alone (no split noise between floor cells).
# =====================================================================

# Sourced AFTER codes/estimators.R (which defines the building blocks and loads
# SuperLearner/survey/tmle). We do not re-library here to avoid masking.

# ---- primary-arm runner with a tunable propensity/EIF floor ----------------
# Args mirror run_estimators(); `g_floor` replaces the fixed g_oof_bound = 0.05
# AND sets the EIF gbound. inpop is REQUIRED (NHANES domain runs only).
# Returns: list(row  = 1-row data.frame of per-split records,
#               eif  = the targeted domain EIF (for deff_clust),
#               g_raw = the RAW pre-clipping OOF ghat over the domain rows).
run_cf_floor <- function(obs, learners, V_cf = 5L, g_floor = 0.05,
                         W_cols, nest = TRUE, inpop) {
  stopifnot(!is.null(inpop), g_floor > 0, g_floor < 0.5)
  # ----- domain slice: identical contract to run_estimators() -----
  sub <- which(inpop)
  d   <- obs[sub, , drop = FALSE]
  n   <- nrow(d); w <- d$weight; A <- d$A
  W   <- d[, W_cols, drop = FALSE]
  XA  <- d[, c("A", W_cols), drop = FALSE]
  Xa1 <- data.frame(A = 1, W); Xa0 <- data.frame(A = 0, W)
  # design-SE closure: FULL design subset()-ed to the domain (locked convention)
  se_des <- function(eif, clustered)
    .se_des(eif, obs$strata, obs$cluster, obs$weight, clustered,
            nest = nest, inpop = inpop)

  ## ---- primary Fully-Aware-CF, floor = g_floor --------------------------
  ## Structure byte-identical to the locked CF arm in run_estimators(); the ONLY
  ## changes are g_oof_bound -> g_floor and the EIF gbound -> g_floor.
  fold <- make_cf_folds(d$strata, d$cluster, V_cf)         # FIRST RNG consumer
  Q0o <- Q1o <- g1o <- numeric(n)
  for (v in sort(unique(fold))) {
    tr <- which(fold != v); ho <- which(fold == v)
    qf <- .sl(d$Y[tr], XA[tr, ], weights = NULL, learners) # UNWEIGHTED (ignorability)
    gf <- .sl(d$A[tr], W[tr, ],  weights = NULL, learners)
    Q1o[ho] <- qf(Xa1[ho, ]); Q0o[ho] <- qf(Xa0[ho, ]); g1o[ho] <- gf(W[ho, ])
  }
  g_raw <- g1o                                             # RAW pre-clipping OOF ghat
  g1c   <- pmin(pmax(g_raw, g_floor), 1 - g_floor)         # cell-floor clipping
  Qoo   <- pmin(pmax(cbind(Q0o, Q1o), 1e-3), 1 - 1e-3)     # Q clip EXACTLY as locked
  cf <- tmle(Y = d$Y, A = d$A, W = W, family = "binomial", obsWeights = w,
             Q = Qoo, g1W = g1c)                           # weighted POOLED targeting
  e_cf <- .eif_from_tmle(cf, d$Y, d$A, w, gbound = g_floor)  # EIF gbound = cell floor
  s_cf <- se_des(e_cf$eif, TRUE)

  ## ---- A9b floor/overlap records from the RAW OOF ghat --------------------
  ## All weighted shares use the DOMAIN weights w over the DOMAIN rows.
  ## Lower tail only (rare exposure -> the floor binds from below); the upper
  ## tail is monitored via g_raw_max (and share_clip_hi_w, a free diagnostic).
  sw <- sum(w)
  row <- data.frame(
    method          = "Fully-Aware-CF",
    floor           = g_floor,
    b               = e_cf$psi,
    se              = s_cf$se,
    df              = s_cf$df,
    # --- share-at-floor (cell floor; RAW ghat, pre-clipping) ---
    share_clip_w    = sum(w * (g_raw < g_floor)) / sw,
    share_clip_unw  = mean(g_raw < g_floor),
    # --- FIXED 0.05 reference mass (comparable across floor cells) ---
    mass05_w        = sum(w * (g_raw < 0.05)) / sw,
    # --- weighted share of the EXPOSED below 0.05 ---
    expmass05_w     = sum(w[A == 1 & g_raw < 0.05]) / sum(w[A == 1]),
    # --- raw + as-used overlap ranges ---
    g_raw_min       = min(g_raw),  g_raw_max = max(g_raw),
    g_used_min      = min(cf$g$g1W), g_used_max = max(cf$g$g1W),  # echoes the supplied g1W (clip sanity only)
    share_clip_hi_w = sum(w * (g_raw > 1 - g_floor)) / sw,        # upper-tail monitor
    # --- targeting diagnostics ---
    eps_cf          = max(abs(cf$epsilon)),
    cf_V_eff        = attr(fold, "V_eff"))
  list(row = row, eif = e_cf$eif, g_raw = g_raw)
}
