# =====================================================================
# dr_helpers.R  —  R14_dr_factorial helper: per-nuisance misspecification
# arms for the double-robustness factorial (spec item A6,
# Writing/comments/phase4-arc-sim-specs.md).
#
# Lives in THIS run folder; codes/*.R is sourced READ-ONLY and never edited.
# Reuses the building blocks EXPOSED by sourcing codes/estimators.R:
#   .sl(Y, X, weights, learners)              SuperLearner -> predictor closure
#   .eif_from_tmle(fit, Y, A, w, gbound)      reconstruct per-unit EIF
#   .se_des(eif, strata, cluster, weight, clustered)  Eq-8 design SE + df
#   make_cf_folds(strata, cluster, V)         whole-PSU-within-stratum folds
# and codes/dgp.R:
#   apply_L(C1, C2, C3, C4, "complex")        the Kang-Schafer transforms
# plus the `tmle` package for the (weighted, pooled) targeting step, exactly
# as the engine's CV/CF arms use it (estimators.R lines 177-178, 195-196).
#
# THE COVARIATE DEVICE (locked design):
#   The sample is drawn with model_type = "simple", so obs$L1..L4 carry the
#   RAW C1..C4 (apply_L is the identity for "simple"; dgp.R lines 83-85).
#   * CORRECT frame ("C")  = the raw columns  -> the learner is correctly
#     specified at L1 (A and Y are generated linearly in C1..C4 on the logit
#     scale).
#   * WRONG frame  ("W")   = apply_L(L1..L4, "complex") = the Kang-Schafer
#     transforms, the canonical misspecification device of the headline run.
#   Column names stay L1..L4 in both frames so every learner call is
#   structurally identical; only what the columns MEAN changes. The Q-learner
#   sees q_spec's frame (+A) and the g-learner sees g_spec's frame,
#   INDEPENDENTLY -> a clean 2x2 Q-correct/wrong x g-correct/wrong factorial.
#
# THE THREE ARMS (per rep; run_dr_arms):
#   FA-w  weighted SINGLE-FIT nuisances (.sl(weights = w)) on the spec'd
#         frames -> g floored at g_floor=0.05, Q clamped 1e-3 -> weighted
#         pooled tmle(Q=, g1W=) targeting -> clustered design SE; EIF
#         reconstructed with gbound = 0.05.  (R03's SF-W pattern, harmonized
#         floor, with split design frames.)
#   CF-u  UNWEIGHTED per-fold out-of-fold fits (the paper-default CF choice)
#         over make_cf_folds; OOF g floor 0.05; same weighted pooled
#         targeting; clustered SE; EIF gbound 0.05.
#   CF-w  WEIGHTED per-fold OOF fits on the SAME folds (paired design, as in
#         R07/cf_arms.R) -- the only difference from CF-u is whether the OOF
#         fits receive the survey weights, so the CF-u minus CF-w bias gap is
#         attributable solely to de-weighting (the Gemini mechanism: an
#         unweighted MISSPECIFIED library converges to the SAMPLE L2
#         projection, which under informative sampling differs from the
#         population projection).
#
# DIVERGENCE GUARD (from R03/estimators_isolation.R): weighted nuisance fits
# can separate and the pooled targeting then diverges (eps -> ~1e13); such a
# rep would wreck a 1000-rep mean, so reps with max|epsilon| > 20 are marked
# NON-CONVERGENT (b/se = NA). The aggregator counts them (n_diverged).
# =====================================================================

.R14_Q_CLAMP <- 1e-3   # Q truncation handed to tmle() (engine convention)
.R14_G_FLOOR <- 0.05   # harmonized propensity floor: nuisance preds AND EIF gbound

# ---- the two design frames ("C" = correct/raw, "W" = wrong/Kang-Schafer) ----
make_dr_frames <- function(obs) {
  list(
    C = obs[, c("L1", "L2", "L3", "L4"), drop = FALSE],
    W = apply_L(obs$L1, obs$L2, obs$L3, obs$L4, "complex")
  )
}

# ---- shared targeting step: pre-fit (Q, g1) -> weighted pooled tmle -> SE ----
.dr_target <- function(d, W_targ, Q, g1, method_label, g_floor, V_eff) {
  w <- d$weight
  fit <- tmle(Y = d$Y, A = d$A, W = W_targ, family = "binomial", obsWeights = w,
              Q = Q, g1W = g1)                       # weighted POOLED targeting
  e <- .eif_from_tmle(fit, d$Y, d$A, w, gbound = g_floor)
  s <- .se_des(e$eif, d$strata, d$cluster, w, clustered = TRUE)
  eps_max  <- max(abs(fit$epsilon))
  diverged <- !is.finite(eps_max) || eps_max > 20    # R03 divergence guard
  list(
    row  = data.frame(method = method_label,
                      b  = if (diverged) NA_real_ else e$psi,
                      se = if (diverged) NA_real_ else s$se,
                      df = s$df, diverged = diverged,
                      stringsAsFactors = FALSE),
    drow = data.frame(method = method_label, eps = eps_max,
                      g_min = min(g1), g_max = max(g1),
                      V_eff = V_eff, stringsAsFactors = FALSE),
    eif  = e$eif
  )
}

# =====================================================================
# run_dr_arms(obs, learners, q_spec, g_spec, V_cf = 5L, g_floor = 0.05)
#   q_spec, g_spec in {"C","W"}: which design frame each nuisance learner sees.
#   Returns list(results, diagnostics) matching the engine row schema
#   (method, b, se, df [, diverged]); diagnostics$eif_faw is the FA-w EIF for
#   the deff_clust audit (mirrors how run_sim.R uses eif_fa).
# =====================================================================
run_dr_arms <- function(obs, learners, q_spec, g_spec,
                        V_cf = 5L, g_floor = .R14_G_FLOOR) {
  stopifnot(q_spec %in% c("C", "W"), g_spec %in% c("C", "W"))
  d <- obs
  n <- nrow(d); w <- d$weight
  fr <- make_dr_frames(d)
  WQ <- fr[[q_spec]]                       # frame the Q-learner sees
  Wg <- fr[[g_spec]]                       # frame the g-learner sees (independent)
  XAq <- data.frame(A = d$A, WQ)
  Xa1 <- data.frame(A = 1, WQ); Xa0 <- data.frame(A = 0, WQ)
  # tmle() requires a W argument, but with Q= and g1W= both supplied it does
  # NOT use it to fit nuisances (engine pattern, estimators.R lines 177/195);
  # pass the raw frame for determinism.
  W_targ <- fr$C

  ## ---- (1) FA-w: weighted SINGLE-FIT nuisances on the spec'd frames --------
  qf <- .sl(d$Y, XAq, weights = w, learners)
  gf <- .sl(d$A, Wg,  weights = w, learners)
  g1 <- pmin(pmax(gf(Wg), g_floor), 1 - g_floor)
  Qm <- pmin(pmax(cbind(qf(Xa0), qf(Xa1)), .R14_Q_CLAMP), 1 - .R14_Q_CLAMP)
  fa <- .dr_target(d, W_targ, Qm, g1, "FA-w", g_floor, NA_integer_)

  ## ---- shared PSU-level folds: built ONCE, used by BOTH CF arms (paired) ---
  fold  <- make_cf_folds(d$strata, d$cluster, V_cf)
  V_eff <- attr(fold, "V_eff")
  oof <- function(weighted) {
    Q0o <- Q1o <- g1o <- numeric(n)
    for (v in sort(unique(fold))) {
      tr <- which(fold != v); ho <- which(fold == v)
      # >>> THE ONE CONTESTED LINE (R07 device): OOF fits weighted vs not <<<
      wtr <- if (weighted) w[tr] else NULL
      qfv <- .sl(d$Y[tr], XAq[tr, ], weights = wtr, learners)
      gfv <- .sl(d$A[tr], Wg[tr, , drop = FALSE], weights = wtr, learners)
      Q1o[ho] <- qfv(Xa1[ho, ]); Q0o[ho] <- qfv(Xa0[ho, ])
      g1o[ho] <- gfv(Wg[ho, , drop = FALSE])
    }
    list(Q  = pmin(pmax(cbind(Q0o, Q1o), .R14_Q_CLAMP), 1 - .R14_Q_CLAMP),
         g1 = pmin(pmax(g1o, g_floor), 1 - g_floor))
  }
  nu_u <- oof(FALSE)                       # CF-u: unweighted OOF (paper default)
  nu_w <- oof(TRUE)                        # CF-w: weighted OOF (same folds)
  cfu <- .dr_target(d, W_targ, nu_u$Q, nu_u$g1, "CF-u", g_floor, V_eff)
  cfw <- .dr_target(d, W_targ, nu_w$Q, nu_w$g1, "CF-w", g_floor, V_eff)

  results <- rbind(fa$row, cfu$row, cfw$row)
  drow    <- rbind(fa$drow, cfu$drow, cfw$drow)
  rownames(results) <- NULL; rownames(drow) <- NULL
  list(results = results,
       diagnostics = list(eif_faw = fa$eif, drow = drow))
}
