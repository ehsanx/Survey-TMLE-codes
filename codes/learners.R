# =====================================================================
# learners.R  —  the SuperLearner library LADDER (Donsker-complexity design)
# See plan-phase2.md sec 11. The ladder is ordered by nuisance-class COMPLEXITY
# (entropy / Donsker), NOT learner count: smooth/parametric (single-fit valid)
# -> deep trees / boosting that overfit in-sample (single-fit under-covers, CF
# rescues). The aggressive rungs use a THIN/SINGLE data-adaptive learner so that
# SuperLearner's internal CV cannot hide overfitting behind a smooth ensemble
# member (which would defeat the demonstration).
# =====================================================================

suppressMessages({ library(SuperLearner) })

# ---- custom aggressive (overfitting / non-Donsker) wrappers -----------------
# Default random forest, single-threaded (L3 "realistic ensemble" member).
SL.ranger.t1 <- function(Y, X, newX, family, obsWeights, ...) {
  SuperLearner::SL.ranger(Y = Y, X = X, newX = newX, family = family,
                          obsWeights = obsWeights, num.threads = 1)
}

# Deep random forest: full-depth trees (min.node.size=1) using all features
# (mtry=p) -> interpolates the training data -> strong in-sample overfit.
SL.ranger.deep <- function(Y, X, newX, family, obsWeights, ...) {
  SuperLearner::SL.ranger(Y = Y, X = X, newX = newX, family = family,
                          obsWeights = obsWeights,
                          num.trees = 500, min.node.size = 1, mtry = ncol(X),
                          num.threads = 1)   # 1 thread: parLapply already parallelizes reps
}

# ---- the ladder (named; each entry is a SuperLearner SL.library vector) ------
# Validated by the local pilot (40 reps, standard/complex): single-fit Fully-Aware
# coverage 0.95 -> 1.00 -> 0.85 -> 0.175 across L1->L4 while cross-fitting (CF)
# holds at ~0.95-1.00 -- empirically tracing the Donsker boundary. Boosting was
# dropped: L4 already gives a decisive (and cheaper) demonstration, deep boosting
# was redundant + the most expensive rung. (SL.xgb.deep wrapper kept available
# below, commented, if a non-RF robustness rung is wanted later.)
SL_LADDER <- list(
  L1_param      = "SL.glm",                              # Donsker (parametric) -> single-fit valid (0.95)
  L2_smooth     = c("SL.glm", "SL.gam", "SL.earth"),     # smooth semiparametric -> ~Donsker, valid (1.00)
  L3_adaptive   = c("SL.glm", "SL.earth", "SL.ranger.t1"), # realistic ensemble + default RF -> boundary (0.85)
  L4_aggressive = "SL.ranger.deep"                        # deep RF (single) -> non-Donsker, single-fit 0.175
)

# packages the ladder needs on each parallel worker
SL_LADDER_PKGS <- c("SuperLearner", "gam", "earth", "ranger")

# Optional non-RF robustness rung (errored in the pilot via the SL.xgboost path;
# not needed given L4). Uncomment + add to SL_LADDER to include:
# SL.xgb.deep <- function(Y, X, newX, family, obsWeights, ...)
#   SuperLearner::SL.xgboost(Y = Y, X = X, newX = newX, family = family,
#                            obsWeights = obsWeights, ntrees = 500, max_depth = 8,
#                            shrinkage = 0.1, minobspernode = 1)
