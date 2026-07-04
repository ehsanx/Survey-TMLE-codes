# =====================================================================
# render_fig_deployable_ladder.R
#   Web Appendix D figure: the unified coverage ladder L1-L7 x design{A,B} x
#   estimator{TMLE,AIPW} x protocol{single-fit, cluster-CV, cross-fit}, stitched
#   from the locked runs (+ R22 when present). Lines BREAK across unmeasured
#   rungs; cells never run are drawn as N/A.
#
#   Sources (all locked CSVs):
#     L1-L4 TMLE  : results/sim_full_summary.csv
#     L1-L4 AIPW  : results/arc/R15_aipw_benchmark_summary.csv  (se_type='lin'; SF/CF)
#     L5    TMLE  : results/arc/R20_cv_vs_cf_nondonsker_summary.csv
#     L6,L7 both  : results/arc/R21_deployable_cvcf_summary.csv
#     L2/L3 AIPW-CV + L5 AIPW-{SF,CV,CF} : results/arc/R22_aipw_ladder_complete_summary.csv
#         (the alpha completion; if ABSENT these cells render as N/A -- the
#          gap-free version appears automatically once R22 lands)
#   L1/L4 internal-CV is N/A BY CONSTRUCTION (a lone learner has nothing to
#   weight -> internal CV degenerates to single-fit).
#
#   Output: Writing/Appendix/images/fig_deployable_ladder.png  (coverage)
#           Writing/Appendix/images/fig_deployable_ladder_seratio.png  (se/sd)
# =====================================================================
RES <- "results"; ARC <- RES   # public: summaries are flat under results/
IMG <- "outputs/figures"; if (!dir.exists(IMG)) dir.create(IMG, recursive = TRUE)
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL

simfull <- rd(file.path(RES, "sim_full_summary.csv"))
r15     <- rd(file.path(ARC, "R15_aipw_benchmark_summary.csv"))
r20     <- rd(file.path(ARC, "R20_cv_vs_cf_nondonsker_summary.csv"))
r21     <- rd(file.path(ARC, "R21_deployable_cvcf_summary.csv"))
r22     <- rd(file.path(ARC, "R22_aipw_ladder_complete_summary.csv"))   # may be NULL (pending)
stopifnot(!is.null(simfull), !is.null(r15), !is.null(r20), !is.null(r21))
if (!is.null(r15) && "se_type" %in% names(r15)) {
  r15 <- r15[r15$se_type == (if ("lin" %in% r15$se_type) "lin" else r15$se_type[1]), ]
}
cat(if (is.null(r22)) "R22 ABSENT -> AIPW@L5 + AIPW-CV@L2/L3 render as N/A (pending completion)\n"
    else "R22 present -> figure is gap-free\n")

rung_map <- c(L1_param=1, L2_smooth=2, L3_adaptive=3, L4_aggressive=4,
              L5_nondonsker=5, L6_deployable=6, L7_hal=7)
des_map  <- c(standard="A", R1="B")
prot_tmle <- c("Fully-Aware"="SF", "Fully-Aware-CV"="CV", "Fully-Aware-CF"="CF")
prot_aipw <- c("AIPW-SF"="SF", "AIPW-CV"="CV", "AIPW-CF"="CF")
grab <- function(df, est, pmap) {
  if (is.null(df)) return(NULL)
  d <- df[df$method %in% names(pmap), ]
  if (!nrow(d)) return(NULL)
  data.frame(rung = rung_map[d$rung], design = des_map[d$scenario], estimator = est,
             protocol = unname(pmap[d$method]), coverage = d$coverage,
             se_ratio = d$se_ratio, stringsAsFactors = FALSE)
}
A <- rbind(
  grab(simfull, "TMLE", prot_tmle),                              # L1-L4 TMLE
  grab(r15,     "AIPW", prot_aipw),                              # L1-L4 AIPW (SF/CF)
  grab(r20,     "TMLE", prot_tmle),                              # L5 TMLE
  grab(r21[r21$rung=="L6_deployable",], "TMLE", prot_tmle),
  grab(r21[r21$rung=="L6_deployable",], "AIPW", prot_aipw),
  grab(r21[r21$rung=="L7_hal",],        "TMLE", prot_tmle),
  grab(r21[r21$rung=="L7_hal",],        "AIPW", prot_aipw),
  grab(r22, "AIPW", prot_aipw)                                   # L2/L3 AIPW-CV + L5 AIPW-{SF,CV,CF}; NULL-safe
)
A <- A[!is.na(A$rung), ]
A <- A[!duplicated(A[, c("rung","design","estimator","protocol")]), ]   # locked runs take precedence

# ---- plotting ---------------------------------------------------------------
col <- c(SF="#E69F00", CV="#CC79A7", CF="#0072B2")          # orange / purple / blue (Okabe-Ito)
pch <- c(SF=15, CV=17, CF=19); lty <- c(SF=2, CV=3, CF=1); lwd <- c(SF=2.0, CV=2.0, CF=2.8)
RLAB <- c("L1\nGLM","L2\n+smooth","L3\n+RF","L4\ndeep RF",
          "L5\n(3 lrn)","L6\n(5 lrn)","L7\n(6 lrn)")

panel <- function(est, des, yvar, ylim, href, ylab, main) {
  plot(NA, xlim=c(0.7,7.3), ylim=ylim, xaxt="n", xlab="", ylab=ylab, main=main, cex.main=1.0)
  axis(1, at=1:7, labels=RLAB, padj=0.5, cex.axis=0.7)
  abline(v=4.5, col="grey70", lty=2)
  abline(h=href, col="grey55", lwd=1)
  rect(0.7, ylim[1], 4.5, ylim[2], col="#00000007", border=NA)
  s <- A[A$estimator==est & A$design==des, ]
  for (p in c("SF","CV","CF")) {
    sp <- s[s$protocol==p, ]; sp <- sp[order(sp$rung), ]
    if (!nrow(sp)) next
    g <- cumsum(c(1, as.integer(diff(sp$rung) != 1)))       # break lines across gaps
    for (run in split(sp, g))
      if (nrow(run) > 1) lines(run$rung, run[[yvar]], col=col[p], lty=lty[p], lwd=lwd[p])
    points(sp$rung, sp[[yvar]], col=col[p], pch=pch[p], cex=1.25, lwd=2)
  }
}

draw <- function(yvar, ylim, href, ylab, file) {
  png(file, width=2100, height=1720, res=210)
  op <- par(mfrow=c(2,2), mar=c(3.8,4.0,2.6,0.8), mgp=c(2.4,0.7,0), las=1, oma=c(4.2,0,1.2,0))
  panel("TMLE","A", yvar, ylim, href, ylab, "TMLE - Design A")
  legend("bottomright", c("single-fit","cluster/internal-CV","cross-fit (CF)"),
         col=col, lty=lty, pch=pch, lwd=lwd, bty="n", cex=0.8)
  panel("TMLE","B", yvar, ylim, href, ylab, "TMLE - Design B")
  panel("AIPW","A", yvar, ylim, href, ylab, "AIPW - Design A")
  panel("AIPW","B", yvar, ylim, href, ylab, "AIPW - Design B")
  mtext("lone-learner complexity ladder            |            deployable ensembles",
        side=1, outer=TRUE, cex=0.72, col="grey35", line=1.3)
  na_note <- if (is.null(r22))
    "L5/L6/L7: GLM,MARS,deep RF (+glmnet,XGBoost; +HAL).  AIPW @ L5 and AIPW cluster-CV @ L2/L3 pending (shown absent); L1/L4 cluster-CV N/A by construction."
  else
    "L5/L6/L7: GLM,MARS,deep RF (+glmnet,XGBoost; +HAL).  L1/L4 cluster-CV N/A by construction (a lone learner has nothing to weight)."
  mtext(na_note, side=1, outer=TRUE, cex=0.58, col="grey45", line=2.6)
  par(op); dev.off(); cat("wrote", file, "\n")
}

draw("coverage", c(0,1), 0.95, "Empirical 95% CI coverage",
     file.path(IMG, "fig_deployable_ladder.png"))
draw("se_ratio", c(0,1.6), 1.0, "mean SE / empirical SD",
     file.path(IMG, "fig_deployable_ladder_seratio.png"))
