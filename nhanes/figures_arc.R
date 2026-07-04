# =====================================================================
# nhanes/figures_arc.R  —  figures from the aggregated ARC results
# Reads results/nhanes_combined.rds and draws, per example:
#   (1) 5-arm forest plot at a representative rung (default L3_adaptive)
#   (2) Fully-Aware (single-fit) vs Fully-Aware-CF estimate +/- CI across L1->L4
#       (the APPLIED echo of Figure 1: single-fit drift/instability vs CF stability)
#   (3) propensity overlap histograms (single-fit vs cross-fitted) at the deep rung
# Outputs to nhanes/figures/.
# =====================================================================

RES <- Sys.getenv("NH_RESULTS", "nhanes/nhanes_output/results")
FIG <- Sys.getenv("NH_FIG", "nhanes/figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
cb  <- readRDS(file.path(RES, "nhanes_combined.rds"))
S   <- cb$summary; cells <- cb$cells
EX  <- unique(S$example); RUNGS <- c("L1_param","L2_smooth","L3_adaptive","L4_aggressive")
ARMS <- c("Non-Aware","Partially-Aware","Fully-Aware","Fully-Aware-CV","Fully-Aware-CF")
PCH  <- c(1, 2, 0, 5, 19)
getcell <- function(ex, rung) Filter(function(o) o$example==ex && o$rung==rung, cells)[[1]]

## (1) 5-arm forest per example at a representative rung
forest_rung <- "L3_adaptive"
png(file.path(FIG, "forest_all_rungL3.png"), 2600, 2000, res = 300); par(mfrow = c(2,2))
for (ex in EX) {
  d <- S[S$example==ex & S$rung==forest_rung, ]; d <- d[match(ARMS, d$method), ]; d <- d[!is.na(d$method),]
  k <- nrow(d); yp <- rev(seq_len(k)); xl <- range(c(d$lcl,d$ucl,0),na.rm=TRUE)+c(-1,1)*0.05*diff(range(c(d$lcl,d$ucl,0)))
  par(mar=c(4,8.5,2.4,1)); plot(NA,xlim=xl,ylim=c(.5,k+.5),yaxt="n",xlab="Risk difference (95% CI)",ylab="",
       main=paste0(ex,": ",d$label[1]),cex.main=.9)
  abline(v=0,lty=2,col="grey55"); axis(2,at=yp,labels=d$method,las=1,cex.axis=.75)
  for(i in seq_len(k)){segments(d$lcl[i],yp[i],d$ucl[i],yp[i],lwd=2)
    points(d$b[i],yp[i],pch=PCH[match(d$method[i],ARMS)],cex=1.2,lwd=1.6,bg="white")}
}
dev.off()

## (2) Fully-Aware (single-fit) vs Fully-Aware-CF across the ladder, per example
png(file.path(FIG, "ladder_FA_vs_CF.png"), 2600, 2000, res = 300); par(mfrow = c(2,2))
for (ex in EX) {
  fa <- S[S$example==ex & S$method=="Fully-Aware", ];  fa <- fa[match(RUNGS,fa$rung),]
  cf <- S[S$example==ex & S$method=="Fully-Aware-CF", ]; cf <- cf[match(RUNGS,cf$rung),]
  x <- seq_along(RUNGS); yl <- range(c(fa$lcl,fa$ucl,cf$lcl,cf$ucl,0),na.rm=TRUE)
  par(mar=c(4.2,4.2,2.4,1)); plot(NA,xlim=c(.8,4.2),ylim=yl,xaxt="n",xlab="SL library (Donsker ladder)",
       ylab="Risk difference (95% CI)",main=paste0(ex,": ",fa$label[1]),cex.main=.9)
  axis(1,at=x,labels=c("L1","L2","L3","L4")); abline(h=0,lty=2,col="grey55")
  off <- 0.08
  segments(x-off,fa$lcl,x-off,fa$ucl,lwd=2,col="grey35"); points(x-off,fa$b,pch=1,lwd=1.6,col="grey35")
  segments(x+off,cf$lcl,x+off,cf$ucl,lwd=2); points(x+off,cf$b,pch=19)
  if(ex==EX[1]) legend("topleft",c("Fully-Aware (single-fit)","Fully-Aware-CF"),pch=c(1,19),
                       col=c("grey35","black"),bty="n",cex=.75)
}
dev.off()

## (3) propensity overlap (single-fit vs CF) at the deep rung, per example
png(file.path(FIG, "overlap_L4.png"), 2600, 2000, res = 300); par(mfrow = c(2,2))
for (ex in EX) {
  o <- tryCatch(getcell(ex, "L4_aggressive"), error=function(e) tryCatch(getcell(ex,"L3_adaptive"),error=function(e2) NULL))
  if (is.null(o)) next
  par(mar=c(4,4,2.4,1))
  hist(o$diagnostics$g_fa, breaks=40, col=rgb(.6,.6,.6,.5), border=NA, xlim=c(0,1), freq=FALSE,
       main=paste0(ex,": propensity overlap"), xlab="estimated P(exposed | covariates)", cex.main=.9)
  hist(o$diagnostics$g_cf, breaks=40, col=rgb(0,0,0,.35), border=NA, add=TRUE, freq=FALSE)
  legend("top",c("single-fit","cross-fitted"),fill=c(rgb(.6,.6,.6,.5),rgb(0,0,0,.35)),border=NA,bty="n",cex=.7)
}
dev.off()

cat("Figures written to", FIG, ":\n - forest_all_rungL3.png\n - ladder_FA_vs_CF.png\n - overlap_L4.png\n")
