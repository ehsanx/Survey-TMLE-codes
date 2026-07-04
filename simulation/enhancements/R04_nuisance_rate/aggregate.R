# =====================================================================
# aggregate.R  —  combine R04_nuisance_rate per-task RDS into the rate table
# Usage:  Rscript aggregate.R   (reads R04_OUT, writes results/R04_nuisance_rate_summary.csv)
#
# Produces the per-rung x per-scenario nuisance-rate table:
#   mean L2 error for Q (u-integrated target) and g (uA-integrated target),
#   their product, and product*sqrt(m) -- the DEMONSTRATION object. Expected:
#   prod_sqrtm shrinks L1->L3 (rate met), blows up / stays large at L4 (violated).
# =====================================================================
CODE <- Sys.getenv("SIM_CODE", "")
if (!nzchar(CODE) || !file.exists(file.path(CODE, "config.R"))) {
  cand <- c("R", "R")
  hit  <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("cannot find config.R; set SIM_CODE or run from repo root")
  CODE <- hit[1]
}
source(file.path(CODE, "config.R"))

OUT <- Sys.getenv("R04_OUT", file.path(DATA_ROOT, "arc_runs", "R04_nuisance_rate"))
files <- list.files(OUT, pattern = "^r04_.*\\.rds$", full.names = TRUE)  # skip SMOKE_*/manifest_*
if (!length(files)) stop("no r04_*.rds found in ", OUT)
cat("aggregating", length(files), "task files from", OUT, "\n")

per_rep <- do.call(rbind, lapply(files, function(f) readRDS(f)$per_rep))

# ladder order for sorting
rung_lv <- c("L1_param", "L2_smooth", "L3_adaptive", "L4_aggressive")
per_rep$rung <- factor(per_rep$rung, levels = rung_lv)

agg <- do.call(rbind, by(per_rep, list(per_rep$scenario, per_rep$rung), function(d) {
  if (is.null(d) || !nrow(d)) return(NULL)
  data.frame(
    scenario = d$scenario[1], rung = as.character(d$rung[1]), n_reps = nrow(d),
    mean_n = mean(d$n), mean_m_psu = mean(d$m_psu),
    # PRIMARY (u-integrated Q, uA-integrated g) -- the rate object
    mean_eQ = mean(d$eQ_int), mean_eg = mean(d$eg_int),
    mean_prod = mean(d$prod_int),
    mean_prod_sqrtm = mean(d$prod_int_sqrtm),
    sd_prod_sqrtm   = sd(d$prod_int_sqrtm),
    mean_prod_sqrtn = mean(d$prod_int_sqrtn),
    # SECONDARY (realized-u Q incl uY, realized-u g incl uA)
    mean_eQ_real = mean(d$eQ_real), mean_eg_real = mean(d$eg_real),
    mean_prod_real = mean(d$prod_real),
    mean_prod_real_sqrtm = mean(d$prod_real_sqrtm),
    truth_join = paste(unique(d$truth_join), collapse = ","),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$scenario, match(agg$rung, rung_lv)), ]
rownames(agg) <- NULL

# decision-rule helper: the ratio of prod_sqrtm at each rung vs L1 (per scenario)
agg$prod_sqrtm_vs_L1 <- ave(agg$mean_prod_sqrtm, agg$scenario,
  FUN = function(v) v / v[1])

num <- sapply(agg, is.numeric); ap <- agg; ap[num] <- round(ap[num], 5)
print(ap, row.names = FALSE)

csv <- file.path(RESULTS_DIR, "arc", "R04_nuisance_rate_summary.csv")
if (!dir.exists(dirname(csv))) dir.create(dirname(csv), recursive = TRUE, showWarnings = FALSE)
write.csv(ap, csv, row.names = FALSE)
saveRDS(list(per_rep = per_rep, summary = agg),
        file.path(OUT, "R04_nuisance_rate_combined.rds"))
cat("\nSaved:\n -", csv, "\n -", file.path(OUT, "R04_nuisance_rate_combined.rds"), "\n")
