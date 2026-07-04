# =====================================================================
# run_smoke.R  —  Phase-2 smoke checkpoint (plan-phase2.md sec 6)
# Two tiers:
#   GLM tier  -> model_type="simple"  (GLM correctly specified): validates the
#                DGP, the unequal-weight design, the EIF design effect, and the
#                SE plumbing; checks CF ~= single-fit (a correctness check).
#   SL  tier  -> model_type="complex" (needs gam/earth): the realistic
#                misspecification where cross-fitting (CF) should beat single-fit.
# Large per-rep output -> survey-tmle2-data ; summaries -> repo/results.
# =====================================================================

CODE <- normalizePath(Sys.getenv("SIM_CODE", "R"))   # absolute, so parallel workers can source it
source(file.path(CODE, "config.R"))
source(file.path(CODE, "dgp.R"))
source(file.path(CODE, "estimators.R"))
source(file.path(CODE, "diagnostics.R"))
suppressMessages({ library(parallel); library(survey) })

z <- qnorm(0.975)

# ---- one replication (top-level so parLapply can ship it) -------------------
smoke_one_rep <- function(i, pop, learners, model_type, seed_base, draw_args) {
  obs <- do.call(draw_sample, c(list(population = pop, sample_seed = seed_base + i,
                                     model_type = model_type), draw_args))
  est <- run_estimators(obs, learners = learners)
  dd  <- deff_clust(est$diagnostics$eif_fa, obs$strata, obs$cluster, obs$weight)
  list(results  = est$results,
       deff     = dd$deff_clust,
       icc_eif  = dd$icc_eif,
       checks   = attr(obs, "checks"),
       cf_veff  = est$diagnostics$cf_V_eff)
}

# ---- run one (scenario, model_type, learners) cell --------------------------
run_smoke_cell <- function(scenario, model_type, learners, n_reps,
                           dgp_args = list(), draw_args = list(),
                           n_cores = N_CORES_SMOKE, seed_base = SAMPLE_SEED_BASE,
                           label = NULL) {
  label <- label %||% sprintf("%s/%s/%s", scenario, model_type,
                              if (length(learners) == 1) learners else paste0("SL", length(learners)))
  cat(sprintf("\n========== smoke cell: %s  (%d reps) ==========\n", label, n_reps))
  pop <- do.call(make_population, c(list(scenario = scenario, model_type = model_type), dgp_args))
  print(pop)
  Psi <- pop$truth$psi
  pa  <- audit_population(pop)
  cat("population ICC:", sprintf("Y=%.3f A=%.3f C1=%.3f  mean_A=%.3f mean_Y=%.3f\n",
                                 pa$icc_Y, pa$icc_A, pa$icc_C1, pa$mean_A, pa$mean_Y))

  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "CODE", envir = environment())
  clusterEvalQ(cl, {
    source(file.path(CODE, "config.R"))
    source(file.path(CODE, "dgp.R"))
    source(file.path(CODE, "estimators.R"))
    source(file.path(CODE, "diagnostics.R"))
  })
  clusterSetRNGStream(cl, iseed = seed_base)
  reps <- parLapply(cl, seq_len(n_reps), smoke_one_rep,
                    pop = pop, learners = learners, model_type = model_type,
                    seed_base = seed_base, draw_args = draw_args)

  # ---- aggregate ----
  res <- do.call(rbind, Map(function(r, i) cbind(rep = i, r$results),
                            reps, seq_along(reps)))
  agg <- do.call(rbind, lapply(split(res, res$method), function(d) {
    crit <- qt(0.975, pmax(1, d$df))          # t-reference with each arm's design df
    p <- mean(abs(d$b - Psi) <= crit * d$se)
    data.frame(method = d$method[1],
               bias    = mean(d$b) - Psi,
               emp_sd  = sd(d$b),
               mean_se = mean(d$se),
               se_ratio = mean(d$se) / sd(d$b),
               coverage = p,
               mcse_cov = sqrt(p * (1 - p) / nrow(d)))
  }))
  rownames(agg) <- NULL
  ord <- c("Fully-Aware", "Fully-Aware-CF", "Fully-Aware-CV", "Partially-Aware", "Non-Aware")
  agg <- agg[order(match(agg$method, ord)), ]

  deff   <- mean(vapply(reps, `[[`, numeric(1), "deff"))
  icceif <- mean(vapply(reps, `[[`, numeric(1), "icc_eif"))
  chk    <- reps[[1]]$checks
  cf_veff <- reps[[1]]$cf_veff
  cat(sprintf("DEFF_clust(eif)=%.2f  icc_eif=%.3f  | n=%d sumw/N=%.3f minPSU=%d df=%d wCV=%.2f | CF V_eff=%d\n",
              deff, icceif, chk$n, chk$sumw_over_N, chk$min_psu_str, chk$df_design, chk$w_cv, cf_veff))
  print(round_df(agg))

  list(label = label, scenario = scenario, model_type = model_type, Psi = Psi,
       agg = agg, deff = deff, icc_eif = icceif, checks = chk, cf_veff = cf_veff,
       reps_raw = res)
}

round_df <- function(d, k = 3) { d[] <- lapply(d, function(x) if (is.numeric(x)) round(x, k) else x); d }

# ---- PASS/FAIL gate (MCSE-aware) --------------------------------------------
check_pass <- function(cell) {
  g <- function(m, col) cell$agg[cell$agg$method == m, col]
  und <- function(m) (g(m,"coverage") + 2*g(m,"mcse_cov")) < 0.92   # under-covers
  rec <- function(m) (g(m,"coverage") - 2*g(m,"mcse_cov")) > 0.90   # recovers
  list(
    DEFF_gt_1.3   = cell$deff > 1.3,
    sumw_ok       = abs(cell$checks$sumw_over_N - 1) < 0.05,
    no_lonely_psu = cell$checks$min_psu_str >= 2,
    FA_recovers   = rec("Fully-Aware"),
    CF_recovers   = if ("Fully-Aware-CF" %in% cell$agg$method) rec("Fully-Aware-CF") else NA,
    PA_undercover = und("Partially-Aware"),
    NA_undercover = und("Non-Aware"),
    FA_unbiased   = abs(g("Fully-Aware","bias")) < 2 * g("Fully-Aware","emp_sd") / sqrt(nrow(cell$reps_raw)/5)
  )
}

# =====================================================================
# MAIN  (edit DGP_ARGS once the RC-2 sweep picks the variance components)
# =====================================================================
if (sys.nframe() == 0L) {
  # Winning DGP config from Phase-2 tuning (clustering + informative weights, stable).
  # These match the dgp.R defaults; stated explicitly for the record.
  DGP_ARGS  <- list(sigma2_C = 0.3, sigma2_A = 0.8, sigma2_Y = 1.5, beta_strat = 1.0,
                    alpha_g = log(1.3), p_treat_target = 0.35, truth_M = 2e6L)
  DRAW_ARGS <- list(alpha_strat = 2.0)   # n0-oversampling; scenario sets m_h/n0 defaults
  GLM_LIB <- "SL.glm"
  SL_LIB  <- c("SL.glm", "SL.gam", "SL.earth")
  TIER    <- Sys.getenv("SMOKE_TIER", "all")     # "glm", "sl", or "all"
  N_GLM   <- 50L; N_SL <- 50L

  cells <- list()
  if (TIER %in% c("glm", "all")) {
    cells$std_glm <- run_smoke_cell("standard", "simple", GLM_LIB, N_GLM, DGP_ARGS, DRAW_ARGS)
    cells$r1_glm  <- run_smoke_cell("R1",       "simple", GLM_LIB, N_GLM, DGP_ARGS, DRAW_ARGS)
  }
  if (TIER %in% c("sl", "all")) {
    cells$std_sl  <- run_smoke_cell("standard", "complex", SL_LIB, N_SL, DGP_ARGS, DRAW_ARGS)
    cells$r1_sl   <- run_smoke_cell("R1",       "complex", SL_LIB, N_SL, DGP_ARGS, DRAW_ARGS)
  }

  cat("\n===================== PASS/FAIL =====================\n")
  for (nm in names(cells)) { cat("\n--", cells[[nm]]$label, "--\n"); print(unlist(check_pass(cells[[nm]]))) }

  saveRDS(cells, file.path(DATA_RESULTS, paste0("smoke_cells_", TIER, ".rds")))
  summ <- do.call(rbind, lapply(cells, function(c) cbind(cell = c$label, DEFF = round(c$deff,2), round_df(c$agg))))
  write.csv(summ, file.path(RESULTS_DIR, paste0("smoke_summary_", TIER, ".csv")), row.names = FALSE)
  cat("\nSaved:", file.path(RESULTS_DIR, paste0("smoke_summary_", TIER, ".csv")), "\n")
}
