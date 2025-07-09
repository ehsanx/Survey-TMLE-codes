library(survey)
library(rsimsum)
library(parallel)

# Step 1: Define the data simulation function
create.data <- function(n_target_sample = 1000, te_log_odds = 1.5, seed = NULL, model_type = "complex"){
  if (!is.null(seed)) { set.seed(seed) }
  N_pop <- 50000 
  C1 <- rnorm(N_pop, 0, 1); C2 <- rnorm(N_pop, 0, 1); C3 <- rnorm(N_pop, 0, 1); C4 <- rnorm(N_pop, 0, 1)
  strata <- rep(1:4, each = N_pop / 4)
  cluster <- rep(1:80, each = N_pop / 80)
  sampling_prob <- plogis(-2.5 + 0.6 * C1)
  selected_indices <- sample(1:N_pop, size = n_target_sample, prob = sampling_prob, replace = FALSE)
  C1_s <- C1[selected_indices]; C2_s <- C2[selected_indices]; C3_s <- C3[selected_indices]; C4_s <- C4[selected_indices]
  strata_s <- strata[selected_indices]; cluster_s <- cluster[selected_indices]
  weights_s <- 1 / sampling_prob[selected_indices]
  pscore <- plogis(-1 + log(1.75) * (C1_s + C2_s + C3_s + C4_s))
  A <- rbinom(n_target_sample, 1, pscore)
  get.prob_Y <- function(A_in, C1_in, C2_in, C3_in, C4_in, TE = te_log_odds) {
    plogis(-2 + TE * A_in + 0.5 * (C1_in + C2_in + C3_in + C4_in))
  }
  observed_prob_Y <- get.prob_Y(A, C1_s, C2_s, C3_s, C4_s)
  Y <- rbinom(n_target_sample, 1, observed_prob_Y)
  if (model_type == "complex") {
    L1 <- exp(C1_s / 2); L2 <- C2_s / (1 + exp(C1_s)) + 10; L3 <- (C1_s * C3_s / 25 + 0.6)^3; L4 <- (C2_s + C4_s + 20)^2
  } else {
    L1 <- C1_s; L2 <- C2_s; L3 <- C3_s; L4 <- C4_s
  }
  prob_Y1_pop <- get.prob_Y(1, C1, C2, C3, C4); prob_Y0_pop <- get.prob_Y(0, C1, C2, C3, C4)
  true_ate_rd <- mean(prob_Y1_pop) - mean(prob_Y0_pop)
  obs <- data.frame(L1, L2, L3, L4, A, Y, strata = strata_s, cluster = cluster_s, weight = weights_s)
  return(list(observed = obs, ate = true_ate_rd))
}

# --- Main simulation function ---
run_one_simulation_svyglm <- function(i, n_sample, model_type) {
  sim_data_list <- create.data(n_target_sample = n_sample, seed = i, model_type = model_type)
  observed_data <- sim_data_list$observed
  true_ate <- sim_data_list$ate
  
  model_formula <- as.formula(Y ~ A + L1 + L2 + L3 + L4)
  data_a1 <- observed_data; data_a1$A <- 1
  data_a0 <- observed_data; data_a0$A <- 0
  
  # --- Method 1: Fully-Aware Analysis ---
  design_full <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = observed_data)
  fit_full_logit <- svyglm(model_formula, design = design_full, family = quasibinomial())
  mm_a1 <- model.matrix(model_formula, data = data_a1)
  mm_a0 <- model.matrix(model_formula, data = data_a0)
  Q_1_W_full <- plogis(mm_a1 %*% coef(fit_full_logit))
  Q_0_W_full <- plogis(mm_a0 %*% coef(fit_full_logit))
  psi_m_rd_full <- svymean(~I(Q_1_W_full - Q_0_W_full), design_full)[1]
  
  rep_design_full <- as.svrepdesign(design_full, type = "bootstrap", replicates = 200)
  boot_var_full <- withReplicates(rep_design_full, 
                                  theta = function(weights, data) {
                                    boot_design <- svydesign(ids = ~cluster, strata = ~strata, weights = weights, data = data)
                                    fit_boot <- svyglm(model_formula, design = boot_design, family = quasibinomial())
                                    Q_1_W_boot <- plogis(mm_a1 %*% coef(fit_boot))
                                    Q_0_W_boot <- plogis(mm_a0 %*% coef(fit_boot))
                                    return(weighted.mean(Q_1_W_boot - Q_0_W_boot, weights))
                                  })
  se_m_rd_full <- sqrt(attr(boot_var_full, "var"))
  res1a <- data.frame(model = model_type, method = "Fully-Aware-svyGLM", estimand = "Marginal-RD", true = true_ate, b = psi_m_rd_full, se = se_m_rd_full)
  
  summary_logit_full <- summary(fit_full_logit)
  log_or_full <- summary_logit_full$coefficients['A', 'Estimate']; se_log_or_full <- summary_logit_full$coefficients['A', 'Std. Error']
  psi_c_or_full <- exp(log_or_full); se_c_or_full <- psi_c_or_full * se_log_or_full
  res1b <- data.frame(model = model_type, method = "Fully-Aware-svyGLM", estimand = "Conditional-OR", true = NA, b = psi_c_or_full, se = se_c_or_full)
  
  fit_full_linear <- svyglm(model_formula, design = design_full, family = gaussian(link="identity"))
  summary_linear_full <- summary(fit_full_linear)
  psi_c_rd_full <- summary_linear_full$coefficients['A', 'Estimate']; se_c_rd_full <- summary_linear_full$coefficients['A', 'Std. Error']
  res1c <- data.frame(model = model_type, method = "Fully-Aware-svyGLM", estimand = "Conditional-RD", true = true_ate, b = psi_c_rd_full, se = se_c_rd_full)
  
  fit_full_log <- svyglm(model_formula, design = design_full, family = quasipoisson(link="log"))
  summary_log_full <- summary(fit_full_log)
  log_rr_full <- summary_log_full$coefficients['A', 'Estimate']; se_log_rr_full <- summary_log_full$coefficients['A', 'Std. Error']
  psi_c_rr_full <- exp(log_rr_full); se_c_rr_full <- psi_c_rr_full * se_log_rr_full
  res1d <- data.frame(model = model_type, method = "Fully-Aware-svyGLM", estimand = "Conditional-RR", true = NA, b = psi_c_rr_full, se = se_c_rr_full)
  
  # --- Method 2: Partially-Aware Analysis ---
  design_partial <- svydesign(ids = ~1, weights = ~weight, data = observed_data)
  fit_partial_logit <- svyglm(model_formula, design = design_partial, family = quasibinomial())
  Q_1_W_partial <- plogis(mm_a1 %*% coef(fit_partial_logit))
  Q_0_W_partial <- plogis(mm_a0 %*% coef(fit_partial_logit))
  psi_m_rd_partial <- svymean(~I(Q_1_W_partial - Q_0_W_partial), design_partial)[1]
  
  rep_design_partial <- as.svrepdesign(design_partial, type = "bootstrap", replicates = 200)
  boot_var_partial <- withReplicates(rep_design_partial, 
                                     theta = function(weights, data) {
                                       boot_design <- svydesign(ids = ~1, weights = weights, data = data)
                                       fit_boot <- svyglm(model_formula, design = boot_design, family = quasibinomial())
                                       Q_1_W_boot <- plogis(mm_a1 %*% coef(fit_boot))
                                       Q_0_W_boot <- plogis(mm_a0 %*% coef(fit_boot))
                                       return(weighted.mean(Q_1_W_boot - Q_0_W_boot, weights))
                                     })
  se_m_rd_partial <- sqrt(attr(boot_var_partial, "var"))
  res2a <- data.frame(model = model_type, method = "Partially-Aware-svyGLM", estimand = "Marginal-RD", true = true_ate, b = psi_m_rd_partial, se = se_m_rd_partial)
  
  summary_logit_partial <- summary(fit_partial_logit)
  log_or_partial <- summary_logit_partial$coefficients['A', 'Estimate']; se_log_or_partial <- summary_logit_partial$coefficients['A', 'Std. Error']
  psi_c_or_partial <- exp(log_or_partial); se_c_or_partial <- psi_c_or_partial * se_log_or_partial
  res2b <- data.frame(model = model_type, method = "Partially-Aware-svyGLM", estimand = "Conditional-OR", true = NA, b = psi_c_or_partial, se = se_c_or_partial)
  
  fit_partial_linear <- svyglm(model_formula, design = design_partial, family = gaussian(link="identity"))
  summary_linear_partial <- summary(fit_partial_linear)
  psi_c_rd_partial <- summary_linear_partial$coefficients['A', 'Estimate']; se_c_rd_partial <- summary_linear_partial$coefficients['A', 'Std. Error']
  res2c <- data.frame(model = model_type, method = "Partially-Aware-svyGLM", estimand = "Conditional-RD", true = true_ate, b = psi_c_rd_partial, se = se_c_rd_partial)
  
  fit_partial_log <- svyglm(model_formula, design = design_partial, family = quasipoisson(link="log"))
  summary_log_partial <- summary(fit_partial_log)
  log_rr_partial <- summary_log_partial$coefficients['A', 'Estimate']; se_log_rr_partial <- summary_log_partial$coefficients['A', 'Std. Error']
  psi_c_rr_partial <- exp(log_rr_partial); se_c_rr_partial <- psi_c_rr_partial * se_log_rr_partial
  res2d <- data.frame(model = model_type, method = "Partially-Aware-svyGLM", estimand = "Conditional-RR", true = NA, b = psi_c_rr_partial, se = se_c_rr_partial)
  
  # --- Method 3: Non-Aware Analysis ---
  design_none <- svydesign(ids = ~1, weights = ~1, data = observed_data)
  fit_none_logit <- svyglm(model_formula, design = design_none, family = quasibinomial())
  Q_1_W_none <- plogis(mm_a1 %*% coef(fit_none_logit))
  Q_0_W_none <- plogis(mm_a0 %*% coef(fit_none_logit))
  psi_m_rd_none <- mean(Q_1_W_none - Q_0_W_none)
  
  rep_design_none <- as.svrepdesign(design_none, type = "bootstrap", replicates = 200)
  boot_var_none <- withReplicates(rep_design_none, 
                                  theta = function(weights, data) {
                                    boot_design <- svydesign(ids = ~1, weights = weights, data = data)
                                    fit_boot <- svyglm(model_formula, design = boot_design, family = quasibinomial())
                                    Q_1_W_boot <- plogis(mm_a1 %*% coef(fit_boot))
                                    Q_0_W_boot <- plogis(mm_a0 %*% coef(fit_boot))
                                    return(weighted.mean(Q_1_W_boot - Q_0_W_boot, weights))
                                  })
  se_m_rd_none <- sqrt(attr(boot_var_none, "var"))
  res3a <- data.frame(model = model_type, method = "Non-Aware-svyGLM", estimand = "Marginal-RD", true = true_ate, b = psi_m_rd_none, se = se_m_rd_none)
  
  summary_logit_none <- summary(fit_none_logit)
  log_or_none <- summary_logit_none$coefficients['A', 'Estimate']; se_log_or_none <- summary_logit_none$coefficients['A', 'Std. Error']
  psi_c_or_none <- exp(log_or_none); se_c_or_none <- psi_c_or_none * se_log_or_none
  res3b <- data.frame(model = model_type, method = "Non-Aware-svyGLM", estimand = "Conditional-OR", true = NA, b = psi_c_or_none, se = se_c_or_none)
  
  fit_none_linear <- svyglm(model_formula, design = design_none, family = gaussian(link="identity"))
  summary_linear_none <- summary(fit_none_linear)
  psi_c_rd_none <- summary_linear_none$coefficients['A', 'Estimate']; se_c_rd_none <- summary_linear_none$coefficients['A', 'Std. Error']
  res3c <- data.frame(model = model_type, method = "Non-Aware-svyGLM", estimand = "Conditional-RD", true = true_ate, b = psi_c_rd_none, se = se_c_rd_none)
  
  fit_none_log <- svyglm(model_formula, design = design_none, family = quasipoisson(link="log"))
  summary_log_none <- summary(fit_none_log)
  log_rr_none <- summary_log_none$coefficients['A', 'Estimate']; se_log_rr_none <- summary_log_none$coefficients['A', 'Std. Error']
  psi_c_rr_none <- exp(log_rr_none); se_c_rr_none <- psi_c_rr_none * se_log_rr_none
  res3d <- data.frame(model = model_type, method = "Non-Aware-svyGLM", estimand = "Conditional-RR", true = NA, b = psi_c_rr_none, se = se_c_rr_none)
  
  # --- Return all 12 result data frames ---
  return(rbind(res1a, res1b, res1c, res1d, 
               res2a, res2b, res2c, res2d,
               res3a, res3b, res3c, res3d))
}

# --- Simulation Parameters ---
RUN_IN_PARALLEL <- TRUE 
n_reps <- 1000           
sample_size <- 500
model_choices <- c("complex", "simple") 

# --- Main simulation loop ---
all_model_results <- list()

for (current_model in model_choices) {
  cat("=========================================================\n")
  cat("Running: Model Type =", toupper(current_model), "with svyGLM\n")
  
  results_dir <- file.path("data", "intermediate_results_svyglm_boot", current_model)
  if (!dir.exists(results_dir)) {dir.create(results_dir, recursive = TRUE)}
  
  if (RUN_IN_PARALLEL) {
    # --- Parallel Path (for Linux cluster) ---
    no_cores <- detectCores() - 1 
    cl <- makeCluster(no_cores)
    clusterEvalQ(cl, library(survey))
    clusterExport(cl, varlist=c("create.data", "run_one_simulation_svyglm", "sample_size", "current_model"))
    cat("=========================================================\n")
    cat("Running", n_reps, "simulations in PARALLEL on", no_cores, "cores...\n")
    results_list <- parLapply(cl, 1:n_reps, function(i) {
      try(run_one_simulation_svyglm(i, n_sample = sample_size, model_type = current_model))
    })
    stopCluster(cl)
  } else {
    # --- Sequential Path (for Windows debugging) ---
    cat("=========================================================\n")
    cat("Running", n_reps, "simulations SEQUENTIALLY for debugging...\n")
    results_list <- list()
    for (i in 1:n_reps) {
      results_list[[i]] <- run_one_simulation_svyglm(i, n_sample = sample_size, model_type = current_model)
    }
  }
  
  cat("Simulations for this scenario complete.\n")
  
  # --- Aggregate results ---
  successful_runs <- results_list[!sapply(results_list, inherits, "try-error")]
  if (length(successful_runs) > 0) {
    all_model_results[[current_model]] <- do.call(rbind, successful_runs)
    # Save intermediate files from successful runs
    for (i in seq_along(successful_runs)) {
      iteration_number <- if(RUN_IN_PARALLEL) which(!sapply(results_list, inherits, "try-error"))[i] else i
      file_name <- sprintf("sim_svyglm_boot_%04d.rds", iteration_number)
      saveRDS(successful_runs[[i]], file = file.path(results_dir, file_name))
    }
  }
}

# --- Combine and save final results ---
if(length(all_model_results) > 0) {
  successful_results <- do.call(rbind, all_model_results)
  rownames(successful_results) <- NULL
  results_filename <- "data/simulation_results_svyglm_boot.rds"
  saveRDS(successful_results, file = results_filename)
  cat("\n---------------------------------------------------------\n")
  cat("Aggregated results for the svyGLM simulation saved to:\n")
  cat("-", results_filename, "\n")
  cat("---------------------------------------------------------\n")
}
