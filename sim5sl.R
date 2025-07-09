# Step 0: Install and load the necessary libraries
# install.packages(c("SuperLearner", "earth", "survey", "rsimsum", "parallel", "gam", "glmnet", "surveyCV"))
library(SuperLearner)
library(earth)
library(survey)
library(rsimsum)
library(parallel)
library(gam)
library(glmnet)
library(surveyCV)

# Step 1: Define the data simulation function with a model_type argument
create.data <- function(n_target_sample = 1000, te_log_odds = 1.5, seed = NULL, model_type = "complex"){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
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
  } else if (model_type == "simple") {
    L1 <- C1_s; L2 <- C2_s; L3 <- C3_s; L4 <- C4_s
  } else {
    stop("model_type must be either 'complex' or 'simple'")
  }
  
  prob_Y1_pop <- get.prob_Y(1, C1, C2, C3, C4); prob_Y0_pop <- get.prob_Y(0, C1, C2, C3, C4)
  true_ate_rd <- mean(prob_Y1_pop) - mean(prob_Y0_pop)
  obs <- data.frame(L1, L2, L3, L4, A, Y, strata = strata_s, cluster = cluster_s, weight = weights_s)
  return(list(observed = obs, ate = true_ate_rd))
}


# Step 2: Create a function to run one full simulation iteration
run_one_simulation <- function(i, n_sample, model_type, learners) {
  library(SuperLearner); library(earth); library(survey); library(gam); library(glmnet); library(surveyCV)
  
  sim_data_list <- create.data(n_target_sample = n_sample, seed = i, model_type = model_type)
  observed_data <- sim_data_list$observed
  true_ate <- sim_data_list$ate
  W_confounders <- observed_data[, c("L1", "L2", "L3", "L4")]
  
  lib_size <- length(learners) 
  
  # --- Method 1: Fully-Aware Analysis ---
  sl_q_w <- SuperLearner(Y = observed_data$Y, X = observed_data[, c("A", "L1", "L2", "L3", "L4")], family = binomial(), SL.library = learners, obsWeights = observed_data$weight)
  data_a1 <- data.frame(A = 1, W_confounders); data_a0 <- data.frame(A = 0, W_confounders)
  Q_1_W_w <- predict(sl_q_w, newdata = data_a1)$pred; Q_0_W_w <- predict(sl_q_w, newdata = data_a0)$pred; Q_A_W_w <- predict(sl_q_w, newdata = observed_data[, c("A", "L1", "L2", "L3", "L4")])$pred
  sl_g_w <- SuperLearner(Y = observed_data$A, X = W_confounders, family = binomial(), SL.library = learners, obsWeights = observed_data$weight)
  g1_W_w <- predict(sl_g_w, newdata = W_confounders)$pred; g1_W_w <- pmax(pmin(g1_W_w, 0.999), 0.001)
  H_A_W_w <- (observed_data$A / g1_W_w) - ((1 - observed_data$A) / (1 - g1_W_w))
  fluc_model_w <- glm(Y ~ -1 + H_A_W_w, data = observed_data, family = binomial(), offset = qlogis(Q_A_W_w), weights = weight)
  epsilon_w <- coef(fluc_model_w)
  Q_1_W_star_w <- plogis(qlogis(Q_1_W_w) + epsilon_w * (1 / g1_W_w)); Q_0_W_star_w <- plogis(qlogis(Q_0_W_w) + epsilon_w * (-1 / (1 - g1_W_w)))
  psi_w <- weighted.mean(Q_1_W_star_w, w = observed_data$weight) - weighted.mean(Q_0_W_star_w, w = observed_data$weight)
  Q_A_W_star_w <- ifelse(observed_data$A == 1, Q_1_W_star_w, Q_0_W_star_w)
  eif_w <- H_A_W_w * (observed_data$Y - Q_A_W_star_w) + (Q_1_W_star_w - Q_0_W_star_w) - psi_w
  svy_design <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = observed_data)
  svy_design <- update(svy_design, eif = eif_w)
  se_w <- as.numeric(SE(svymean(~eif, svy_design)))
  res1 <- data.frame(lib_size = lib_size, model = model_type, method = "Fully-Aware", true = true_ate, b = psi_w, se = se_w)
  
  # --- Method 2: Partially-Aware ---
  se_p <- sqrt(var(eif_w) / n_sample)
  res2 <- data.frame(lib_size = lib_size, model = model_type, method = "Partially-Aware", true = true_ate, b = psi_w, se = se_p)
  
  # --- Method 3: Non-Aware ---
  sl_q_uw <- SuperLearner(Y = observed_data$Y, X = observed_data[, c("A", "L1", "L2", "L3", "L4")], family = binomial(), SL.library = learners)
  Q_1_W_uw <- predict(sl_q_uw, newdata = data_a1)$pred; Q_0_W_uw <- predict(sl_q_uw, newdata = data_a0)$pred; Q_A_W_uw <- predict(sl_q_uw, newdata = observed_data[, c("A", "L1", "L2", "L3", "L4")])$pred
  sl_g_uw <- SuperLearner(Y = observed_data$A, X = W_confounders, family = binomial(), SL.library = learners)
  g1_W_uw <- predict(sl_g_uw, newdata = W_confounders)$pred; g1_W_uw <- pmax(pmin(g1_W_uw, 0.999), 0.001)
  H_A_W_uw <- (observed_data$A / g1_W_uw) - ((1 - observed_data$A) / (1 - g1_W_uw))
  fluc_model_uw <- glm(Y ~ -1 + H_A_W_uw, data = observed_data, family = binomial(), offset = qlogis(Q_A_W_uw))
  epsilon_uw <- coef(fluc_model_uw)
  Q_1_W_star_uw <- plogis(qlogis(Q_1_W_uw) + epsilon_uw * (1 / g1_W_uw)); Q_0_W_star_uw <- plogis(qlogis(Q_0_W_uw) + epsilon_uw * (-1 / (1 - g1_W_uw)))
  psi_uw <- mean(Q_1_W_star_uw) - mean(Q_0_W_star_uw)
  Q_A_W_star_uw <- ifelse(observed_data$A == 1, Q_1_W_star_uw, Q_0_W_star_uw)
  eif_uw <- H_A_W_uw * (observed_data$Y - Q_A_W_star_uw) + (Q_1_W_star_uw - Q_0_W_star_uw) - psi_uw
  se_uw <- sqrt(var(eif_uw) / n_sample)
  res3 <- data.frame(lib_size = lib_size, model = model_type, method = "Non-Aware", true = true_ate, b = psi_uw, se = se_uw)
  
  # --- Method 4: Fully-Aware with Survey-Aware Cross-Validation ---
  n_folds <- 5
  fold_IDs <- folds.svy(observed_data, nfolds = n_folds, clusterID = "cluster")
  validation_rows_list <- split(seq_len(nrow(observed_data)), fold_IDs)
  sl_control <- list(V = n_folds, validRows = validation_rows_list)
  sl_q_cv <- SuperLearner(Y = observed_data$Y, X = observed_data[, c("A", "L1", "L2", "L3", "L4")], family = binomial(), SL.library = learners, obsWeights = observed_data$weight, cvControl = sl_control)
  Q_1_W_cv <- predict(sl_q_cv, newdata = data_a1)$pred; Q_0_W_cv <- predict(sl_q_cv, newdata = data_a0)$pred; Q_A_W_cv <- predict(sl_q_cv, newdata = observed_data[, c("A", "L1", "L2", "L3", "L4")])$pred
  sl_g_cv <- SuperLearner(Y = observed_data$A, X = W_confounders, family = binomial(), SL.library = learners, obsWeights = observed_data$weight, cvControl = sl_control)
  g1_W_cv <- predict(sl_g_cv, newdata = W_confounders)$pred; g1_W_cv <- pmax(pmin(g1_W_cv, 0.999), 0.001)
  H_A_W_cv <- (observed_data$A / g1_W_cv) - ((1 - observed_data$A) / (1 - g1_W_cv))
  fluc_model_cv <- glm(Y ~ -1 + H_A_W_cv, data = observed_data, family = binomial(), offset = qlogis(Q_A_W_cv), weights = weight)
  epsilon_cv <- coef(fluc_model_cv)
  Q_1_W_star_cv <- plogis(qlogis(Q_1_W_cv) + epsilon_cv * (1 / g1_W_cv)); Q_0_W_star_cv <- plogis(qlogis(Q_0_W_cv) + epsilon_cv * (-1 / (1 - g1_W_cv)))
  psi_cv <- weighted.mean(Q_1_W_star_cv, w = observed_data$weight) - weighted.mean(Q_0_W_star_cv, w = observed_data$weight)
  Q_A_W_star_cv <- ifelse(observed_data$A == 1, Q_1_W_star_cv, Q_0_W_star_cv)
  eif_cv <- H_A_W_cv * (observed_data$Y - Q_A_W_star_cv) + (Q_1_W_star_cv - Q_0_W_star_cv) - psi_cv
  svy_design_cv <- svydesign(ids = ~cluster, strata = ~strata, weights = ~weight, data = observed_data)
  svy_design_cv <- update(svy_design_cv, eif = eif_cv)
  se_cv <- as.numeric(SE(svymean(~eif, svy_design_cv)))
  res4 <- data.frame(lib_size = lib_size, model = model_type, method = "Fully-Aware-CV", true = true_ate, b = psi_cv, se = se_cv)
  
  # --- Diagnostics ---
  diagnostics <- list(lib_size = lib_size, model_type = model_type,
                      fully_aware = list(q_coef = sl_q_w$coef, g_coef = sl_g_w$coef, epsilon = epsilon_w, g_summary = summary(g1_W_w)),
                      non_aware = list(q_coef = sl_q_uw$coef, g_coef = sl_g_uw$coef, epsilon = epsilon_uw, g_summary = summary(g1_W_uw)),
                      fully_aware_cv = list(q_coef = sl_q_cv$coef, g_coef = sl_g_cv$coef, epsilon = epsilon_cv, g_summary = summary(g1_W_cv)))
  
  return(list(results = rbind(res1, res2, res3, res4), diagnostics = diagnostics))
}


# Step 3: Define the simulation parameters
n_reps <- 1000 
sample_size <- 500
model_choices <- c("complex", "simple")

# You can easily change this single vector to run a new simulation.
learner_library <- c("SL.glm", "SL.gam", "SL.earth")

# Get the size of the library for naming files and folders
lib_size <- length(learner_library)


all_model_results <- list()
all_model_diagnostics <- list()

for (current_model in model_choices) {
  cat("=========================================================\n")
  cat("Running: Library Size =", lib_size, "| Model Type =", toupper(current_model), "\n")
  
  results_dir <- file.path("data", "intermediate_results", current_model, paste0("lib_size_", lib_size))
  if (!dir.exists(results_dir)) {dir.create(results_dir, recursive = TRUE)}
  
  no_cores <- detectCores() - 1 
  cl <- makeCluster(no_cores)
  
  clusterExport(cl, varlist=c("create.data", "run_one_simulation", "sample_size", "current_model", "learner_library", "results_dir"))
  
  cat("=========================================================\n")
  cat("Running", n_reps, "simulations in parallel on", no_cores, "cores...\n")
  
  parLapply(cl, 1:n_reps, function(i) {
    result_obj <- try(run_one_simulation(i, n_sample = sample_size, model_type = current_model, learners = learner_library))
    if (!inherits(result_obj, "try-error")) {
      file_name <- sprintf("sim_%04d.rds", i)
      saveRDS(result_obj, file = file.path(results_dir, file_name))
    }
  })
  
  stopCluster(cl)
  cat("Simulations for this scenario complete.\n")
  
  # Aggregate results for this specific model type and library size
  result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(result_files) > 0) {
    loaded_objects <- lapply(result_files, readRDS)
    all_model_results[[current_model]] <- do.call(rbind, lapply(loaded_objects, `[[`, "results"))
    all_model_diagnostics[[current_model]] <- lapply(loaded_objects, `[[`, "diagnostics")
  }
}

if(length(all_model_results) > 0) {
  successful_results <- do.call(rbind, all_model_results)
  successful_diagnostics <- do.call(c, all_model_diagnostics)
  rownames(successful_results) <- NULL
  
  results_filename <- paste0("data/simulation_results_", lib_size, ".rds")
  diagnostics_filename <- paste0("data/simulation_diagnostics_", lib_size, ".rds")
  
  saveRDS(successful_results, file = results_filename)
  saveRDS(successful_diagnostics, file = diagnostics_filename)
  
  cat("\n---------------------------------------------------------\n")
  cat("Aggregated results for library size", lib_size, "saved to:\n")
  cat("-", results_filename, "\n")
  cat("-", diagnostics_filename, "\n")
  cat("---------------------------------------------------------\n")
}

