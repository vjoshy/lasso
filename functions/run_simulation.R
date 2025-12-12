# ==============================================================================
# Description: Loads datasets and benchmarks Lasso algorithms (not parallel).
# ==============================================================================

source("lasso.R") 

# 1. Load Data
# ------------------------------------------------------------------------------
input_file <- "data/simulation_datasets.rds"
all_datasets <- readRDS(input_file)
results_list <- list()
counter <- 1

cat(sprintf("Loaded %d datasets. Starting benchmarking...\n", length(all_datasets)))

# 2. Methods
# ------------------------------------------------------------------------------
methods_to_test <- list(
  list(name="CD_Basic",   method="cd",  screen=FALSE),
  list(name="CD_Screen",  method="cd",  screen=TRUE),
  list(name="PGD_Basic",  method="pgd", screen=FALSE),
  list(name="PGD_Screen", method="pgd", screen=TRUE),
  list(name="Ground Truth", method="glmnet", screen=FALSE)
)

# 3. Main Loop
# ------------------------------------------------------------------------------
for(k in 1:length(all_datasets)) {

  sim_obj <- all_datasets[[k]]
  
  X <- sim_obj$data$X
  y <- sim_obj$data$y
  true_active <- sim_obj$data$active_indices
  lambda_seq <- sim_obj$lambda_seq
  p_info <- sim_obj$params

  if(k %% 10 == 0 || k == 1) {
    cat(sprintf("Processing Dataset %d/%d (n=%d, p=%d)...\n", 
                k, length(all_datasets), p_info$n, p_info$p))
  }

  for(m in methods_to_test) {
    start_time <- Sys.time()
    cv_res <- cv_lasso(X, y, lambda_seq, k_folds = 5, method = m$method, screen = m$screen)
    
    runtime <- as.numeric(Sys.time() - start_time, units = "secs")
    final_path <- solve_lasso_path(X, y, lambda_seq, method = m$method, screen = m$screen)
    best_idx <- which(lambda_seq == cv_res$best_lambda)
    beta_final <- final_path[, best_idx]

    tol <- 1e-4
    est_active <- which(abs(beta_final[-1]) > tol) 
    
    TP <- length(intersect(est_active, true_active))
    FP <- length(setdiff(est_active, true_active))
    FN <- length(setdiff(true_active, est_active))
    
    precision <- TP / max(1, (TP + FP))
    recall <- TP / max(1, (TP + FN))
    f1 <- 2 * (precision * recall) / max(1e-10, (precision + recall))

    results_list[[counter]] <- data.frame(
      Dataset_ID = sim_obj$id,
      Scenario_ID = sim_obj$scenario_id,
      Replication = sim_obj$replication,
      n = p_info$n,
      p = p_info$p,
      sparsity = p_info$sparsity,
      corr_type = p_info$corr_type,
      snr = p_info$snr,
      Method = m$method,
      Screen = m$screen,
      Runtime = runtime,
      Selected_Lambda = cv_res$best_lambda,
      MSE = cv_res$min_mse,
      Precision = precision,
      Recall = recall,
      F1_Score = f1,
      Num_Selected = length(est_active),
      Num_True = length(true_active)
    )
    counter <- counter + 1
  }
}

# 4. Save Results
# ------------------------------------------------------------------------------
final_results_df <- do.call(rbind, results_list)
saveRDS(final_results_df, "final_simulation_results.rds")
