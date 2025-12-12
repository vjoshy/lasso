# ==============================================================================
# Description: Loads datasets and benchmarks Lasso algorithms in parallel.
# ==============================================================================

library(foreach)
library(doParallel)
library(glmnet) 

source("lasso.R") 

input_file <- "data/simulation_datasets.rds"

all_datasets <- readRDS(input_file)
cat(sprintf("Loaded %d datasets.\n", length(all_datasets)))

# 2. Setup Parallel Backend
num_cores <- parallel::detectCores() - 1 
registerDoParallel(cores = num_cores)

cat(sprintf("Starting parallel simulation on %d cores...\n", num_cores))

# 3. Methods
# ------------------------------------------------------------------------------
methods_to_test <- list(
  list(name="CD_Basic",     method="cd",     screen=FALSE),
  list(name="CD_Screen",    method="cd",     screen=TRUE),
  list(name="PGD_Basic",    method="pgd",    screen=FALSE),
  list(name="PGD_Screen",   method="pgd",    screen=TRUE),
  list(name="Ground Truth", method="glmnet", screen=FALSE)
)

# 4. Main Parallel Loop
# ------------------------------------------------------------------------------

start <- Sys.time()
final_results_df <- foreach(k = 1:length(all_datasets), 
                            .combine = rbind, 
                            .packages = c("glmnet")) %dopar% {
                              sim_obj <- all_datasets[[k]]

                              X <- sim_obj$data$X
                              y <- sim_obj$data$y
                              true_active <- sim_obj$data$active_indices
                              lambda_seq <- sim_obj$lambda_seq
                              p_info <- sim_obj$params
                              
                              dataset_results <- list()
                              local_counter <- 1
                              
                              for(m in methods_to_test) {
                                
                                start_time <- Sys.time()
                                cv_res <- cv_lasso(X, y, lambda_seq, k_folds = 5, 
                                                   method = m$method, screen = m$screen)
                                
                                runtime <- as.numeric(Sys.time() - start_time, units = "secs")
                                final_path <- solve_lasso_path(X, y, lambda_seq, 
                                                               method = m$method, screen = m$screen)
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
                                
                                dataset_results[[local_counter]] <- data.frame(
                                  Dataset_ID = sim_obj$id,
                                  Scenario_ID = sim_obj$scenario_id,
                                  Replication = sim_obj$replication,
                                  n = p_info$n,
                                  p = p_info$p,
                                  sparsity = p_info$sparsity,
                                  corr_type = p_info$corr_type,
                                  snr = p_info$snr,
                                  Method_Name = m$name,      
                                  Method_Type = m$method,
                                  Screening = m$screen,
                                  Runtime = runtime,
                                  Selected_Lambda = cv_res$best_lambda,
                                  MSE = cv_res$min_mse,
                                  Precision = precision,
                                  Recall = recall,
                                  F1_Score = f1,
                                  Num_Selected = length(est_active),
                                  Num_True = length(true_active)
                                )
                                local_counter <- local_counter + 1
                              }

                              do.call(rbind, dataset_results)
                            }

# 5. Save Final Results
# ------------------------------------------------------------------------------
stopImplicitCluster()

saveRDS(final_results_df, "final_simulation_results.rds")
end <- Sys.time()
print(end - start)