# ==============================================================================
# Description: Cross-validation and regularization path functions
# ==============================================================================

library(glmnet)

solve_lasso_path <- function(X, y, lambda_seq, method = "cd", screen = FALSE, eta = 0.0001) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (method == "glmnet") {

    X_features <- X[, -1] 

    fit <- glmnet(X_features, y, lambda = lambda_seq/n, alpha = 1, 
                  intercept = TRUE, standardize = FALSE)
    
    return(as.matrix(coef(fit)))
  }

  beta_mat <- matrix(0, nrow = p, ncol = length(lambda_seq))
  beta_prev <- rep(0, p) 
  
  for(i in 1:length(lambda_seq)) {
    l <- lambda_seq[i]
    

    if (screen) {
      if(i == 1) {
        keep_idx <- screen_variables(X[,-1], y, lambda = l) + 1
      } else {
        keep_idx <- screen_variables_sequential(X, y, l, lambda_seq[i-1], beta_prev)
      }
      active_set <- sort(unique(c(1, keep_idx))) 
    } else {
      active_set <- 1:p 
    }
    
    X_sub <- X[, active_set, drop = FALSE]
    if(i == 1) {
      beta_start <- NULL 
    } else {
      beta_start <- beta_prev[active_set]
    }
    
    if (method == "cd") {
      beta_sub <- cd_lasso(X_sub, y, lambda = l)
    } else if (method == "pgd") {
      beta_sub <- pgd_lasso(X_sub, y, lambda = l, eta = eta)
    }

    beta_full <- rep(0, p)
    beta_full[active_set] <- beta_sub
    
    beta_mat[, i] <- beta_full
    beta_prev <- beta_full 
  }
  
  return(beta_mat)
}

cv_lasso <- function(X, y, lambda_seq, k_folds = 5, method = "cd", screen = TRUE, eta = 0.0001) {
  n <- nrow(X)

  set.seed(123) 
  folds <- sample(rep(1:k_folds, length.out = n))
  
  mse_matrix <- matrix(NA, nrow = k_folds, ncol = length(lambda_seq))
  
  for (k in 1:k_folds) {

    test_idx <- which(folds == k)
    X_train <- X[-test_idx, ]; y_train <- y[-test_idx]
    X_test <- X[test_idx, ];   y_test <- y[test_idx]
    beta_path <- solve_lasso_path(X_train, y_train, lambda_seq, method, screen, eta)
    
    preds <- X_test %*% beta_path
    residuals <- sweep(preds, 1, y_test, "-")
    mse_matrix[k, ] <- colMeans(residuals^2)
  }

  mean_mse <- colMeans(mse_matrix)
  best_lambda_idx <- which.min(mean_mse)
  best_lambda <- lambda_seq[best_lambda_idx]
  
  return(list(
    best_lambda = best_lambda,
    min_mse = min(mean_mse),
    mean_mse_path = mean_mse
  ))
}