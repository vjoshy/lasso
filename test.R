set.seed(123)
n <- 100
p <- 5
X <- scale(matrix(rnorm(n * p), n, p))
X <- cbind(1, X)
beta_true <- c(1, 3, -1.5, 0, 0, 2)
y <- X %*% beta_true + rnorm(n, 0, 1)


beta_cd <- cd_lasso(X, y, lambda = 5)
beta_pgd <- pgd_lasso(X, y, lambda = 5, eta = 0.01)



# trry ----------------------------------------------------------------------------------------

set.seed(123)
n <- 100
p <- 10 

X_preds <- scale(matrix(rnorm(n * p), n, p)) 
X <- cbind(1, X_preds) 

# True beta: Intercept + 5 real signals
beta_true <- c(1, 3, -1.5, 2, 0, 0, rep(0, p-5)) 
y <- X %*% beta_true + rnorm(n)

corrs <- abs(t(X[,-1]) %*% y)
lambda_max <- max(corrs)
lambda_seq <- exp(seq(log(lambda_max), log(lambda_max * 0.05), length.out = 100))

# Store results here
results_beta <- matrix(0, nrow = ncol(X), ncol = 100)
kept_counts <- numeric(100) 

for(i in 1:length(lambda_seq)){
  l <- lambda_seq[i]
  
  if(i == 1){
    keep_indices <- screen_variables(X[,-1, drop = FALSE], y, lambda = l) +1
  } else {
    keep_indices <- screen_variables_sequential(X, y, l, lambda_seq[i-1], results_beta[,i-1]) 
  }
  
  active_set <- sort(unique(c(1, keep_indices)))
  subset_beta <- pgd_lasso(X[, active_set, drop = FALSE], y, lambda = l, eta = 0.01)
  results_beta[active_set, i] <- subset_beta
  
  kept_counts[i] <- length(active_set)
}
