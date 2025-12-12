# ==============================================================================
# Description: LASSO estimation methods
# ==============================================================================

soft_operator <- function(rho_j, lambda){
  
  sign(rho_j) * pmax(abs(rho_j) - lambda, 0)

}

# coordinate descent
cd_lasso <- function(x, y, lambda, beta_init = NULL, tol = 1e-6, maxiter = 1000){
  
  n <- nrow(x)
  p <- ncol(x)
  beta <- rep(0, p)
  
  # warm start
  if(is.null(beta_init)) {
    beta <- rep(0, p) 
  } else {
    beta <- beta_init 
  }
  
  r <- y - x %*% beta 
  z <- colSums(x^2)
  
  for(iter in 1:maxiter){
    beta_old <- beta
    
    for(j in 1:p){
      
      if(j == 1){
        lam_j <- 0
      } else {
        lam_j <- lambda
      }
      
      rho <- sum(x[,j] * r) + beta[j] * z[j]
      
      beta_new_j <- soft_operator(rho, lam_j)/z[j]
      
      r <- r - (beta_new_j - beta[j]) * x[, j]
      
      beta[j] <- beta_new_j
    }
    if (sum((beta - beta_old)^2) < tol) {
      break
    }
  }
  return(beta)
}

# proximal gradient descent
pgd_lasso <- function(x, y, lambda, eta, beta_init = NULL, tol = 1e-6, maxiter = 1000){
  n <- nrow(x)
  p <- ncol(x)
  
  # warm start
  if(is.null(beta_init)) {
    beta <- rep(0, p)
  } else {
    beta <- beta_init
  }

  for(iter in 1:maxiter){
    beta_old <- beta
    gradient <- t(x) %*% (x %*% beta_old - y)
    u <- beta_old - eta * gradient
    
    # Proximal Step 
    beta[1] <- soft_operator(u[1], 0) # intercept (no penalty)
    if(p > 1){
      beta[2:p] <- soft_operator(u[2:p], lambda * eta)
    }
    
    if (sum((beta - beta_old)^2) < tol) {
      break
    }
  }
  return(beta)
}


inner_prods <- function(x, y){
  p <- ncol(x)
  corr_vec <- numeric(p)
  for(j in 1:p){
    corr_vec[j] <- sum(x[,j] * y) 
  }
  return(abs(corr_vec))
}

# screening rules
screen_variables <- function(x, y, lambda){
  p  <- ncol(x)
  cor_vec <- inner_prods(x, y)
  lambda_max <- max(cor_vec)
  
  idx <- which(cor_vec >= (2*lambda - lambda_max), arr.ind = TRUE)
  
  return(idx)
}

screen_variables_sequential <- function(x, y, lambda, lambda_prev, beta_prev){
  p  <- ncol(x)
  res <- y - (x %*% beta_prev)
  cor_vec <- inner_prods(x, res)
  
  idx <- which(cor_vec >= (2*lambda - lambda_prev), arr.ind = TRUE)
  
  return(idx)
  
}
