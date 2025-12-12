# ==============================================================================
# Description: Generates synthetic datasets for Lasso simulation and saves to disk.
# ==============================================================================

# function to generate data
library(MASS)
generate_data <- function(n, p, sparsity, correlation_type = "independent", rho = 0.5, snr = 3) {
  
  
  k <- floor(p * sparsity)  
  beta <- rep(0, p)
  
  active_indices <- sample(1:p, k)
  beta[active_indices] <- runif(k, min = 0.5, max = 2) * sample(c(-1, 1), k, replace = TRUE)
  
  
  if (correlation_type == "independent") {
    Sigma <- diag(1, p)
  } else if (correlation_type == "block_diagonal") {
    
    block_size = 10
    Sigma <- outer(1:p, 1:p, function(i, j) {
      same_block <- ceiling(i / block_size) == ceiling(j / block_size)
      ifelse(same_block, rho^abs(i - j), 0)
    })
  } else {
    stop("Unknown correlation type")
  }
  
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)
  signal <- X %*% beta
  signal_variance <- var(as.vector(signal))
  noise_variance <- signal_variance / snr
  sigma_epsilon <- sqrt(noise_variance)
  y <- signal + rnorm(n, mean = 0, sd = sigma_epsilon)
  
  return(list(
    X = cbind(1, X),
    y = y,
    beta_true = beta,
    active_indices = sort(active_indices),
    sigma_epsilon = sigma_epsilon
  ))
}

# 1. Define Grid
param_grid <- expand.grid(
  n = c(100, 300, 500),
  p = c(10, 50, 200),
  sparsity = c(0.2, 0.4, 0.8),
  corr_type = c("independent", "block_diagonal"),
  snr = c(3, 1.0), 
  stringsAsFactors = FALSE
)


n_reps <- 50 
output_file <- "data/simulation_datasets.rds"
all_datasets <- list()
counter <- 1

set.seed(2025)

# 2. Generation Loop
# ------------------------------------------------------------------------------
total_scenarios <- nrow(param_grid)

for(i in 1:total_scenarios) {
  params <- param_grid[i, ]
  
  if(i %% 100 == 0 || i == 1) {
  cat(sprintf("\nGenerating Scenario %d/%d: n=%d, p=%d, s=%.2f, corr=%s, snr = %.2f\n", 
              i, total_scenarios, params$n, params$p, params$sparsity, params$corr_type, param_grid$snr))
  }
  for(rep in 1:n_reps) {

    data <- generate_data(
      n = params$n, 
      p = params$p, 
      sparsity = params$sparsity, 
      correlation_type = params$corr_type, 
      rho = 0.2, 
      snr = params$snr
    )

    X_no_int <- data$X[, -1] 
    lambda_max <- max(abs(t(X_no_int) %*% data$y))
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_max * 0.01), length.out = 50))

    all_datasets[[counter]] <- list(
      id = counter,
      scenario_id = i,
      replication = rep,
      params = params,
      data = data,       
      lambda_seq = lambda_seq
    )
    
    counter <- counter + 1
  }
}

# 3. Save
# ------------------------------------------------------------------------------
saveRDS(all_datasets, output_file)
