#define functions
"%^%" <- function(x, n){with(eigen(x), vectors %*% (values^n * t(vectors)))}
huber_func<-function(vector,c){
  for(i in 1:length(vector)){
    if(vector[i]>c){vector[i]<-c}
    if(vector[i]<(-c)){vector[i]<--c}
  }
  return(vector)
}



simulate_logistic_estimates <- function(n_k, p, beta_true, K) {
  estimates <- matrix(NA, nrow = K, ncol = p)
  
  for (sim in 1:K) {
    # Simulate data
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    eta <- X %*% beta_true
    pi <- 1 / (1 + exp(-eta))
    y <- rbinom(n, size = 1, prob = pi)
    
    # Fit logistic regression
    data <- data.frame(y = y, X)
    fit <- glm(y ~ . -1, data = data, family = binomial())
    
    # Store estimated coefficients (skip intercept)
    estimates[sim, ] <- t(vcov(fit)%^%(-1/2)%*%(coef(fit)-beta_true))
    
  }
  
  #return(list(estimates = estimates, cov_matrices = vcov(fit)))
  return(estimates)
}




# Function as defined
compute_tau <- function(e_i, c) {
  n <- length(e_i)
  
  inside <- abs(e_i) <= c
  outside <- !inside
  
  b_hat <- sum(inside)
  numerator <- (b_hat)^2
  denom <- n * sum(inside * (e_i^2) + outside * (c^2))
  
  tau_hat <- numerator / denom
  return(tau_hat)
}













set.seed(123)
n_k <- 1000
p <- 3
beta_true <- c(0.5, -1, 2)
K <- 50
n <- 100  # sample size
estimates <- simulate_logistic_estimates(n_k, p, beta_true, K)




# Parameters
set.seed(123)
n_sim <- 100
n <- 100  # sample size
c_values <- seq(0, 3, length.out = 10)
avg_tau <- numeric(length(c_values))
# Run simulations with just normal
for (j in seq_along(c_values)) {
  c <- c_values[j]
  tau_sims <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    estimates <- simulate_logistic_estimates(n_k, p, beta_true, K)
    #covar <- fit$cov_matrices
    #estimates <- fit$estimates
    #estimates[28:30,1:3] <- 10
    e <- estimates
    tau_values <- compute_tau(as.vector(e),c)
    tau_sims[i] <- tau_values
  }
  
  avg_tau[j] <- mean(tau_sims)
}


# Plot
plot(c_values, avg_tau, type = "l", lwd = 2,
     xlab = "c", ylab = "Average tau_hat",
     main = "Average tau_hat vs c over 100 simulations under M-estimator without contamination")
grid()









# Parameters
set.seed(123)
n_sim <- 100
n <- 100  # sample size
c_values <- seq(0, 3, length.out = 10)
avg_tau <- numeric(length(c_values))
# Run simulations with just normal
for (j in seq_along(c_values)) {
  c <- c_values[j]
  tau_sims <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    estimates <- simulate_logistic_estimates(n_k, p, beta_true, K)
    #covar <- fit$cov_matrices
    #estimates <- fit$estimates
    #estimates[28:30,1:3] <- 10
    e <- estimates
    e[48:50,1:3] <- 1000
    tau_values <- compute_tau(as.vector(e),c)
    tau_sims[i] <- tau_values
  }
  
  # Append the new rows to the existing CSV
  write.table(t(data.frame(c(c,tau_sims))), "my_data.csv", 
              append = TRUE, 
              sep = ",", 
              col.names = T, 
              row.names = T)
  avg_tau[j] <- mean(tau_sims)
}


# Plot
plot(c_values, avg_tau, type = "l", lwd = 2,
     xlab = "c", ylab = "Average tau_hat",
     main = "Average tau_hat vs c over 100 simulations under M-estimator with contamination of 1000")
grid()











