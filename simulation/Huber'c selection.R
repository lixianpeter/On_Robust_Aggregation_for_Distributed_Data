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

# Parameters
set.seed(123)
n_sim <- 1000
n <- 100  # sample size
c_values <- seq(0, 3, length.out = 50)
avg_tau <- numeric(length(c_values))

# Run simulations with just normal
for (j in seq_along(c_values)) {
  c <- c_values[j]
  tau_sims <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    e_i <- rnorm(n)
    tau_sims[i] <- compute_tau(e_i, c)
  }
  
  avg_tau[j] <- mean(tau_sims)
}

# Plot
plot(c_values, avg_tau, type = "l", lwd = 2,
     xlab = "c", ylab = "Average tau_hat",
     main = "Average tau_hat vs c over 1000 simulations under normal without contamination")
grid()




# Run simulations with contamination
epsilon <- 0.1
contam_idx <- sample(1:n, size = floor(epsilon * n))
for (j in seq_along(c_values)) {
  c <- c_values[j]
  tau_sims <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    e_i <- rnorm(n)
    e_i[contam_idx] <- rnorm(length(contam_idx), mean = 0, sd = 10)
    tau_sims[i] <- compute_tau(e_i, c)
  }
  
  avg_tau[j] <- mean(tau_sims)
}

# Plot
plot(c_values, avg_tau, type = "l", lwd = 2,
     xlab = "c", ylab = "Average tau_hat",
     main = "Average tau_hat vs c over 1000 simulations 90% from N(0,1) and 10% from N(0,10)")
grid()


# Run simulations with contamination
epsilon <- 0.1
contam_idx <- sample(1:n, size = floor(epsilon * n))
for (j in seq_along(c_values)) {
  c <- c_values[j]
  tau_sims <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    e_i <- rnorm(n)
    e_i[contam_idx] <- 1000
    tau_sims[i] <- compute_tau(e_i, c)
  }
  
  avg_tau[j] <- mean(tau_sims)
}

# Plot
plot(c_values, avg_tau, type = "l", lwd = 2,
     xlab = "c", ylab = "Average tau_hat",
     main = "Average tau_hat vs c over 1000 simulations 90% from N(0,1) and 10% from values of 1000 (Extreme value)")
grid()

# Run simulations with contamination
epsilon <- 0.4
contam_idx <- sample(1:n, size = floor(epsilon * n))
for (j in seq_along(c_values)) {
  c <- c_values[j]
  tau_sims <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    e_i <- rnorm(n)
    e_i[contam_idx] <- 1000
    tau_sims[i] <- compute_tau(e_i, c)
  }
  
  avg_tau[j] <- mean(tau_sims)
}

# Plot
plot(c_values, avg_tau, type = "l", lwd = 2,
     xlab = "c", ylab = "Average tau_hat",
     main = "Average tau_hat vs c over 1000 simulations 60% from N(0,1) and 40% from values of 1000 (Extreme value)")
grid()






