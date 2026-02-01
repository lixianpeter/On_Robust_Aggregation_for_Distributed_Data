set.seed(42)

# Parameters
m <- 100
epsilon <- 0.05
delta <- 0.10
n_clean <- floor((1 - epsilon) * m)
n_contam <- m - n_clean

# Generate data
clean_data <- rnorm(n_clean, mean = 0, sd = 1)

# Place contamination exactly at 1 (your original code says rep(1, n_contam), not 50)
contam_data <- rep(1, n_contam)

# Combine data
X <- c(clean_data, contam_data)

# Compute distance matrix
D <- as.matrix(dist(X, diag = TRUE, upper = TRUE))

# Compute R_i
k <- ceiling((1 - delta) * m)
R <- apply(D, 1, function(d) sort(d)[k + 1])

# Find i* with smallest R
i_star <- which.min(R)
R_star <- R[i_star]

# Points within this radius
V_indices <- which(D[i_star, ] <= R_star)

# Data points in neighborhood
X_V <- X[V_indices]

# Huber loss function for location theta and data y
huber_loss <- function(theta, y, delta = 1.0) {
  r <- y - theta
  loss <- ifelse(abs(r) <= delta,
                 0.5 * r^2,
                 delta * (abs(r) - 0.5 * delta))
  sum(loss)
}

# Find Huber estimator on X_V with delta tuned (e.g., 1)
delta_huber <- 1.345
huber_est <- optim(par = 0, fn = huber_loss, y = X, delta = delta_huber, method = "BFGS")$par
huber_est
# Ordinary mean (full data)
mean_est <- mean(X)

# Report
cat(sprintf("Huber estimator:    %.3f\n", huber_est))
cat(sprintf("Ordinary mean:      %.3f\n", mean_est))
cat(sprintf("True center:        0.000\n"))
cat(sprintf("Number in V:        %d\n", length(V_indices)))
cat(sprintf("Center point X[i*]: %.3f\n", X[i_star]))
n_contam_in_V <- sum(X[V_indices] == 1)
cat(sprintf("Contaminated points in V: %d\n", n_contam_in_V))

# Plot
plot(X, rep(0, length(X)), pch=19, 
     col=ifelse(1:length(X) %in% V_indices, "red", 
                ifelse(X == 1, "darkgreen", "black")),
     xlab="X values", ylab="", yaxt='n', 
     main="Huber estimator on selected neighborhood")
points(X[i_star], 0, pch=19, col="blue", cex=2)
legend("topright", legend=c("Selected neighborhood", "Contamination cluster", "Center point"), 
       col=c("red", "darkgreen", "blue"), pch=19)
