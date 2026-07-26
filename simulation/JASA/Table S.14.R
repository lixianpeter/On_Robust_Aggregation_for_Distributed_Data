#########################################################---
# 1. Required packages ----
#########################################################---

library("data.table")
library("nleqslv")





#########################################################---
# 2. Simulation settings ----
#########################################################---

#settings
K <- 60
n <- 5000
N <- n*K
R <- 1000 # number of replication
theta <- 0

#number of contaminated machines
number_attacked <- floor(K^(1/4))

#huber tuning parameters
huber_c1 <- 0.9818
huber_c2 <- 1.345

#quantiles used by the averaged-quantile aggregation
quantile_probabilities <- c(
  0.15,
  0.25,
  0.50,
  0.75,
  0.85
)





#########################################################---
# 3. Functions used in the simulation ----
#########################################################---

#define the elementwise huber function
huber_func <- function(vector,c){
  
  for(j in 1:length(vector)){
    
    if(vector[j]>c){
      vector[j] <- c
    }
    
    if(vector[j]<(-c)){
      vector[j] <- -c
    }
  }
  
  return(vector)
}


#define the one-dimensional huber aggregation equation
robust_aggregation <- function(
    theta,
    theta_vector,
    n,
    sigma,
    c
){
  
  standardised_difference <- sqrt(n)*
    (theta-theta_vector)/sqrt(sigma)
  
  result <- sum(
    sqrt(n)*
      huber_func(
        standardised_difference,
        c
      )
  )
  
  return(result)
}


#define the averaged-quantile aggregation
averaged_quantile_aggregation <- function(
    theta_vector,
    quantile_probabilities
){
  
  selected_quantiles <- quantile(
    theta_vector,
    probs=quantile_probabilities,
    names=FALSE,
    type=7
  )
  
  averaged_quantile <- mean(
    selected_quantiles
  )
  
  return(averaged_quantile)
}





#########################################################---
# 4. Output and simulation setup ----
## 4.1 Set the output folder ----
#########################################################---

output_folder <- "output result"

if(!dir.exists(output_folder)){
  dir.create(output_folder)
}


#########################################################---
## 4.2 Set the output file ----
#########################################################---

output_file <- file.path(
  output_folder,
  paste0(
    "table_S14_simulation",
    "_K=",K,
    "_n=",n,
    ".csv"
  )
)

if(file.exists(output_file)){
  file.remove(output_file)
}





#########################################################---
# 5. Main simulation loop ----
#########################################################---

for(i in 1:R){
  
  set.seed(i)
  
  
  #########################################################---
  ## 5.1 Generate the local sample means and variances ----
  #########################################################---
  
  local_estimates <- rep(
    0,
    K
  )
  
  local_sigma_estimates <- rep(
    0,
    K
  )
  
  for(k in 1:K){
    
    Y <- rnorm(
      n,
      mean=theta,
      sd=1
    )
    
    local_estimates[k] <- mean(
      Y
    )
    
    local_sigma_estimates[k] <- mean(
      (Y-local_estimates[k])^2
    )
  }
  
  
  #########################################################---
  ## 5.2 Calculate the spatial median variance estimator ----
  #########################################################---
  
  #in one dimension, the spatial median is the ordinary median
  full_sigma <- median(
    local_sigma_estimates
  )
  
  
  #########################################################---
  ## 5.3 Create the contaminated local estimates ----
  #########################################################---
  
  contaminated_estimates <- local_estimates
  
  contaminated_estimates[
    1:number_attacked
  ] <- -10^6
  
  
  #########################################################---
  ## 5.4 Calculate the clean weighted average benchmark ----
  #########################################################---
  
  weighted_average_clean <- mean(
    local_estimates
  )
  
  
  #########################################################---
  ## 5.5 Calculate the median aggregations ----
  #########################################################---
  
  median_clean <- median(
    local_estimates
  )
  
  median_contaminated <- median(
    contaminated_estimates
  )
  
  
  #########################################################---
  ## 5.6 Calculate the averaged-quantile aggregations ----
  #########################################################---
  
  quantile_clean <- averaged_quantile_aggregation(
    theta_vector=local_estimates,
    quantile_probabilities=quantile_probabilities
  )
  
  quantile_contaminated <- averaged_quantile_aggregation(
    theta_vector=contaminated_estimates,
    quantile_probabilities=quantile_probabilities
  )
  
  
  #########################################################---
  ## 5.7 Calculate the Huber aggregations with c=0.9818 ----
  #########################################################---
  
  huber_c1_clean <- nleqslv(
    theta,
    robust_aggregation,
    theta_vector=local_estimates,
    n=n,
    sigma=full_sigma,
    c=huber_c1,
    control=list(
      allowSingular=TRUE
    )
  )$x
  
  huber_c1_contaminated <- nleqslv(
    theta,
    robust_aggregation,
    theta_vector=contaminated_estimates,
    n=n,
    sigma=full_sigma,
    c=huber_c1,
    control=list(
      allowSingular=TRUE
    )
  )$x
  
  
  #########################################################---
  ## 5.8 Calculate the Huber aggregations with c=1.345 ----
  #########################################################---
  
  huber_c2_clean <- nleqslv(
    theta,
    robust_aggregation,
    theta_vector=local_estimates,
    n=n,
    sigma=full_sigma,
    c=huber_c2,
    control=list(
      allowSingular=TRUE
    )
  )$x
  
  huber_c2_contaminated <- nleqslv(
    theta,
    robust_aggregation,
    theta_vector=contaminated_estimates,
    n=n,
    sigma=full_sigma,
    c=huber_c2,
    control=list(
      allowSingular=TRUE
    )
  )$x
  
  
  #########################################################---
  ## 5.9 Save the simulation result ----
  #########################################################---
  
  estimation <- data.frame(
    weighted_average_clean,
    median_clean,
    quantile_clean,
    huber_c1_clean=as.numeric(huber_c1_clean),
    huber_c2_clean=as.numeric(huber_c2_clean),
    median_contaminated,
    quantile_contaminated,
    huber_c1_contaminated=as.numeric(huber_c1_contaminated),
    huber_c2_contaminated=as.numeric(huber_c2_contaminated)
  )
  
  fwrite(
    estimation,
    output_file,
    append=TRUE,
    col.names=FALSE
  )
  
  
  #########################################################---
  ## 5.10 Print the simulation progress ----
  #########################################################---
  
  cat(
    "Replication",
    i,
    "out of",
    R,
    "\n"
  )
}





#########################################################---
# 6. End of simulation script ----
#########################################################---

cat("\nSimulation completed.\n")
cat("Output file:",output_file,"\n")