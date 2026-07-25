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
R <- 100 #number of replications
theta <- c(2,1)

#choose one data contamination setting:
#"Uncontaminated" or "Covariate"
attack <- "Covariate"

#number of contaminated machines
number_attacked <- floor(K^(1/4))

#huber tuning parameter
huber_c <- 1.345

#asymptotic relative efficiency
b_c <- 2*pnorm(huber_c)-1
sigma_c_squared <- (1-b_c)*huber_c^2-2*dnorm(huber_c)*huber_c+b_c
tau_c <- b_c^2/sigma_c_squared





#########################################################---
# 3. Functions used in the simulation ----
#########################################################---

#define matrix power function
"%^%" <- function(x,n){
  with(eigen(x),vectors%*%(values^n*t(vectors)))
}


#define elementwise huber function
huber_func <- function(vector,c){
  
  for(i in 1:length(vector)){
    if(vector[i]>c){vector[i] <- c}
    if(vector[i]<(-c)){vector[i] <- -c}
  }
  
  return(vector)
}


#define the robust aggregation from K machines
robust_aggregation <- function(theta,theta_vector,n,sigma,c=3){
  
  result <- rep(0,length(theta))
  
  for(k in 1:dim(theta_vector)[1]){
    
    standardised_difference <- sqrt(n)*sigma%^%(-1/2)%*%
      (theta-theta_vector[k,])
    
    result <- result+sqrt(n)*
      huber_func(as.numeric(standardised_difference),c)
  }
  
  return(as.numeric(result))
}


#define the local asymptotic covariance estimator
#this is for linear regression only
local_sigma_estimator <- function(theta,X,Y){
  
  local_n <- dim(X)[1]
  
  residual <- as.numeric(Y-X%*%theta)
  score <- X*residual
  mean_score <- colMeans(score)
  
  mean_score_matrix <- matrix(
    mean_score,
    nrow=local_n,
    ncol=dim(X)[2],
    byrow=TRUE
  )
  
  centred_score <- score-mean_score_matrix
  
  U_hat <- t(X)%*%X/local_n
  V_hat <- t(centred_score)%*%centred_score/local_n
  
  sigma_hat <- solve(U_hat)%*%V_hat%*%solve(U_hat)
  
  #cap each matrix entry between -10^3 and 10^3
  #otherwise numerical optimization could handle too large value
  sigma_hat[sigma_hat>10^3] <- 10^3
  sigma_hat[sigma_hat<(-10^3)] <- -10^3
  
  return(sigma_hat)
}


#define the half vectorization of a symmetric matrix
vector_half_matrix <- function(sigma){
  
  sigma_vector <- sigma[lower.tri(sigma,diag=TRUE)]
  
  return(sigma_vector)
}


#define a symmetric matrix recovered from its lower triangular vector
symmetric_matrix <- function(sigma_vector,matrix_dimension){
  
  sigma <- matrix(
    0,
    nrow=matrix_dimension,
    ncol=matrix_dimension
  )
  
  sigma[lower.tri(sigma,diag=TRUE)] <- sigma_vector
  
  sigma <- sigma+t(sigma)
  diag(sigma) <- diag(sigma)/2
  
  return(sigma)
}


#define the spatial median L2 loss function
spatial_median_loss <- function(spatial_median,sigma_vector_list,n){
  
  loss <- 0
  
  #go through each machine's sigma_k
  for(k in 1:dim(sigma_vector_list)[1]){
    
    difference <- sigma_vector_list[k,]-spatial_median
    loss <- loss+sqrt(n)*sqrt(sum(difference^2))
  }
  
  return(loss)
}


#define the spatial median of the covariance matrices
spatial_median_sigma <- function(sigma_list,n){
  
  number_matrices <- length(sigma_list)
  matrix_dimension <- dim(sigma_list[[1]])[1]
  vector_dimension <- matrix_dimension*(matrix_dimension+1)/2
  
  sigma_vector_list <- matrix(
    0,
    nrow=number_matrices,
    ncol=vector_dimension
  )
  
  for(k in 1:number_matrices){
    sigma_vector_list[k,] <- vector_half_matrix(sigma_list[[k]])
  }
  
  initial_value <- colMeans(sigma_vector_list)
  
  spatial_median_vector <- optim(
    initial_value,
    spatial_median_loss,
    sigma_vector_list=sigma_vector_list,
    n=n,
    method="Nelder-Mead"
  )$par
  
  spatial_median_matrix <- symmetric_matrix(
    spatial_median_vector,
    matrix_dimension
  )
  
  return(spatial_median_matrix)
}


#define the contamination detection method
detect_contaminated_machines <- function(
    theta_vector,
    robust_estimate,
    n,
    sigma,
    alpha=0.05
){
  
  number_machines <- dim(theta_vector)[1]
  parameter_dimension <- length(robust_estimate)
  
  chi_squared_cutoff <- qchisq(
    1-alpha,
    df=parameter_dimension
  )
  
  sigma_inverse <- solve(sigma)
  detected <- rep(FALSE,number_machines)
  
  for(k in 1:number_machines){
    
    difference <- as.numeric(
      theta_vector[k,]-robust_estimate
    )
    
    distance_squared <- n*as.numeric(
      t(difference)%*%sigma_inverse%*%difference
    )
    
    detected[k] <- distance_squared>chi_squared_cutoff
  }
  
  detected_machines <- which(detected)
  
  return(detected_machines)
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
    "model=linear_robust",
    "_K=",K,
    "_n=",n,
    "_c=",huber_c,
    "_attack=",attack,
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
  ## 5.1 Generate the local datasets and estimators ----
  #########################################################---
  
  theta_list_attack <- NULL
  X_list <- list()
  Y_list <- list()
  
  for(k in 1:K){
    
    #generate the clean covariates
    X <- rnorm(n*2)
    dim(X) <- c(n,2)
    
    #generate the response using the clean covariates
    e <- rnorm(n)
    Y <- as.numeric(X%*%theta+e)
    
    #replace both covariates on each attacked machine
    #the response generated above remains unchanged
    if(attack=="Covariate" && k<=number_attacked){
      X[,1] <- rchisq(n,df=30)
      X[,2] <- rchisq(n,df=30)
    }
    
    simu_data <- data.frame(
      Y=Y,
      X1=X[,1],
      X2=X[,2]
    )
    
    model_sub <- lm(
      Y~.-1,
      data=simu_data
    )
    
    theta_list_attack <- rbind(
      theta_list_attack,
      model_sub$coefficients
    )
    
    X_list[[k]] <- X
    Y_list[[k]] <- Y
  }
  
  
  #########################################################---
  ## 5.2 Calculate the local covariance estimators ----
  #########################################################---
  
  sigma_list <- list()
  
  for(k in 1:K){
    
    sigma_list[[k]] <- local_sigma_estimator(
      theta=theta_list_attack[k,],
      X=X_list[[k]],
      Y=Y_list[[k]]
    )
  }
  
  
  #########################################################---
  ## 5.3 Calculate the spatial median covariance matrix ----
  #########################################################---
  
  full_sigma <- spatial_median_sigma(
    sigma_list=sigma_list,
    n=n
  )
  
  
  #########################################################---
  ## 5.4 Calculate the robust aggregated estimator ----
  #########################################################---
  
  #use the true parameter as the initial value
  theta_start <- theta
  
  robust_estimate <- nleqslv(
    theta_start,
    robust_aggregation,
    theta_vector=theta_list_attack,
    n=n,
    sigma=full_sigma,
    c=huber_c,
    control=list(allowSingular=TRUE)
  )$x
  
  
  #########################################################---
  ## 5.5 Detect contamination and calculate hit rate ----
  #########################################################---
  
  detected_machines <- detect_contaminated_machines(
    theta_vector=theta_list_attack,
    robust_estimate=robust_estimate,
    n=n,
    sigma=full_sigma,
    alpha=0.05
  )
  
  if(attack=="Uncontaminated"){
    
    hit_rate <- NA_real_
    
  }else{
    
    contaminated_machines <- 1:number_attacked
    
    correctly_detected <- intersect(
      detected_machines,
      contaminated_machines
    )
    
    hit_rate <- length(correctly_detected)/number_attacked
  }
  
  
  #########################################################---
  ## 5.6 Calculate standard errors and confidence intervals ----
  #########################################################---
  
  robust_SE <- sqrt(diag(full_sigma)/(tau_c*N))
  
  robust_lower <- robust_estimate-1.96*robust_SE
  robust_upper <- robust_estimate+1.96*robust_SE
  
  CI <- (robust_lower<theta)&(robust_upper>theta)
  
  
  #########################################################---
  ## 5.7 Save the simulation result ----
  #########################################################---
  
  estimation <- data.frame(
    t(
      data.frame(
        c(
          as.vector(
            cbind(
              robust_estimate,
              robust_SE,
              CI
            )
          ),
          hit_rate
        )
      )
    )
  )
  
  fwrite(
    estimation,
    output_file,
    append=TRUE,
    col.names=FALSE
  )
  
  
  #########################################################---
  ## 5.8 Print the simulation progress ----
  #########################################################---
  
  cat("Replication",i,"out of",R,"\n")
  flush.console()
}





#########################################################---
# 6. End of simulation script ----
#########################################################---

cat("\nSimulation completed.\n")
cat("Output file:",output_file,"\n")