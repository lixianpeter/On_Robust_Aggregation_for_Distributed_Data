#########################################################---
# 1. Required packages ----
#########################################################---

library("data.table")





#########################################################---
# 2. Simulation settings ----
#########################################################---

#settings
K <- 60
n <- 5000
N <- n*K
R <- 100 #number of replications
theta <- c(2,1)

#choose one contamination setting:
#"Uncontaminated", "Omniscient", "Gaussian",
#"Bit-flip", or "Covariate"
attack <- "Uncontaminated"

#number of contaminated machines
number_attacked <- floor(K^(1/4))





#########################################################---
# 3. Functions used in the simulation ----
#########################################################---

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
  
  return(sigma_hat)
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
    "model=linear_weighted_average",
    "_K=",K,
    "_n=",n,
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
  
  theta_list <- NULL
  X_list <- list()
  Y_list <- list()
  
  for(k in 1:K){
    
    #generate the clean covariates
    X <- rnorm(n*2)
    dim(X) <- c(n,2)
    
    #generate the response using the clean covariates
    e <- rnorm(n)
    Y <- as.numeric(X%*%theta+e)
    
    #replace both covariates on the attacked machines
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
    
    theta_list <- rbind(
      theta_list,
      model_sub$coefficients
    )
    
    X_list[[k]] <- X
    Y_list[[k]] <- Y
  }
  
  
  #########################################################---
  ## 5.2 Apply the estimate contamination setting ----
  #########################################################---
  
  theta_list_attack <- theta_list
  
  if(attack=="Omniscient"){
    theta_list_attack[1:number_attacked,] <- -10^6
  }
  
  if(attack=="Gaussian"){
    
    theta_list_attack[1:number_attacked,] <- matrix(
      rnorm(
        number_attacked*length(theta),
        mean=0,
        sd=sqrt(200)
      ),
      nrow=number_attacked,
      ncol=length(theta)
    )
  }
  
  if(attack=="Bit-flip"){
    theta_list_attack[1:number_attacked,] <-
      theta_list_attack[1:number_attacked,]*(-1)
  }
  
  
  #########################################################---
  ## 5.3 Calculate the local covariance estimators ----
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
  ## 5.4 Calculate the weighted average estimator ----
  #########################################################---
  
  weighted_estimate <- colMeans(theta_list_attack)
  
  
  #########################################################---
  ## 5.5 Calculate the weighted covariance estimator ----
  #########################################################---
  
  weighted_sigma <- matrix(
    0,
    nrow=length(theta),
    ncol=length(theta)
  )
  
  for(k in 1:K){
    weighted_sigma <- weighted_sigma+(n/N)*sigma_list[[k]]
  }
  
  
  #########################################################---
  ## 5.6 Calculate standard errors and confidence intervals ----
  #########################################################---
  
  weighted_SE <- sqrt(diag(weighted_sigma)/N)
  
  weighted_lower <- weighted_estimate-1.96*weighted_SE
  weighted_upper <- weighted_estimate+1.96*weighted_SE
  
  CI <- (weighted_lower<theta)&(weighted_upper>theta)
  
  
  #########################################################---
  ## 5.7 Save the simulation result ----
  #########################################################---
  
  estimation <- data.frame(
    t(
      data.frame(
        as.vector(
          cbind(
            weighted_estimate,
            weighted_SE,
            CI
          )
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
}





#########################################################---
# 6. End of simulation script ----
#########################################################---

cat("\nSimulation completed.\n")
cat("Output file:",output_file,"\n")