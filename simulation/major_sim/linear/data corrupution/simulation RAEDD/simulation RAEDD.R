library("data.table")
library("rstudioapi")
library("parallel")
library("doParallel")
library("nleqslv")
library("pROC")





#settings
K_list<-c(60,80,100)
n_list<-c(5000,10000,15000)
K=60
n=5000
theta<-c(2,1)




#define functions
"%^%" <- function(x, n){with(eigen(x), vectors %*% (values^n * t(vectors)))}
huber_func<-function(vector,c){
  for(i in 1:length(vector)){
    if(vector[i]>c){vector[i]<-c}
    if(vector[i]<(-c)){vector[i]<--c}
  }
  return(vector)
}
robust_aggregation<-function(theta,theta_vector,n_vector,sigma,c=3){
  result<-0
  for(k in 1:dim(theta_vector)[1]){
    #print(sqrt(n_vector[k])*sigma%^%(-1/2)%*%(theta-theta_vector[k,]))
    result=result+sqrt(n_vector[k])*huber_func(sqrt(n_vector[k])*sigma%^%(-1/2)%*%(theta-theta_vector[k,]),c)
  }
  return(result)
}



#tau
constant=1.345
b=pnorm(constant,0,1)-pnorm(-constant,0,1)
sigma_psi_squared=((1-(constant*dnorm(constant,0,1)-(-constant)*dnorm(-constant,0,1))/(pnorm(constant,0,1)-pnorm(-constant,0,1)))*(pnorm(constant,0,1)-pnorm(-constant,0,1))+constant^2*(1-b))
tau=b^2/sigma_psi_squared
tau
1/tau




#parallel computing
totalCores = detectCores()
cluster <- makeCluster(totalCores-2) 
registerDoParallel(cluster)



combo<-expand.grid(K = K_list, n = n_list)


#data generation

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path("simulation data"))
location=dirname(rstudioapi::getActiveDocumentContext()$path)
foreach(index = 1:dim(combo)[1])%dopar%{
  library("data.table")
  library("rstudioapi")
  library("parallel")
  library("doParallel")
  library("nleqslv")
  library("pROC")
  "%^%" <- function(x, n){with(eigen(x), vectors %*% (values^n * t(vectors)))}
  huber_func<-function(vector,c){
    for(i in 1:length(vector)){
      if(vector[i]>c){vector[i]<-c}
      if(vector[i]<(-c)){vector[i]<--c}
    }
    return(vector)
  }
  robust_aggregation<-function(theta,theta_vector,n_vector,sigma,c=3){
    result<-0
    for(k in 1:dim(theta_vector)[1]){
      #print(sqrt(n_vector[k])*sigma%^%(-1/2)%*%(theta-theta_vector[k,]))
      result=result+sqrt(n_vector[k])*huber_func(sqrt(n_vector[k])*sigma%^%(-1/2)%*%(theta-theta_vector[k,]),c)
    }
    return(result)
  }
  
  K=combo[index,]$K
  n=combo[index,]$n
  folder=paste0("K=",K,"_n=",n)
  #dir.create(file.path(folder))
  start_row=1
  if(file.exists(paste0(location,"//K=",K,"_n=",n,"_AEDD",".csv"))){start_row=dim(read.csv(paste0(location,"//K=",K,"_n=",n,"_AEDD",".csv")))[1]}
  for(i in start_row:1000){
    set.seed(i)
    for(attack in c("X","Y")){#"Omniscient","Gaussian",
      for(k in 1:K){
        X<-rnorm(n*2)
        if(attack=="X"&(k>(K-floor(K^(1/4))))){x<- rchisq(n*2,30)}
        dim(X)<-c(n,2)
        e<-rnorm(n)
        Y<-X%*%theta+e
        if(attack=="Y"&(k>(K-floor(K^(1/4))))){
          Y<-Y-0.5
          Y<-Y*(rbinom(length(Y),1,0.8)-0.5)*2
          Y<-Y+0.5
        }
        simu_data<-data.frame(cbind(Y,X))
        fwrite(simu_data,paste0(folder,'//',k,".csv"))
        
      }
    #FSE
    ################################################################
    # full_sample<-NULL
    # for(k in 1:K){
    #   full_sample<-rbind(full_sample,read.csv(paste0(folder,'//',k,".csv")))
    # }
    # FSE<-glm(Y ~.-1,family=binomial(link='logit'),data=full_sample)
    # FSE_estimate<-FSE$coefficients
    # FSE_SE<-sqrt(diag(vcov(FSE)))
    # FSE_lower<- FSE_estimate-1.96*FSE_SE
    # FSE_upper<- FSE_estimate+1.96*FSE_SE
    # CI<-(FSE_lower<theta)&(FSE_upper>theta)
    # estimation<-data.frame(t(data.frame(as.vector( cbind(FSE_estimate,FSE_SE,CI)))))
    # fwrite(estimation,paste0(location,"//K=",K,"_n=",n,"_FSE",".csv"),append=TRUE,col.names=FALSE)
    ################################################################
    #RAEDD
    ################################################################
    theta_list<-NULL
    sigma_list<-NULL
    for(k in 1:K){
      subsample<-read.csv(paste0(folder,'//',k,".csv"))
      model_sub<-lm(X1 ~.-1,data=subsample)
      theta_list<-rbind(theta_list,model_sub$coefficients)
      sigma_list<-rbind(sigma_list,vcov(model_sub))
    }
    full_sigma=sigma_list[1:2,1:2]
    RAEDD_estimate<-nleqslv(c(2,1),robust_aggregation,theta_vector=theta_list,n_vector=rep(n,K),sigma=full_sigma*n,c=constant,
                            control=list(allowSingular=TRUE))$x
    RAEDD_SE<-sqrt(diag(full_sigma)/K/tau)
    RAEDD_lower<- RAEDD_estimate-1.96*RAEDD_SE
    RAEDD_upper<- RAEDD_estimate+1.96*RAEDD_SE
    CI<-(RAEDD_lower<theta)&(RAEDD_upper>theta)
    estimation<-data.frame(t(data.frame(as.vector( cbind(RAEDD_estimate,RAEDD_SE,CI)))))
    fwrite(estimation,paste0(location,"//K=",K,"_n=",n,"_RAEDD_",attack,".csv"),append=TRUE,col.names=FALSE)
    }
    ################################################################
  }
}








