library("data.table")
library("rstudioapi")
library("parallel")
library("doParallel")
library("nleqslv")
library("pROC")
library("Matrix")
library("caret")

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







#data generation

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

location=dirname(rstudioapi::getActiveDocumentContext()$path)
#########################################################
theta_list <- read.csv("server_carrier.csv",header=TRUE)
cov_list <- read.csv("server_carrier_cov.csv",header=TRUE)





full_sigma<-as.matrix((matrix(as.numeric(read.csv("server_carrier.csv",header=FALSE)[1,c(-1,-2)]),
                              nrow=sqrt(dim(read.csv("servers_cov.csv",header=FALSE))[2]-2))))
theta_list<-as.matrix(na.omit(read.csv("server_carrier.csv"))[c(-1,-2)])
n_list<-as.vector(na.omit(read.csv("servers.csv"))[2])$n
RAEDD_estimate<-nleqslv(as.vector(as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])[1,]),
                        robust_aggregation,theta_vector=theta_list,n_vector=n_list,sigma=full_sigma*n_list[1],c=constant,
                        control=list(allowSingular=TRUE))$x
RAEDD_SE<-sqrt(diag(full_sigma)/tau)
AEDD_estimate<-nleqslv(as.vector(as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])[1,]),
                       robust_aggregation,theta_vector=theta_list,n_vector=n_list,sigma=full_sigma*n_list[1],c=3000,
                        control=list(allowSingular=TRUE))$x

saveRDS(RAEDD_estimate,"RAEDD_estimate")
saveRDS(AEDD_estimate,"AEDD_estimate")
saveRDS(RAEDD_SE,"RAEDD_SE")



# for(attack in c("Omniscient","Gaussian","Bit-flip","None")){
#   theta_list<-as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])
#   K=dim(theta_list)[1]
#   if(attack=="Omniscient"){theta_list[1:floor(K^(1/4)),]<--10^10}
#   if(attack=="Gaussian"){theta_list[1:floor(K^(1/4)),]<-rnorm(floor(K^(1/4))*dim(theta_list)[2],0,200^2)}
#   if(attack=="Bit-flip"){theta_list[1:floor(K^(1/4)),]<-theta_list[1:floor(K^(1/4)),]*(-1)}
#   RAEDD_estimate<-nleqslv(as.vector(as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])[1,]),robust_aggregation,theta_vector=theta_list,n_vector=n_list,sigma=full_sigma*n_list[1],c=constant,
#                           control=list(allowSingular=TRUE))$x
#   RAEDD_SE<-sqrt(diag(full_sigma)/tau)
#   estimation<-data.frame(t(data.frame(as.vector( cbind(RAEDD_estimate,RAEDD_SE)))))
#   estimation<-cbind(attack,constant,estimation)
#   fwrite(estimation,paste0("servers_result.csv"),append=TRUE,col.names=FALSE)
# }





# 
# f = list.files("servers")[1]
# subsample<-read.csv(paste0("servers//",f))
# ArrDelay<-subsample$ArrDelay
# subsample$ArrDelay[ArrDelay>15]<-1
# subsample$ArrDelay[ArrDelay<=15]<-0
# subsample$DayOfWeek<-as.character(subsample$DayOfWeek)
# subsample$Month<-as.character(subsample$Month)
# DepTime<-subsample$DepTime
# subsample$DepTime[600<=DepTime&DepTime<1200]<-"Morning"
# subsample$DepTime[1200<=DepTime&DepTime<1800]<-"Afternoon"
# subsample$DepTime[1800<=DepTime&DepTime<2200]<-"Evening"
# subsample$DepTime[subsample$DepTime!="Morning"&subsample$DepTime!="Afternoon"&subsample$DepTime!="Evening"]<-"Night"
# CRSDepTime<-subsample$CRSDepTime
# subsample$CRSDepTime[600<=CRSDepTime&DepTime<1200]<-"Morning"
# subsample$CRSDepTime[1200<=CRSDepTime&CRSDepTime<1800]<-"Afternoon"
# subsample$CRSDepTime[1800<=CRSDepTime&CRSDepTime<2200]<-"Evening"
# subsample$CRSDepTime[subsample$CRSDepTime!="Morning"&subsample$CRSDepTime!="Afternoon"&subsample$CRSDepTime!="Evening"]<-"Night"
# CRSArrTime<-subsample$CRSArrTime
# subsample$CRSArrTime[600<=CRSArrTime&CRSArrTime<1200]<-"Morning"
# subsample$CRSArrTime[1200<=CRSArrTime&CRSArrTime<1800]<-"Afternoon"
# subsample$CRSArrTime[1800<=CRSArrTime&CRSArrTime<2200]<-"Evening"
# subsample$CRSArrTime[subsample$CRSArrTime!="Morning"&subsample$CRSArrTime!="Afternoon"&subsample$CRSArrTime!="Evening"]<-"Night"
# subsample[names(select_if(subsample, is.numeric))[-1]]<-data.frame(scale(subsample[names(select_if(subsample, is.numeric))[-1]]))
# model_sub<-glm(formula = ArrDelay ~ Year + Month + DayOfWeek + DepTime + CRSDepTime + CRSArrTime  + Distance ,family = binomial(link = "logit"), data = subsample)
# 
# model_sub$coefficients<-as.vector(as.vector(as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])[1,]))
# model_sub$coefficients<--10^10
# # 1. Logistic regression
# logit_P = predict(model_sub, newdata =  subsample,type = 'response' )
# dmy <- dummyVars(" ~ .", data = subsample[ c("Year","Month","DayOfWeek","DepTime","CRSDepTime","CRSArrTime", "Distance")])
# trsf <- data.frame(predict(dmy, newdata = subsample[ c("Year","Month","DayOfWeek","DepTime","CRSDepTime","CRSArrTime", "Distance")]))
# X<-cbind(1,as.matrix(trsf))
# n <- nrow(X)
# w<-logit_P*(1-logit_P)
# v<- Diagonal(n, w)
# var_b<-solve(t(X)%*%v%*%X)





# for(attack in c("Omniscient","Gaussian","Bit-flip","None")){
#   theta_list<-as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])
#   K=dim(theta_list)[1]
#   if(attack=="Omniscient"){theta_list[1:floor(K^(1/4)),]<--10^10}
#   if(attack=="Gaussian"){theta_list[1:floor(K^(1/4)),]<-rnorm(floor(K^(1/4))*dim(theta_list)[2],0,200^2)}
#   if(attack=="Bit-flip"){theta_list[1:floor(K^(1/4)),]<-theta_list[1:floor(K^(1/4)),]*(-1)}
#   RAEDD_estimate<-nleqslv(as.vector(as.matrix(na.omit(read.csv("servers.csv"))[c(-1,-2)])[1,]),robust_aggregation,theta_vector=theta_list,n_vector=n_list,sigma=full_sigma*n_list[1],c=3000,
#                           control=list(allowSingular=TRUE))$x
#   RAEDD_SE<-sqrt(diag(cov(theta_list)))
#   estimation<-data.frame(t(data.frame(as.vector( cbind(RAEDD_estimate,RAEDD_SE)))))
#   estimation<-cbind(attack,3000,estimation)
#   fwrite(estimation,paste0("servers_result.csv"),append=TRUE,col.names=FALSE)
# }
################################################################
################################################################
################################################################


