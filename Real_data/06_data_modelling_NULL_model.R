library(data.table)
library(readr)
library(ggplot2)
library(fst)
library(tidyr)
library(dplyr)
library("data.table")
library("rstudioapi")
library("parallel")
library("doParallel")
library("nleqslv")
library("pROC")
library("pracma")
library("ks")
#change working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
location=dirname(rstudioapi::getActiveDocumentContext()$path)

#read data
df_cleaned<-read_fst("us_airline_cleaned.fst")
df_cleaned <- df_cleaned %>% mutate(Month = as.factor(Month),
                                    DayOfWeek = as.factor(DayOfWeek),
                                    schedule_Dep = as.factor(schedule_Dep)) %>% 
  mutate(schedule_Dep = relevel(schedule_Dep, ref = "Morning"))

server_result<-NULL
theta_list<-NULL
sigma_list<-NULL


for(carrier in unique(df_cleaned$UniqueCarrier)){
  model_sub <- glm(delayed ~ 1,
                   data = (df_cleaned %>%filter(UniqueCarrier==carrier)), family = binomial)
  # Extract components
  # chi_sq <- model_sub$null.deviance - model_sub$deviance
  # df <- model_sub$df.null - model_sub$df.residual
  # p_val <- pchisq(chi_sq, df, lower.tail = FALSE)
  # theta <- c(carrier,model_sub$null.deviance, model_sub$deviance, chi_sq, df, p_val)
  # names(theta) <- c("carrier","Null Deviance", "Residual Deviance", "test stat", "df",  "p-value")
  # write.table(data.frame(t(theta)),  "server_carrier_deviance.csv", sep = ",", col.names = !file.exists("server_carrier.csv"), row.names = FALSE, append = T)
  theta <- c(carrier,dim(df_cleaned %>%filter(UniqueCarrier==carrier))[1],model_sub$coefficients)
  # names(theta)[c(1,2)] <- c("carrier", "n_k")
  # sigma <- c(carrier, dim(df_cleaned %>%filter(UniqueCarrier==carrier))[1]*vech(vcov(model_sub)))
  # names(sigma)[c(1)] <- c("carrier")
  write.table(data.frame(t(theta)),  "server_carrier_NULL_model.csv", sep = ",", col.names = !file.exists("server_carrier_NULL_model.csv"), row.names = FALSE, append = T)
  #write.table(data.frame(t(sigma)),  "server_carrier_cov.csv", sep = ",", col.names = !file.exists("server_carrier_cov.csv"), row.names = FALSE, append = T)
}








