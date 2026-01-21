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


#NULL model
eta <- -1.428596425
eps <- 1e-15                         # to avoid log(0)
p_hat <- 1 / (1 + exp(-eta))
y <- df_cleaned$delayed
y <- y[!is.na(y)]
deviance <- (
  y * log(p_hat) + (1 - y) * log(1 - p_hat)
)


# Residual model
theta_hat <- c(
  "(Intercept)" = -1.234,
  "Year" = -0.018,
  "Duration" = 0.114,
  
  # Day of week
  "Tuesday"    = -0.046,
  "Wednesday"  = 0.050,
  "Thursday"   = 0.203,
  "Friday"     = 0.252,
  "Saturday"   = -0.159,
  "Sunday"     = -0.020,
  
  # Departure time (Morning is reference)
  "DepTAfternoon" = 0.573,
  "DepTEvening"   = 0.633,
  "DepTNight"     = -0.245,
  
  # Month (January is reference)
  "February"  = -0.054,
  "March"     = -0.125,
  "April"     = -0.332,
  "May"       = -0.335,
  "June"      = -0.033,
  "July"      = -0.094,
  "August"    = -0.141,
  "September" = -0.507,
  "October"   = -0.347,
  "November"  = -0.248,
  "December"  = 0.150
)

# X <- df_cleaned %>%
#   select(delayed, Year, Month, DayOfWeek, ActualElapsedTime, schedule_Dep)
df_no_na<- na.omit(df_cleaned)
X <- model.matrix( ~ Year+Month+DayOfWeek+ActualElapsedTime+schedule_Dep,
         data = (df_no_na))

eta <- X %*% theta_hat  
eps <- 1e-15                         # to avoid log(0)
p_hat <- 1 / (1 + exp(-eta))
p_hat <- p_hat[!is.na(p_hat)]
y <- df_no_na$delayed
#y <- y[!is.na(y)]
deviance <- -2 * sum(
  y * log(p_hat) + (1 - y) * log(1 - p_hat)
)
deviance <- as.numeric(deviance)
deviance <- deviance[!is.na(deviance)]
deviance <- -2*sum(deviance)
-2*sum(deviance[1:529059])


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








