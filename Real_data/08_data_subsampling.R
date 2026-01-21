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



N = nrow(df_cleaned)
# For 1% subsampling
set.seed(123)
df_sub_1 <- df_cleaned[sample(N, size = as.integer(0.01*N), replace = F), ]

model_sub_1 <- glm(delayed ~ Year+Month+DayOfWeek+ActualElapsedTime+schedule_Dep,
                 data = df_sub_1, family = binomial)
options(scipen = 999)
summary(model_sub_1)

# For 10% subsampling
set.seed(123)
df_sub_10 <- df_cleaned[sample(N, size = as.integer(0.1*N), replace = F), ]

model_sub_10 <- glm(delayed ~ Year+Month+DayOfWeek+ActualElapsedTime+schedule_Dep,
                   data = df_sub_10, family = binomial)
summary(model_sub_10)





