library(data.table)
library(readr)
library(ggplot2)
library(fst)
library(tidyr)
library(dplyr)
df<-read_fst("us_airline.fst")
df_cleaned <- df%>%
  mutate(
    delayed = ifelse(ArrDelay>15,1,0),
    schedule_Dep = case_when(
      CRSDepTime >= 600 & CRSDepTime < 1200 ~ "Morning",
      CRSDepTime >= 1200 & CRSDepTime < 1800 ~ "Afternoon",
      CRSDepTime >= 1800 & CRSDepTime < 2400 ~ "Evening",
      .default = "Night"
    ),
    schedule_Arr = case_when(
      CRSArrTime >= 600 & CRSArrTime < 1200 ~ "Morning",
      CRSArrTime >= 1200 & CRSArrTime < 1800 ~ "Afternoon",
      CRSArrTime >= 1800 & CRSArrTime < 2400 ~ "Evening",
      .default = "Night"
    ),
    #if the value greater than 0 then flag is 1
    security_delay_flag = ifelse(SecurityDelay>0, 1, 0),
    weather_delay_flag = ifelse(WeatherDelay>0, 1, 0),
    carrier_delay_flag = ifelse(CarrierDelay>0, 1, 0),
    security_delay_flag = ifelse(is.na(security_delay_flag), 0, SecurityDelay),
    weather_delay_flag = ifelse(is.na(weather_delay_flag ), 0, WeatherDelay),
    carrier_delay_flag = ifelse(is.na(carrier_delay_flag), 0, CarrierDelay)
  )
#Sanity check 
# min(df_cleaned %>% 
#       select(CRSDepTime,schedule_Dep) %>%
#       filter(schedule_Dep=="Evening") %>% 
#   select(CRSDepTime) 
# )
# max(df_cleaned %>% 
#       select(CRSDepTime,schedule_Dep) %>%
#       filter(schedule_Dep=="Evening") %>% 
#       select(CRSDepTime) 
# )
# 
# df_cleaned %>% 
#     select(CRSDepTime,CRSArrTime, DepTime,ArrTime, ArrDelay,DepDelay,CarrierDelay, WeatherDelay, NASDelay, SecurityDelay,delayed)  %>% 
#    filter(SecurityDelay>0)
write_fst(df_cleaned, "us_airline_cleaned.fst")






