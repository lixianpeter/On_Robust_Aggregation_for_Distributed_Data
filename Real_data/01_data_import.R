library(data.table)
library(readr)
library(ggplot2)
library(fst)


df<-fread("combined_csv.csv")
write_fst(df, "us_airline.fst") 
