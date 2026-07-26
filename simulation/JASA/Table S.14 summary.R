#########################################################---
# 1. Required package ----
#########################################################---

library("data.table")





#########################################################---
# 2. Summary settings ----
#########################################################---

#simulation settings
K <- 60
n <- 5000

#folder containing the simulation results
input_folder <- "output result"

#individual replication results
input_file <- file.path(
  input_folder,
  paste0(
    "table_S14_simulation",
    "_K=",K,
    "_n=",n,
    ".csv"
  )
)

#Table S.14 summary results
output_file <- file.path(
  input_folder,
  paste0(
    "table_S14_summary",
    "_K=",K,
    "_n=",n,
    ".csv"
  )
)






#########################################################---
# 3. Function used in the summary ----
#########################################################---

#define the relative efficiency
calculate_efficiency <- function(
    benchmark_estimates,
    method_estimates
){
  
  efficiency <- 100*
    var(benchmark_estimates)/
    var(method_estimates)
  
  return(efficiency)
}





#########################################################---
# 4. Read the individual replication results ----
#########################################################---

simulation_result <- fread(
  input_file,
  header=FALSE
)

names(simulation_result) <- c(
  "weighted_average_clean",
  "median_clean",
  "quantile_clean",
  "huber_c1_clean",
  "huber_c2_clean",
  "median_contaminated",
  "quantile_contaminated",
  "huber_c1_contaminated",
  "huber_c2_contaminated"
)





#########################################################---
# 5. Calculate RE under no contamination ----
#########################################################---

RE_median <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$median_clean
)

RE_quantile <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$quantile_clean
)

RE_huber_c1 <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$huber_c1_clean
)

RE_huber_c2 <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$huber_c2_clean
)





#########################################################---
# 6. Calculate RE dagger under Omniscient Contamination ----
#########################################################---

RE_dagger_median <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$median_contaminated
)

RE_dagger_quantile <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$quantile_contaminated
)

RE_dagger_huber_c1 <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$huber_c1_contaminated
)

RE_dagger_huber_c2 <- calculate_efficiency(
  benchmark_estimates=simulation_result$weighted_average_clean,
  method_estimates=simulation_result$huber_c2_contaminated
)





#########################################################---
# 7. Create the Table S.14 summary ----
#########################################################---

summary_result <- data.table(
  
  setting=c(
    "RE under contamination-free setting",
    "RE dagger under Omniscient Contamination"
  ),
  
  K=c(
    K,
    K
  ),
  
  n=c(
    n,
    n
  ),
  
  median=c(
    RE_median,
    RE_dagger_median
  ),
  
  averaged_quantile=c(
    RE_quantile,
    RE_dagger_quantile
  ),
  
  huber_c_0.9818=c(
    RE_huber_c1,
    RE_dagger_huber_c1
  ),
  
  huber_c_1.345=c(
    RE_huber_c2,
    RE_dagger_huber_c2
  )
)





#########################################################---
# 8. Round and save the summary ----
#########################################################---

result_columns <- c(
  "median",
  "averaged_quantile",
  "huber_c_0.9818",
  "huber_c_1.345"
)

summary_result[
  ,
  (result_columns) := lapply(
    .SD,
    round,
    digits=1
  ),
  .SDcols=result_columns
]

fwrite(
  summary_result,
  output_file
)





#########################################################---
# 9. End of summary script ----
#########################################################---

cat("\nSummary completed.\n")
cat("Output file:",output_file,"\n")