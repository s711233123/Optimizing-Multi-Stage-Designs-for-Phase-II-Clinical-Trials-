setwd("D:/NTPU/PSO/kStageSimon-main/")
library(globpso)
source("util.R")
source("kStageP2A_Objective.R")

# Clinical Trial Requirements
cliRequirement <- list(
  p0 = 0.1,    # response rate in the null hypothesis
  p1 = 0.3,    # response rate in the alternative hypothesis
  alpha = 0.1, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

nStage <- 3 # Number of stages
n1Min <- 1  # Minimal sample size at the first stage
nrMin <- 1  # Minimal sample size at stages after the 1st one
nMaxRange <- c(10, 70) # Total sample size range
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1)) # Minimal sample sizes at each stage

# Constraints for the objective function
upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))

# Initialize storage for results
results <- list()
seeds <- 1:50

# Loop over random seeds
for (seed in seeds) {
  cat(sprintf("Running for seed: %d\n", seed))
  
  # --- Single PSO Run ---
  single_algSetting <- getPSOInfo(nSwarm = 256, maxIter = 400, psoType = "basic")
  single_result <- globpso(
    objFunc = kStageMinMaxObj, PSO_INFO = single_algSetting,
    lower = lower, upper = upper, seed = seed, verbose = FALSE,
    nMin = nMinEachInterim, cliRequirement = cliRequirement
  )
  
  # --- Two-Stage PSO (100 + 300) ---
  stage1_algSetting <- getPSOInfo(nSwarm = 256, maxIter = 100, psoType = "basic")
  stage1_result <- globpso(
    objFunc = kStageMinMaxObj, PSO_INFO = stage1_algSetting,
    lower = lower, upper = upper, seed = seed, verbose = FALSE,
    nMin = nMinEachInterim, cliRequirement = cliRequirement
  )
  
  upper2 <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  lower2 <- c(stage1_result$par[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  stage2_algSetting <- getPSOInfo(nSwarm = 256, maxIter = 300, psoType = "basic")
  stage2_result <- globpso(
    objFunc = kStageOptimObj, PSO_INFO = stage2_algSetting,
    lower = lower2, upper = upper2, seed = seed, verbose = FALSE,
    nMin = nMinEachInterim, cliRequirement = cliRequirement
  )
  
  # --- Two-Stage PSO (100 + 100) ---
  stage2_100_algSetting <- getPSOInfo(nSwarm = 256, maxIter = 100, psoType = "basic")
  stage2_100_result <- globpso(
    objFunc = kStageOptimObj, PSO_INFO = stage2_100_algSetting,
    lower = lower2, upper = upper2, seed = seed, verbose = FALSE,
    nMin = nMinEachInterim, cliRequirement = cliRequirement
  )
  
  # Save results for this seed
  results[[as.character(seed)]] <- list(
    single = list(val = single_result$val, cputime = single_result$cputime),
    two_stage_100_300 = list(val = stage2_result$val, cputime = stage1_result$cputime + stage2_result$cputime),
    two_stage_100_100 = list(val = stage2_100_result$val, cputime = stage1_result$cputime + stage2_100_result$cputime)
  )
}

# Summarize results
summary_df <- data.frame(
  Seed = seeds,
  SingleTime = sapply(results, function(x) x$single$cputime),
  TwoStageTime_100_300 = sapply(results, function(x) x$two_stage_100_300$cputime),
  TwoStageTime_100_100 = sapply(results, function(x) x$two_stage_100_100$cputime),
  Success_Single = sapply(results, function(x) x$single$val < 17.80),
  Success_TwoStage_100_300 = sapply(results, function(x) x$two_stage_100_300$val < 17.80),
  Success_TwoStage_100_100 = sapply(results, function(x) x$two_stage_100_100$val < 17.80)
)

# Calculate success rates and average computing times
success_rates <- colMeans(summary_df[, grep("Success_", names(summary_df))])
avg_times <- colMeans(summary_df[, grep("Time", names(summary_df))])

# Print summary
cat("Success Rates:\n")
print(success_rates)
cat("\nAverage Computing Times:\n")
print(avg_times)
