# Load required packages and files
library(globpso)
source("util.R")
source("kStageP2A_Objective.R")

# Set Requirements of the Clinical Trial
cliRequirement <- list(
  p0 = 0.05,    # response rate in the null hypothesis
  p1 = 0.25,    # response rate in the alternative hypothesis
  alpha = 0.05, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

# Set the required number of stages
nStage <- 2

# Set Constraints of the Clinical Trial Design
n1Min <- 1 # minimal sample size at the first stage
nrMin <- 1 # minimal sample size at stages after the 1st one
nMaxRange <- c(10, 70) # the range of the total sample size
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1))

# Set the PSO configuration
algSetting <- getPSOInfo(
  nSwarm = 300,     # swarm size 
  maxIter = 200,    # number of iterations
  psoType = "basic"
)

# Set the lower and upper bounds for PSO search
upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))

# Initialize vectors to store results from 20 runs
num_runs <- 20
results <- vector("list", num_runs)

# Perform the PSO process 20 times
for (i in 1:num_runs) {
  # Run PSO for minimax design
  minMaxRes <- globpso(
    objFunc = kStageMinMaxObj, PSO_INFO = algSetting, 
    lower = lower, upper = upper, 
    seed = NULL, verbose = FALSE,
    nMin = nMinEachInterim, cliRequirement = cliRequirement
  )
  
  # Transform the PSO outcome into a readable minimax design
  minMaxDesign <- kStageFreqCrit(
    nPolarized = minMaxRes$par[2:nStage], 
    rProportion = minMaxRes$par[(nStage + 1):length(minMaxRes$par)], 
    nMax = minMaxRes$par[1], nMin = nMinEachInterim, cliRequirement
  )
  
  # Store results in the list
  results[[i]] <- list(
    nseq = minMaxDesign$nseq,
    rseq = minMaxDesign$rseq,
    t1e = minMaxDesign$t1e,
    t2e = minMaxDesign$t2e,
    en = minMaxDesign$en,
    pet_seq = minMaxDesign$pet_seq,
    obj_val = minMaxRes$val,
    cpu_time = minMaxRes$cputime
  )
}

# Print out the results of each run
for (i in 1:num_runs) {
  cat("Run", i, "results:\n")
  cat("Sample sizes (nseq):", results[[i]]$nseq, "\n")
  cat("Stopping cutoffs (rseq):", results[[i]]$rseq, "\n")
  cat("Type I error (t1e):", results[[i]]$t1e, "\n")
  cat("Type II error (t2e):", results[[i]]$t2e, "\n")
  cat("Expected sample size (en):", results[[i]]$en, "\n")
  cat("Early termination probabilities (pet_seq):", results[[i]]$pet_seq, "\n")
  cat("Objective function value:", results[[i]]$obj_val, "\n")
  cat("Computing time:", results[[i]]$cpu_time, "seconds\n\n")
}
