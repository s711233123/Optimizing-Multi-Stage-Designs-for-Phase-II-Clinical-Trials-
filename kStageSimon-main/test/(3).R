getwd()
setwd("D:/NTPU/PSO/kStageSimon-main/")
library(globpso)
# Import two local R files "util.R" and "kStageP2A_Objective.R"
source("util.R")
source("kStageP2A_Objective.R")

# Set Requirements of the Clinical Trial
cliRequirement <- list(
  p0 = 0.55,    # response rate in the null hypothesis
  p1 = 0.75,    # response rate in the alternative hypothesis
  alpha = 0.1, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

# Set the required number of stages
nStage <- 3

# Set Constraints of the Clinical Trial Design
n1Min <- 1 # minimal sample size at the first stage
nrMin <- 1  # minimal sample size at stages after the 1st one
nMaxRange <- c(15, 70) # the range of the total sample size
# Generate the constraint vector of minimal sample sizes at each stage
# This is the input of the objective function
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1)) # (do not change)

# Set the PSO configuration
algSetting <- getPSOInfo(
  nSwarm = 300,     # swarm size 
  maxIter = 300,    # number of iterations
  psoType = "basic" # PSO type (one can use "basic" or "quantum" or "cso")
)

# Set Random seed for reproducibility
#pso_seed <- 1

# Set the lower and upper bounds for PSO search
#  particle = (nMax, nPolarized, rProportionEachInterim) 
#  with sizes (1, nStage - 1, nStage)
upper <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0*pi, nStage - 1), rep(0, nStage))

### -------------------------------------------------------- ###
### Find Optimal Design 
### -------------------------------------------------------- ###
# Run PSO for optimal design
optimRes  <- globpso(objFunc = kStageOptimObj, PSO_INFO = algSetting, 
                     lower = lower, upper = upper, 
                     seed = NULL, verbose = TRUE,
                     nMin = nMinEachInterim, cliRequirement = cliRequirement)

# View the optimal design search results
optimRes$val     # Objective function value
optimRes$cputime # computing time

# Transform the PSO outcome into the readable optimal design
optimDesign <- kStageFreqCrit(
  nPolarized = optimRes$par[2:nStage],  
  rProportion = optimRes$par[(nStage + 1):length(optimRes$par)], 
  nMax = optimRes$par[1], nMin = nMinEachInterim, cliRequirement)

# The resulting optimal design
optimDesign$nseq # Sample sizes at each stage (n_1, ..., n_K)
optimDesign$rseq # Stopping cutoff sizes at each stage (r_1, ..., r_K)

# Properties of the resulting optimal design
optimDesign$t1e # Type I error
optimDesign$t2e # Type II error
optimDesign$en  # Expected sample size under null hypothesis
optimDesign$pet_seq # The probabilities of early termination of the trial at each stage

