### ======================================================== ###
### Example Code
### ------------
### Using PSO to generate K-stage clinical trial design
### under Simon's frequentist framework
### ======================================================== ###
setwd("D:/NTPU/PSO/kStageSimon-main/")
# Install and load the "globpso" package by following the 
#  instruction in https://github.com/PingYangChen/globpso
library(globpso)
# Import two local R files "util.R" and "kStageP2A_Objective.R"
source("util.R")
source("kStageP2A_Objective.R")

# Set Requirements of the Clinical Trial
cliRequirement <- list(
  p0 = 0.1,    # response rate in the null hypothesis
  p1 = 0.3,    # response rate in the alternative hypothesis
  alpha = 0.1, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

# Set the required number of stages
nStage <- 3

# Set Constraints of the Clinical Trial Design
n1Min <- 1 # minimal sample size at the first stage
nrMin <- 1  # minimal sample size at stages after the 1st one
nMaxRange <- c(10, 70) # the range of the total sample size
# Generate the constraint vector of minimal sample sizes at each stage
# This is the input of the objective function
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1)) # (do not change)



# 固定總樣本數並繼續尋找其他參數
upper <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0*pi, nStage - 1), rep(0, nStage))


# 設定PSO配置（初步搜索）
initial_algSetting <- getPSOInfo(
  nSwarm = 256,      # 初步設置的粒子數
  maxIter = 100,     # 初步設置的迭代次數
  psoType = "basic" # PSO類型
)

pso_seed <- 2

# 運行PSO以尋找minimax設計（初步搜索）
initial_minMaxRes <- globpso(
  objFunc = kStageMinMaxObj, PSO_INFO = initial_algSetting, 
  lower = lower, upper = upper, 
  seed = pso_seed, verbose = TRUE,
  nMin = nMinEachInterim, cliRequirement = cliRequirement
)

# 查看初步的總樣本數
initial_minMaxRes$par[1]
initial_nMax <- initial_minMaxRes$par[1] 

# 固定總樣本數並繼續尋找其他參數
upper2 <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
lower2 <- c(initial_nMax, rep(0.0*pi, nStage - 1), rep(0, nStage))

# 設定PSO配置（進一步搜索）
final_algSetting <- getPSOInfo(
  nSwarm = 256,     # 增加粒子數
  maxIter = 300,    # 增加迭代次數
  psoType = "basic" 
)

# 運行PSO以尋找optimal設計（進一步搜索）
final_OptimRes <- globpso(
  objFunc = kStageOptimObj, PSO_INFO = final_algSetting, 
  lower = lower2, upper = upper2, 
  seed = pso_seed, verbose = TRUE,
  nMin = nMinEachInterim, cliRequirement = cliRequirement
)

# 查看最終的optimal設計結果
final_OptimRes$val     # Objective function value
final_OptimRes$cputime # Computing time

# 轉換最終的PSO結果為可讀的optimal設計
OptimDesign <- kStageFreqCrit(
  nPolarized = final_OptimRes$par[2:nStage], 
  rProportion = final_OptimRes$par[(nStage + 1):length(final_OptimRes$par)], 
  nMax = final_OptimRes$par[1], nMin = nMinEachInterim, cliRequirement
)

# The resulting optimal design
OptimDesign$nseq # Sample sizes at each stage (n_1, ..., n_K)
OptimDesign$rseq # Stopping cutoff sizes at each stage (r_1, ..., r_K)

# Properties of the resulting optimal design
OptimDesign$t1e     # Type I error  
OptimDesign$t2e     # Type II error
OptimDesign$en      # Expected sample size under null hypothesis
OptimDesign$pet_seq # The probabilities of early termination of the trial at each stage
