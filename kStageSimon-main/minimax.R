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
  p0 = 0.05,    # response rate in the null hypothesis
  p1 = 0.25,    # response rate in the alternative hypothesis
  alpha = 0.05, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

# Set the required number of stages
nStage <- 2

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
  nSwarm = 100,      # 初步設置的粒子數
  maxIter = 50,     # 初步設置的迭代次數
  psoType = "basic" # PSO類型（保持量子模式）
)

# 運行PSO以尋找minimax設計（初步搜索）
initial_minMaxRes <- globpso(
  objFunc = kStageMinMaxObj, PSO_INFO = initial_algSetting, 
  lower = lower, upper = upper, 
  seed = NULL, verbose = TRUE,
  nMin = nMinEachInterim, cliRequirement = cliRequirement
)

# 查看初步的總樣本數
initial_minMaxRes$par[1]
initial_nMax <- initial_minMaxRes$par[1] 

# 固定總樣本數並繼續尋找其他參數
upper2 <- c(initial_nMax, rep(0.5*pi, nStage - 1), rep(1, nStage))
lower2 <- c(initial_nMax, rep(0.0*pi, nStage - 1), rep(0, nStage))

# 設定PSO配置（進一步搜索）
final_algSetting <- getPSOInfo(
  nSwarm = 100,     # 增加粒子數
  maxIter = 100,    # 增加迭代次數
  psoType = "basic" 
)

# 運行PSO以尋找minimax設計（進一步搜索）
final_minMaxRes <- globpso(
  objFunc = kStageMinMaxObj, PSO_INFO = final_algSetting, 
  lower = lower2, upper = upper2, 
  seed = NULL, verbose = TRUE,
  nMin = nMinEachInterim, cliRequirement = cliRequirement
)

# 查看最終的minimax設計結果
final_minMaxRes$val     # Objective function value
final_minMaxRes$cputime # Computing time

# 轉換最終的PSO結果為可讀的minimax設計
minMaxDesign <- kStageFreqCrit(
  nPolarized = final_minMaxRes$par[2:nStage], 
  rProportion = final_minMaxRes$par[(nStage + 1):length(final_minMaxRes$par)], 
  nMax = final_minMaxRes$par[1], nMin = nMinEachInterim, cliRequirement
)

# The resulting minimax design
minMaxDesign$nseq # Sample sizes at each stage (n_1, ..., n_K)
minMaxDesign$rseq # Stopping cutoff sizes at each stage (r_1, ..., r_K)

# Properties of the resulting minimax design
minMaxDesign$t1e     # Type I error
minMaxDesign$t2e     # Type II error
minMaxDesign$en      # Expected sample size under null hypothesis
minMaxDesign$pet_seq # The probabilities of early termination of the trial at each stage
