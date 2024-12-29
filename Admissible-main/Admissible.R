### ======================================================== ###
### Example Code
### ------------
### Using PSO to generate K-stage clinical trial design
### under Simon's frequentist framework
### ======================================================== ###
getwd()
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
nMaxRange <- c(15, 70) # the range of the total sample size
# Generate the constraint vector of minimal sample sizes at each stage
# This is the input of the objective function
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1)) # (do not change)

# Set the PSO configuration
algSetting <- getPSOInfo(
  nSwarm = 300,     # swarm size 
  maxIter = 100,    # number of iterations
  psoType = "basic" # PSO type (one can use "basic" or "quantum" or "cso")
)

# Set Random seed for reproducibility
#pso_seed <- 1

# 設定搜尋範圍
n_range <- 25:32

# 初始化結果存儲
results <- data.frame(
  N = numeric(0),
  EN = numeric(0)
)

# 迭代搜尋不同的 N
for (nMax in n_range) {
  # 更新 PSO 的上下限
  upper <- c(nMax, rep(0.5 * pi, nStage - 1), rep(1, nStage))
  lower <- c(nMax, rep(0.0 * pi, nStage - 1), rep(0, nStage))
  
  # 執行 PSO
  optimRes <- globpso(
    objFunc = kStageOptimObj, PSO_INFO = algSetting,
    lower = lower, upper = upper, 
    seed = NULL, verbose = FALSE,
    nMin = nMinEachInterim, cliRequirement = cliRequirement
  )
  
  # 將最佳設計轉為可讀格式
  optimDesign <- kStageFreqCrit(
    nPolarized = optimRes$par[2:nStage],
    rProportion = optimRes$par[(nStage + 1):length(optimRes$par)],
    nMax = optimRes$par[1], nMin = nMinEachInterim, cliRequirement
  )
  
  # 存儲結果
  results <- rbind(results, data.frame(
    N = nMax,
    EN = optimDesign$en
  ))
}

# 繪製圖表
library(ggplot2)
ggplot(results, aes(x = N, y = EN)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(
    title = "Expected Sample Size under Null Hypothesis (EN) vs N",
    x = "Maximal Sample Size (N)",
    y = "Expected Sample Size (EN)"
  ) +
  theme_minimal()

