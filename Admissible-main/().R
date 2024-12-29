# 安裝和加載 "globpso" 套件，並導入兩個本地 R 文件 "util.R" 和 "kStageP2A_Objective.R"
library(globpso)
source("util.R")
source("kStageP2A_Objective.R")

# 設置臨床試驗的要求
cliRequirement <- list(
  p0 = 0.4,    # 零假設的反應率
  p1 = 0.6,    # 替代假設的反應率
  alpha = 0.1, # I 型錯誤的上限
  beta = 0.1   # II 型錯誤的上限
)

# 設置階段數
nStage <- 4

# 設置臨床試驗設計的約束條件
n1Min <- 1      # 第一階段的最小樣本數
nrMin <- 1      # 之後階段的最小樣本數
nMaxRange <- c(10, 70) # 總樣本數的範圍

# 生成每階段的最小樣本數向量
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1))

# 設置PSO的配置
algSetting <- getPSOInfo(
  nSwarm = 256,     # 蜂群大小
  maxIter = 2000,   # 最大迭代次數
  psoType = "basic" # PSO類型 ("basic"、"quantum" 或 "cso")
)

# 設置PSO搜尋的上下限
upper <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0*pi, nStage - 1), rep(0, nStage))

# 設置重複次數
nRepeat <- 10

# 儲存結果的向量
optim_en_results <- numeric(nRepeat)
minmax_en_results <- numeric(nRepeat)

# 迴圈跑10次
for (i in 1:nRepeat) {
  
  ### -------------------------------------------------------- ###
  ### 尋找最佳設計 (Optimal Design)
  ### -------------------------------------------------------- ###
  optimRes  <- globpso(objFunc = kStageOptimObj, PSO_INFO = algSetting, 
                       lower = lower, upper = upper, 
                       seed = NULL, verbose = TRUE,
                       nMin = nMinEachInterim, cliRequirement = cliRequirement)
  
  # 將最佳設計結果轉換為可讀格式
  optimDesign <- kStageFreqCrit(
    nPolarized = optimRes$par[2:nStage],  
    rProportion = optimRes$par[(nStage + 1):length(optimRes$par)], 
    nMax = optimRes$par[1], nMin = nMinEachInterim, cliRequirement)
  
  # 儲存最佳設計的預期樣本數
  optim_en_results[i] <- optimDesign$en
  
  ### -------------------------------------------------------- ###
  ### 尋找最小最大設計 (Minimax Design)
  ### -------------------------------------------------------- ###
  minMaxRes <- globpso(objFunc = kStageMinMaxObj, PSO_INFO = algSetting, 
                       lower = lower, upper = upper, 
                       seed = NULL, verbose = TRUE,
                       nMin = nMinEachInterim, cliRequirement = cliRequirement)
  
  # 將最小最大設計結果轉換為可讀格式
  minMaxDesign <- kStageFreqCrit(
    nPolarized = minMaxRes$par[2:nStage], 
    rProportion = minMaxRes$par[(nStage + 1):length(minMaxRes$par)], 
    nMax = minMaxRes$par[1], nMin = nMinEachInterim, cliRequirement)
  
  # 儲存最小最大設計的預期樣本數
  minmax_en_results[i] <- minMaxDesign$en
}

# 顯示10次的結果
print("Optimal Design Expected Sample Sizes (EN) in 10 runs:")
print(optim_en_results)

print("Minimax Design Expected Sample Sizes (EN) in 10 runs:")
print(minmax_en_results)
