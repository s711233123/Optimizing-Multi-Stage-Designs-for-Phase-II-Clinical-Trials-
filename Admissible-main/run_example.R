### ======================================================== ###
### Using PSO to generate K-stage clinical trial design
### under Simon's frequentist framework
### ======================================================== ###
# Install and load the "globpso" package
library(globpso)
# Import two local R files "util.R" and "kStageP2A_Objective.R"
source("util.R")
source("kStageP2A_Objective.R")
getwd()
# Set Requirements of the Clinical Trial
cliRequirement <- list(
  p0 = 0.05,    # Response rate in the null hypothesis
  p1 = 0.25,    # Response rate in the alternative hypothesis
  alpha = 0.05, # Upper bound of type I error
  beta = 0.1   # Upper bound of type II error
)

# Set the required number of stages
nStage <- 2

# Set Constraints of the Clinical Trial Design
n1Min <- 1       # Minimal sample size at the first stage
nrMin <- 1       # Minimal sample size at stages after the 1st one
nMaxRange <- c(10, 70) # The range of the total sample size
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1)) # Minimal sample sizes at each stage

# Set the PSO configuration
algSetting <- getPSOInfo(
  nSwarm = 300,     # Swarm size 
  maxIter = 200,    # Number of iterations
  psoType = "basic" # PSO type: "basic", "quantum", or "cso"
)

# Set the lower and upper bounds for PSO search
upper <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0*pi, nStage - 1), rep(0, nStage))


### -------------------------------------------------------- ###
### Find Optimal Design
### -------------------------------------------------------- ###
optimRes <- globpso(
  objFunc = kStageOptimObj, PSO_INFO = algSetting,
  lower = lower, upper = upper, 
  seed = NULL, verbose = TRUE,
  nMin = nMinEachInterim, cliRequirement = cliRequirement
)
optimDesign <- kStageFreqCrit(
  nPolarized = optimRes$par[2:nStage],  
  rProportion = optimRes$par[(nStage + 1):length(optimRes$par)], 
  nMax = optimRes$par[1], nMin = nMinEachInterim, cliRequirement
)

### -------------------------------------------------------- ###
### Find Minimax Design
### -------------------------------------------------------- ###
minMaxRes <- globpso(
  objFunc = kStageMinMaxObj, PSO_INFO = algSetting,
  lower = lower, upper = upper, 
  seed = NULL, verbose = TRUE,
  nMin = nMinEachInterim, cliRequirement = cliRequirement
)
minMaxDesign <- kStageFreqCrit(
  nPolarized = minMaxRes$par[2:nStage], 
  rProportion = minMaxRes$par[(nStage + 1):length(minMaxRes$par)], 
  nMax = minMaxRes$par[1], nMin = nMinEachInterim, cliRequirement
)

### -------------------------------------------------------- ###
### Generate Results for Multiple q Values
### -------------------------------------------------------- ###
# 設定 q 值範圍
q_values <- seq(0, 1, by = 0.1)  
compromise_results_df <- data.frame() 

# 迭代不同的 q 值，生成折衷設計
for (q in q_values) {
  cat("Running for q =", q, "...\n")
  
  compromiseRes <- globpso(
    objFunc = function(par) kStageCompromiseObj(par, nMinEachInterim, cliRequirement, q),
    PSO_INFO = algSetting,
    lower = lower, upper = upper,
    seed = NULL, verbose = FALSE
  )
  
  compromiseDesign <- kStageFreqCrit(
    nPolarized = compromiseRes$par[2:nStage], 
    rProportion = compromiseRes$par[(nStage + 1):length(compromiseRes$par)], 
    nMax = compromiseRes$par[1], nMin = nMinEachInterim, cliRequirement
  )
  
  # 將結果儲存到 data.frame
  compromise_results_df <- rbind(compromise_results_df, data.frame(
    q = q,
    N = compromiseDesign$nseq[length(compromiseDesign$nseq)],
    EN = compromiseDesign$en
  ))
}

# 檢視結果
print(compromise_results_df)



library(ggplot2)


# 繪製 N 隨 q 的變化
ggplot(compromise_results_df, aes(x = q, y = N)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(title = "Maximal Sample Size (N) vs q",
       x = "Weight q", y = "Maximal Sample Size (N)") +
  theme_minimal()

# 繪製 EN 隨 q 的變化
ggplot(compromise_results_df, aes(x = q, y = EN)) +
  geom_line(color = "green", size = 1) +
  geom_point(color = "black", size = 2) +
  labs(title = "Expected Sample Size (EN) vs q",
       x = "Weight q", y = "Expected Sample Size (EN)") +
  theme_minimal()

