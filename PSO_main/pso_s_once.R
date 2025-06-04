setwd("~/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")
library(globpso)
source("run/util.R")
source("run/kStageP2A_Objective.R")
source("run/run_pso_stage.R") 


# 自動建立資料夾
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("csv")) dir.create("csv")

# === 參數設定 ===
nStage <- 2
nMaxRange <- c(10, 70)
seed <- 1
maxIter <- 600

p_combinations <- data.frame(
  p0 = seq(0.05, 0.75, by = 0.05),
  p1 = seq(0.25, 0.95, by = 0.05)
)

results <- list()

# === 主迴圈 ===
for (i in 1:nrow(p_combinations)) {
  p0 <- p_combinations$p0[i]
  p1 <- p_combinations$p1[i]
  cat(sprintf("=== Running PSO Single: p0 = %.2f, p1 = %.2f ===\n", p0, p1))
  cliRequirement <- list(p0 = p0, p1 = p1, alpha = 0.1, beta = 0.1)
  
  # optimal: q = 0
  opt_design <- run_pso_single_stage(
    seed = seed,
    q_val = 0,
    cliRequirement = cliRequirement,
    nStage = nStage,
    maxIter = maxIter,
    nMaxRange = nMaxRange
  )
  
  # minimax: q = 1
  min_design <- run_pso_single_stage(
    seed = seed,
    q_val = 1,
    cliRequirement = cliRequirement,
    nStage = nStage,
    maxIter = maxIter,
    nMaxRange = nMaxRange
  )
  
  results[[i]] <- list(
    p0 = p0,
    p1 = p1,
    optimal = opt_design,
    minimax = min_design
  )
}

# === 動態整理 dataframe ===
opt_df <- do.call(rbind, lapply(results, function(res) {
  row <- data.frame(p0 = res$p0, p1 = res$p1)
  for (k in 1:nStage) {
    row[[paste0("stage", k)]] <- paste0(res$optimal$rseq[k], "/", res$optimal$nseq[k])
  }
  row$EN_opt <- res$optimal$en
  for (k in 1:(nStage - 1)) {
    row[[paste0("pet", k, "_opt")]] <- formatC(res$optimal$pet_seq[k], format = "f", digits = 6)
  }
  row
}))

min_df <- do.call(rbind, lapply(results, function(res) {
  row <- data.frame(p0 = res$p0, p1 = res$p1)
  for (k in 1:nStage) {
    row[[paste0("stage", k)]] <- paste0(res$minimax$rseq[k], "/", res$minimax$nseq[k])
  }
  row$EN_minmax <- res$minimax$en
  for (k in 1:(nStage - 1)) {
    row[[paste0("pet", k, "_min")]] <- formatC(res$minimax$pet_seq[k], format = "f", digits = 6)
  }
  row
}))

print(opt_df)
print(min_df)



# === 輸出 ===
# 動態產出檔名（建議用變數）
method <- "pso"            # 或 "ga"
mode   <- "single"         # 或 "two"
opt_file <- sprintf("csv/opt_df_%s_%s_%dstage.csv", method, mode, nStage)
min_file <- sprintf("csv/min_df_%s_%s_%dstage.csv", method, mode, nStage)

# === 輸出 ===
write.csv(opt_df, opt_file, row.names = FALSE)
write.csv(min_df, min_file, row.names = FALSE)


