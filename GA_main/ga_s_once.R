setwd("~/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")
library(GA)
source("run/util.R")
source("run/kStageP2A_Objective.R")
source("run/run_ga_stage.R")

# === 基本參數設定 ===
nStage <- 3
nMaxRange <- c(10, 70)
seed <- 1
popSize <- 256
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
  cat(sprintf("=== GA Running: p0 = %.2f, p1 = %.2f ===\n", p0, p1))
  
  cliRequirement <- list(p0 = p0, p1 = p1, alpha = 0.1, beta = 0.1)
  
  res_opt <- run_ga_single_stage(
    seed = seed, q_val = 0,
    cliRequirement = cliRequirement,
    nStage = nStage,
    maxIter = maxIter,
    popSize = popSize,
    nMaxRange = nMaxRange
  )
  
  res_min <- run_ga_single_stage(
    seed = seed, q_val = 1,
    cliRequirement = cliRequirement,
    nStage = nStage,
    maxIter = maxIter,
    popSize = popSize,
    nMaxRange = nMaxRange
  )
  
  results[[i]] <- list(
    p0 = p0,
    p1 = p1,
    optimal = res_opt,
    minimax = res_min
  )
}

# === 動態 opt_df ===
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

# === 動態 min_df ===
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


# === 儲存輸出 ===
write.csv(opt_df, "results/opt_df_ga.csv", row.names = FALSE)
write.csv(min_df, "results/min_df_ga.csv", row.names = FALSE)
