setwd("~/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")
library(ABCoptim)
source("run/util.R")
source("run/kStageP2A_Objective.R")
source("run/run_abc_stage.R")  # 包含 run_abc_two_stage_fixed()

# 建立資料夾
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("csv")) dir.create("csv")

# === 參數設定 ===
nStage <- 3
nMaxRange <- c(10, 70)
seed <- 1 
maxIter1 <- 100
maxIter2 <- 500

p_combinations <- data.frame(
  p0 = seq(0.05, 0.75, by = 0.05),
  p1 = seq(0.25, 0.95, by = 0.05)
)


results <- list()

# === 主迴圈 ===
for (i in 1:nrow(p_combinations)) {
  p0 <- p_combinations$p0[i]
  p1 <- p_combinations$p1[i]
  cat(sprintf("=== ABC Running: p0 = %.2f, p1 = %.2f ===\n", p0, p1))
  cliRequirement <- list(p0 = p0, p1 = p1, alpha = 0.1, beta = 0.1)
  
  opt_design <- run_abc_two_stage(
    seed = seed, q_stage1 = 1, q_stage2 = 0,
    cliRequirement = cliRequirement,
    nStage = nStage,
    maxIter1 = maxIter1, maxIter2 = maxIter2,
    foodNumber = 256,
    nMaxRange = nMaxRange,
    upper2_is_fixed = FALSE
  )
  
  min_design <- run_abc_two_stage(
    seed = seed, q_stage1 = 1, q_stage2 = 1,
    cliRequirement = cliRequirement,
    nStage = nStage,
    maxIter1 = maxIter1, maxIter2 = maxIter2,
    foodNumber = 256,
    nMaxRange = nMaxRange,
    upper2_is_fixed = TRUE
  )
  
  if (!is.null(opt_design) && !is.null(min_design)) {
    results[[length(results) + 1]] <- list(
      p0 = p0,
      p1 = p1,
      optimal = opt_design,
      minimax = min_design
    )
  }
}

# === 整理結果（opt_df）===
opt_list <- lapply(results, function(res) {
  if (is.null(res$optimal)) return(NULL)
  row <- data.frame(p0 = res$p0, p1 = res$p1)
  for (k in 1:nStage) {
    row[[paste0("stage", k)]] <- paste0(res$optimal$rseq[k], "/", res$optimal$nseq[k])
  }
  row$EN_opt <- res$optimal$en
  for (k in 1:(nStage - 1)) {
    row[[paste0("pet", k, "_opt")]] <- formatC(res$optimal$pet_seq[k], format = "f", digits = 6)
  }
  row
})
opt_df <- do.call(rbind, Filter(Negate(is.null), opt_list))

# === 整理結果（min_df）===
min_list <- lapply(results, function(res) {
  if (is.null(res$minimax)) return(NULL)
  row <- data.frame(p0 = res$p0, p1 = res$p1)
  for (k in 1:nStage) {
    row[[paste0("stage", k)]] <- paste0(res$minimax$rseq[k], "/", res$minimax$nseq[k])
  }
  row$EN_minmax <- res$minimax$en
  for (k in 1:(nStage - 1)) {
    row[[paste0("pet", k, "_min")]] <- formatC(res$minimax$pet_seq[k], format = "f", digits = 6)
  }
  row
})
min_df <- do.call(rbind, Filter(Negate(is.null), min_list))

# === 印出結果 ===
print(opt_df)
print(min_df)


