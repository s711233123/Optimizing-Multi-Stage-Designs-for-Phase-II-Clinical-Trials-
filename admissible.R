# === 設定 ===
library(globpso)
source("run/util.R")
source("run/kStageP2A_Objective.R")
source("run/run_pso_stage.R")

# 固定臨床參數
cli <- list(p0 = 0.05, p1 = 0.25, alpha = 0.05, beta = 0.1)
nStage <- 2
maxIter1 <- 100
maxIter2 <- 500
nSwarm <- 300
nMaxRange <- c(15, 70)

# 要搜尋的 q 值列表
q_list <- seq(0, 1, by = 0.05)

# 結果儲存
result_df <- data.frame()

# 執行主迴圈
for (q in q_list) {
  cat(sprintf(">>> Running q = %.2f\n", q))
  res <- run_pso_two_stage(
    seed = 123 + round(q * 1000),  # 每個 q 用不同 seed
    q_stage1 = 1,
    q_stage2 = q,
    cliRequirement = cli,
    nStage = nStage,
    maxIter1 = maxIter1,
    maxIter2 = maxIter2,
    nSwarm = nSwarm,
    nMaxRange = nMaxRange,
    upper2_is_fixed = FALSE
  )
  
  row <- data.frame(
    q = q,
    n1 = res$nseq[1],
    r1 = res$rseq[1],
    n = res$nseq[2],
    r = res$rseq[2],
    EN = round(res$en, 2),
    t1e = round(res$t1e, 4),
    power = round(1 - res$t2e, 4),
    PET0 = round(res$pet_seq[1], 3)
  )
  
  result_df <- rbind(result_df, row)
}

# 顯示表格
print(result_df)

# 畫圖：EN vs q
library(ggplot2)
ggplot(result_df, aes(x = q, y = EN)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Expected Sample Size vs q",
       x = "q (Trade-off Weight)",
       y = "Expected Sample Size (EN)") +
  theme_minimal(base_size = 13)
