library(pso)

# 適應度函數：評估設計參數是否滿足約束並計算期望樣本數
evaluate_design <- function(x, p1, p0, alpha, beta) {
  r1 <- round(x[1])
  n1 <- round(x[2])
  r2 <- round(x[3])
  n2 <- round(x[4])
  
  # 計算第一類型錯誤率
  typeI <- pbinom(r1, n1, p1)
  typeII <- pbinom(r1, n1, p0)
  if (r1 + 1 <= min(r2, n1)) {
    for (obs in seq(r1 + 1, min(r2, n1), 1)) {
      typeI <- typeI + dbinom(obs, n1, p1) * pbinom(r2 - obs, n2, p1)
      typeII <- typeII + dbinom(obs, n1, p0) * pbinom(r2 - obs, n2, p0)
    }
  }
  
  # 檢查是否滿足錯誤率約束
  if ((1 - typeI) < alpha & typeII < beta) {
    # 計算早期終止概率
    prob_et_null <- pbinom(r1, n1, p1)
    prob_et_alt <- pbinom(r1, n1, p0)
    # 計算期望樣本數
    en_p0 <- n1 * prob_et_null + (n1 + n2) * (1 - prob_et_null)
    return(en_p0)
  } else {
    return(Inf)  # 不滿足約束條件返回一個非常大的值
  }
}

# PSO適應度函數
pso_fitness <- function(x) {
  return(evaluate_design(x, p1, p0, alpha, beta))
}

# PSO參數
p1 <- 0.2
p0 <- 0.4
alpha <- 0.05
beta <- 0.1
n <- 40

# 搜索範圍
lower_bounds <- c(0, 1, 0, 1)
upper_bounds <- c(n - 1, n - 1, n - 1, n - 1)

# PSO優化
pso_result <- psoptim(par = rep(NA, 4), fn = pso_fitness, lower = lower_bounds, upper = upper_bounds, control = list(maxit = 1000, s = 100))

# 打印結果並四捨五入到整數
rounded_par <- round(pso_result$par)
print(rounded_par)
print(round(pso_result$value, 1))

