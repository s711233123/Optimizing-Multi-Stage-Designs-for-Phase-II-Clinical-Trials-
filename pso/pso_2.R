# 定義累積二項分佈函數
binomial_cdf <- function(r, p, n) {
  return(pbinom(r, n, p))
}

# 定義二項分佈概率質量函數
binomial_pmf <- function(x, p, n) {
  return(dbinom(x, n, p))
}

# 定義早期終止概率函數
PET_function <- function(r1, p, n1) {
  return(binomial_cdf(r1, p, n1))
}

# 定義拒絕藥物的概率函數
reject_probability <- function(r1, r, p, n1, n2) {
  if (r1 >= n1) {
    return(1)  # 若r1超過n1，則概率為1
  }
  term1 <- binomial_cdf(r1, p, n1)
  term2 <- sum(sapply((r1 + 1):min(n1, r), function(x) {
    binomial_pmf(x, p, n1) * binomial_cdf(r - x, p, n2)
  }))
  return(term1 + term2)
}

# 定義目標函數: 預期樣本大小
calculate_expected_sample_size <- function(n1, n2, PET) {
  return(n1 + (1 - PET) * n2)
}

# 定義求解最佳兩階段設計的函數
find_optimal_two_stage_design <- function(p0, p1, alpha, beta) {
  # 計算搜索範圍的上界
  p_bar <- (p0 + p1) / 2
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(1 - beta)
  upper_bound <- p_bar * (1 - p_bar) * ((z_alpha + z_beta) / (p1 - p0))^2
  
  # 初始化最小預期樣本大小和相應的設計參數
  min_expected_sample_size <- Inf
  optimal_design <- NULL
  
  # 遍歷所有可能的總樣本n和第一階段樣本n1
  for (n in seq(ceiling(upper_bound), 100)) {
    for (n1 in 1:(n - 1)) {
      # 對每個可能的總樣本n和第一階段樣本n1，通過枚舉計算最佳的r和r1的整數值
      for (r1 in 0:(n1 - 1)) {
        # 獲取第二階段樣本n2
        n2 <- n - n1
        
        # 計算早期終止概率
        PET_p0 <- PET_function(r1, p0, n1)
        
        # 計算拒絕藥物的概率
        reject_p0 <- reject_probability(r1, r1, p0, n1, n2)
        
        # 檢查是否滿足型I和型II錯誤概率的限制
        if (reject_p0 <= alpha) {
          # 如果滿足限制，計算預期樣本大小
          EN_p0 <- calculate_expected_sample_size(n1, n2, PET_p0)
          
          # 更新最小預期樣本大小和相應的設計參數
          if (EN_p0 < min_expected_sample_size) {
            min_expected_sample_size <- EN_p0
            optimal_design <- c(r1 = r1, n1 = n1, r = r1, n = n)
          }
        }
      }
    }
  }
  
  return(optimal_design)
}

# 示例用法
p0 <- 0.10
p1 <- 0.30
alpha <- 0.10
beta <- 0.10

optimal_design <- find_optimal_two_stage_design(p0, p1, alpha, beta)
cat("Optimal two-stage design parameters(r1, n1, r, n):", optimal_design, "\n")

