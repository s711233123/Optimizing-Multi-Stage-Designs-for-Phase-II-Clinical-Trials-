library(globpso)
library(Rcpp)
library(RcppArmadillo)

# 使用Rcpp和RcppArmadillo定義C++函數
cppFunction('double evaluate_design_cpp(arma::rowvec x, double p1, double p0, double alpha, double beta) {
  int r1 = round(x[0]);
  int n1 = round(x[1]);
  int r2 = round(x[2]);
  int n2 = round(x[3]);
  
  // 計算第一類型錯誤率
  double typeI = R::pbinom(r1, n1, p1, 1, 0);
  double typeII = R::pbinom(r1, n1, p0, 1, 0);
  if (r1 + 1 <= std::min(r2, n1)) {
    for (int obs = r1 + 1; obs <= std::min(r2, n1); obs++) {
      typeI += R::dbinom(obs, n1, p1, 0) * R::pbinom(r2 - obs, n2, p1, 1, 0);
      typeII += R::dbinom(obs, n1, p0, 0) * R::pbinom(r2 - obs, n2, p0, 1, 0);
    }
  }
  
  // 檢查是否滿足錯誤率約束
  if ((1 - typeI) < alpha && typeII < beta) {
    // 計算早期終止概率
    double prob_et_null = R::pbinom(r1, n1, p1, 1, 0);
    double prob_et_alt = R::pbinom(r1, n1, p0, 1, 0);
    // 計算期望樣本數
    double en_p0 = n1 + n2 * (1 - prob_et_null);
    return en_p0;
  } else {
    return R_PosInf;  // 不滿足約束條件返回一個非常大的值
  }
}', depends = "RcppArmadillo")

# PSO適應度函數
pso_fitness <- function(x, p1, p0, alpha, beta) {
  return(evaluate_design_cpp(x, p1, p0, alpha, beta))
}

# PSO參數
p1 <- 0.1
p0 <- 0.3
alpha <- 0.1
beta <- 0.1
n <- 40

# 搜索範圍
lower_bounds <- c(0, 1, 0, 1)
upper_bounds <- c(n - 1, n - 1, n - 1, n - 1)

# 進行優化
# 進行優化
res_c <- globpso(
  objFunc = function(x) pso_fitness(x, p1, p0, alpha, beta), 
  lower = lower_bounds, 
  upper = upper_bounds, 
  PSO_INFO = getPSOInfo(maxIter = 1000, nSwarm = 100)
)

# 四捨五入
rounded_par <- round(res_c$par)
rounded_val <- round(res_c$val,1) 


# 結果
cat("最佳解參數：", rounded_par, "\n")
cat("最小期望樣本數：", rounded_val, "\n")
  
