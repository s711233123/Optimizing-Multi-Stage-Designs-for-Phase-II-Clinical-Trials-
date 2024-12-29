library(globpso)
library(Rcpp)
library(RcppArmadillo)

# 使用Rcpp和RcppArmadillo定義C++函數
cppFunction('double evaluate_design_cpp(arma::rowvec x, double p1, double p0, double alpha, double beta) {
  int r1 = round(x[0]);
  int n1 = round(x[1]);
  int r2 = round(x[2]);
  int n2 = round(x[3]);
  int r3 = round(x[4]);
  int n3 = round(x[5]);
  
  // 計算型 I & II 錯誤機率
  double typeI = R::pbinom(r1, n1, p1, 1, 0);
  double typeII = R::pbinom(r1, n1, p0, 1, 0);

  if (r1 + 1 <= std::min(r2, n1)) {
      for (int x1 = r1 + 1; x1 <= std::min(r2, n1); x1++) {
          typeI += R::dbinom(x1, n1, p1, 0) * R::pbinom(r2 - x1, n2, p1, 1, 0);
          typeII += R::dbinom(x1, n1, p0, 0) * R::pbinom(r2 - x1, n2, p0, 1, 0);
      }
  }

  if (r1 + 1 <= std::min(r3, n1)) {
      for (int x1 = r1 + 1; x1 <= std::min(r3, n1); ++x1) {
          if (r2 + 1 - x1 <= std::min(r3 - x1, n2)) {
              for (int x2 = r2 + 1 - x1; x2 <= std::min(r3 - x1, n2); ++x2) {
                  typeI += R::dbinom(x1, n1, p1, 0) *
                           R::dbinom(x2, n2, p1, 0) *
                           R::pbinom(r3 - x1 - x2, n3, p1, 1, 0);
                  typeII += R::dbinom(x1, n1, p0, 0) *
                            R::dbinom(x2, n2, p0, 0) *
                            R::pbinom(r3 - x1 - x2, n3, p0, 1, 0);
              }
          }
      }
  }

  
  // 檢查是否滿足錯誤率約束
  if ((1 - typeI) < alpha && typeII < beta) {
    // 計算第一階段提前終止機率 PET_1
    double PET_1 = R::pbinom(r1, n1, p1, 1, 0);
    
    // 計算第二階段提前終止機率 PET_2
    double PET_2 = 0.0;
    if (r1 + 1 <= std::min(r2, n1)) {
        for (int x1 = r1 + 1; x1 <= std::min(r2, n1); ++x1) {
            PET_2 += R::dbinom(x1, n1, p1, 0) * R::pbinom(r2 - x1, n2, p1, 1, 0);
        }
    }
    
    // 總提前終止機率 PET_all
    double PET_all = PET_1 + PET_2;
    
    // 計算期望樣本量 EN
    double EN = n1 + (1 - PET_1) * n2 + (1 - PET_all) * n3;
    return EN;
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
lower_bound <- c(0, 1, 0, 1, 1, 1)
upper_bound <- c(n - 1, n - 1, n - 1, n - 1, n - 1, n - 1)

# 進行優化
res_c <- globpso(objFunc = pso_fitness, lower = lower_bound, upper = upper_bound, PSO_INFO = getPSOInfo(maxIter = 1200, nSwarm = 500), p1 = p1, p0 = p0, alpha = alpha, beta = beta)

# 四捨五入
rounded_par <- round(res_c$par)
rounded_val <- round(res_c$val, 2)

# 結果
cat("最佳解參數：", rounded_par, "\n")
cat("最小期望樣本數：", rounded_val, "\n")

