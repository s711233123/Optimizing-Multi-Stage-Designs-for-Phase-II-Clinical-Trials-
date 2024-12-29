library(globpso)
library(Rcpp)
library(RcppArmadillo)

cppFunction('
double evaluate_design_cpp(arma::rowvec x, double p1, double p0, double alpha, double beta, int k) {
  arma::ivec r(k);
  arma::ivec n(k);

  for (int i = 0; i < k; ++i) {
    r[i] = round(x[i*2]);
    n[i] = round(x[i*2 + 1]);
  }

  double typeI = R::pbinom(r[0], n[0], p1, 1, 0);
  double typeII = R::pbinom(r[0], n[0], p0, 1, 0);

  for (int j = 1; j < k; ++j) {
    for (int i = 0; i < j; ++i) {
      for (int x_i = r[i] + 1; x_i <= std::min(r[j], n[i]); ++x_i) {
        typeI += R::dbinom(x_i, n[i], p1, 0) * R::pbinom(r[j] - x_i, n[j], p1, 1, 0);
        typeII += R::dbinom(x_i, n[i], p0, 0) * R::pbinom(r[j] - x_i, n[j], p0, 1, 0);
      }
    }
  }

  if ((1 - typeI) < alpha && typeII < beta) {
    double PET_all = R::pbinom(r[0], n[0], p1, 1, 0);
    for (int i = 1; i < k; ++i) {
      double PET_i = 0.0;
      for (int x_i = r[i-1] + 1; x_i <= std::min(r[i], n[i-1]); ++x_i) {
        PET_i += R::dbinom(x_i, n[i-1], p1, 0) * R::pbinom(r[i] - x_i, n[i], p1, 1, 0);
      }
      PET_all += PET_i;
    }

    double EN = n[0];
    for (int i = 1; i < k; ++i) {
      EN += (1 - PET_all) * n[i];
    }
    return EN;
  } else {
    return R_PosInf;
  }
}', depends = "RcppArmadillo")

# PSO適應度函數
pso_fitness <- function(x, p1, p0, alpha, beta, k) {
  return(evaluate_design_cpp(x, p1, p0, alpha, beta, k))
}

# 廣義k階段設計最佳化函數
optimize_k_stage <- function(k, p1, p0, alpha, beta, n) {
  lower_bounds <- rep(c(0, 1), k)
  upper_bounds <- rep(c(n - 1, n - 1), k)
  
  res_c <- globpso(
    objFunc = pso_fitness, 
    lower = lower_bounds, 
    upper = upper_bounds, 
    PSO_INFO = getPSOInfo(maxIter = 1000 * k, nSwarm = 100 * k), 
    p1 = p1, 
    p0 = p0, 
    alpha = alpha, 
    beta = beta, 
    k = k
  )
  
  rounded_par <- round(res_c$par)
  rounded_val <- round(res_c$val, 1)
  
  cat("最佳解參數：", rounded_par, "\n")
  cat("最小期望樣本數：", rounded_val, "\n")
  return(list(par = rounded_par, val = rounded_val))
}

# 設定參數
k <- 2
p1 <- 0.1
p0 <- 0.3
alpha <- 0.1
beta <- 0.1
n <- 40

# 進行優化
result <- optimize_k_stage(k, p1, p0, alpha, beta, n)
