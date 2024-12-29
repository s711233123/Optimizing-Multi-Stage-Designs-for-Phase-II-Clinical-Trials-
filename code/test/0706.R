# 定義累積二項分佈函數
binomial_cdf <- function(r, p, n) {
  return(pbinom(r, n, p))
}

# 定義二項分佈機率質量函數
binomial_pmf <- function(x, p, n) {
  return(dbinom(x, n, p))
}

# 定義早期終止機率函數
PET_function <- function(r1, p, n1) {
  return(pbinom(r1, n1, p))
}

# 定義拒絕藥物的機率函數
reject_probability <- function(r1, r, p, n1, n) {
  n2 <- n - n1
  term1 <- binomial_cdf(r1, p, n1)
  term2 <- sum(sapply((r1 + 1):min(n1, r), function(x) {
    binomial_pmf(x, p, n1) * binomial_cdf(r - x, p, n2)
  }))
  
  return(term1 + term2)
}

# 定義型I誤差計算函數
calculate_type_I_error <- function(r1, r, p0, n1, n) {
  n2 <- n - n1
  term1 <- sum(sapply((r1 + 1):min(n1, r), function(x) {
    binomial_pmf(x, p0, n1) * (1 - binomial_cdf(r - x, p0, n2))
  }))
  return(term1)
}

# 定義型II誤差計算函數
calculate_type_II_error <- function(r1, r, p1, n1, n) {
  n2 <- n - n1
  term1 <- binomial_cdf(r1, p1, n1)
  term2 <- sum(sapply((r1 + 1):min(n1, r), function(x) {
    binomial_pmf(x, p1, n1) * binomial_cdf(r - x, p1, n2)
  }))
  return(term1 + term2)
}
# 測試函數
# print(PET_function(1, 0.1, 12))
# print(reject_probability(3, 5, 0.5, 10, 10))

# 定義目標函數: 預期樣本大小
calculate_expected_sample_size <- function(n1, n2, PET) {
  return(n1 + (1 - PET) * n2)
}

fitness_function <- function(params, p0, p1, alpha, beta) {
  r1 <- floor(params[1])
  n1 <- floor(params[2])
  r <- floor(params[3])
  n <- floor(params[4])
  n2 <- n - n1
  
  # 確保參數的合理性
  if (n1 <= 0 || n2 <= 0 || r1 < 0 || r < r1 || n < r || n < n1 || n < n2) {
    return(Inf)
  }
  
  # 計算類型I和類型II錯誤
  type_I_error <- calculate_type_I_error(r1, r, p0, n1, n)
  type_II_error <- calculate_type_II_error(r1, r, p1, n1, n)
  
  if (type_I_error > alpha || type_II_error > beta) {
    return(Inf)
  }
  
  # 計算期望樣本量
  PET_p0 <- PET_function(r1, p0, n1)
  EN_p0 <- calculate_expected_sample_size(n1, n2, PET_p0)
  
  return(EN_p0)
}


# PSO算法
pso <- function(num_particles, num_dimensions, max_iter, p0, p1, alpha, beta) {
  # 初始化粒子位置和速度
  n_fixed <- ceiling(((p0 + p1) / 2) * (1 - ((p0 + p1) / 2)) * ((qnorm(1 - alpha) + qnorm(1 - beta)) / (p1 - p0))^2)
  positions <- matrix(runif(num_particles * num_dimensions, 1, n_fixed + 20), num_particles, num_dimensions)
  velocities <- matrix(runif(num_particles * num_dimensions, -1, 1), num_particles, num_dimensions)
  
  # 確保初始位置的合理性
  positions <- t(apply(positions, 1, function(pos) {
    pos[2] <- min(pos[2], pos[4] - 1)  # n1 < n
    pos[3] <- min(pos[3], pos[4])      # r <= n
    pos[1] <- min(pos[1], pos[3])      # r1 <= r
    pos[1] <- min(pos[1], pos[2])      # r1 <= n1
    return(pos)
  }))
  
  # 記錄個體最佳位置和值
  personal_best_positions <- positions
  personal_best_values <- apply(positions, 1, fitness_function, p0 = p0, p1 = p1, alpha = alpha, beta = beta)
  
  # 記錄全局最佳位置和值
  global_best_index <- which.min(personal_best_values)
  global_best_position <- personal_best_positions[global_best_index, ]
  global_best_value <- personal_best_values[global_best_index]
  
  # 開始迭代
  for (iter in 1:max_iter) {
    for (i in 1:num_particles) {
      # 更新粒子速度和位置
      inertia_weight <- runif(1, 0.4, 0.9)
      cognitive_factor <- 2
      social_factor <- 2
      velocity_update <- inertia_weight * velocities[i, ] +
        cognitive_factor * runif(1) * (personal_best_positions[i, ] - positions[i, ]) +
        social_factor * runif(1) * (global_best_position - positions[i, ])
      
      positions[i, ] <- positions[i, ] + velocities[i, ]
      velocities[i, ] <- velocity_update
      
      # 確保位置的合理性
      positions[i, ] <- pmax(pmin(positions[i, ], n_fixed + 20), 1)
      positions[i, 2] <- min(positions[i, 2], positions[i, 4] - 1)  # n1 < n
      positions[i, 3] <- min(positions[i, 3], positions[i, 4])      # r <= n
      positions[i, 1] <- min(positions[i, 1], positions[i, 3])      # r1 <= r
      positions[i, 1] <- min(positions[i, 1], positions[i, 2])      # r1 <= n1
      
      # 更新個體最佳位置和值
      current_value <- fitness_function(positions[i, ], p0, p1, alpha, beta)
      if (current_value < personal_best_values[i]) {
        personal_best_positions[i, ] <- positions[i, ]
        personal_best_values[i] <- current_value
      }
    }
    
    # 更新全局最佳位置和值
    new_global_best_index <- which.min(personal_best_values)
    new_global_best_value <- personal_best_values[new_global_best_index]
    if (new_global_best_value < global_best_value) {
      global_best_index <- new_global_best_index
      global_best_position <- personal_best_positions[global_best_index, ]
      global_best_value <- new_global_best_value
    }
  }
  
  # 返回結果
  best_fit <- list(EN_p0 = global_best_value,
                   PET_p0 = PET_function(floor(global_best_position[1]), p0, floor(global_best_position[2])))
  
  return(list(global_best_position = global_best_position, 
              global_best_value = global_best_value,
              EN_p0 = best_fit$EN_p0,
              PET_p0 = best_fit$PET_p0))
}


# 示例用法
num_particles <- 100
num_dimensions <- 4  # r1, n1, r, n
max_iterations <- 2000
p0 <- 0.10
p1 <- 0.30
alpha <- 0.10
beta <- 0.10

result <- pso(num_particles, num_dimensions, max_iterations, p0, p1, alpha, beta)
cat("Global best position (r1, n1, r, n):", result$global_best_position, "\n")
cat("Global best value (EN_p0):", result$global_best_value, "\n")
cat("PET_p0:", result$PET_p0, "\n")
cat("type_I_error:", calculate_type_I_error(floor(result$global_best_position[1]), floor(result$global_best_position[3]), p0, floor(result$global_best_position[2]), floor(result$global_best_position[4])), "\n")
cat("type_II_error:", calculate_type_II_error(floor(result$global_best_position[1]), floor(result$global_best_position[3]), p1, floor(result$global_best_position[2]), floor(result$global_best_position[4])), "\n")
