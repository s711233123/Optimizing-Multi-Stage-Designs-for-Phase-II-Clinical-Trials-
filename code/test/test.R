# PSO算法
pso <- function(num_particles, num_dimensions, max_iter, p0, p1, alpha, beta) {
  # 初始化粒子位置和速度
  n_fixed <- ceiling(((p0 + p1) / 2) * (1 - ((p0 + p1) / 2)) * ((qnorm(1 - alpha) + qnorm(1 - beta)) / (p1 - p0))^2)
  positions <- matrix(runif(num_particles * num_dimensions, 1, n_fixed + 20), num_particles, num_dimensions)
  velocities <- matrix(runif(num_particles * num_dimensions, -1, 1), num_particles, num_dimensions)
  
  # 確保初始位置的合理性
  positions <- t(apply(positions, 1, function(pos) {
    pos[2] <- min(pos[2], pos[4])      # n1 <= n
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
      positions[i, ] <- pmax(pmin(positions[i, ], n_fixed + 20), 0.1)
      positions[i, 2] <- min(positions[i, 2], positions[i, 4])  # n1 < n
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
  
  # 將最佳位置取整
  
  best_fit <- list(EN_p0 = global_best_value,
                   PET_p0 = PET_function(global_best_position[1], p0, global_best_position[2]))
  
  return(list(global_best_position = global_best_position, 
              global_best_value = global_best_value,
              EN_p0 = best_fit$EN_p0,
              PET_p0 = best_fit$PET_p0))
}

# 適應度函數
fitness_function <- function(position, p0, p1, alpha, beta) {
  r1 <- floor(position[1])
  n1 <- floor(position[2])
  r <- floor(position[3])
  n <- floor(position[4])
  
  typeI <- pbinom(r1, n1, p0)
  typeII <- pbinom(r1, n1, p1)
  
  if (r1 + 1 <= min(r, n1)) {
    for (obs in seq(r1 + 1, min(r, n1), 1)) {
      typeI <- typeI + dbinom(obs, n1, p0) * pbinom(r - obs, n - n1, p0)
      typeII <- typeII + dbinom(obs, n1, p1) * pbinom(r - obs, n - n1, p1)
    }
  }
  
  if ((1 - typeI) < alpha & typeII < beta) {
    return(n1 * pbinom(r1, n1, p0) + n * (1 - pbinom(r1, n1, p0)))
  } else {
    return(Inf)
  }
}

# 提前終止概率（PET）函數
PET_function <- function(r1, p0, n1) {
  return(pbinom(r1, n1, p0))
}

# 計算I型錯誤率
calculate_type_I_error <- function(r1, r, p0, n1, n) {
  typeI <- pbinom(r1, n1, p0)
  if (r1 + 1 <= min(r, n1)) {
    for (obs in seq(r1 + 1, min(r, n1), 1)) {
      typeI <- typeI + dbinom(obs, n1, p0) * pbinom(r - obs, n - n1, p0)
    }
  }
  return(1 - typeI)
}

# 計算II型錯誤率
calculate_type_II_error <- function(r1, r, p1, n1, n) {
  typeII <- pbinom(r1, n1, p1)
  if (r1 + 1 <= min(r, n1)) {
    for (obs in seq(r1 + 1, min(r, n1), 1)) {
      typeII <- typeII + dbinom(obs, n1, p1) * pbinom(r - obs, n - n1, p1)
    }
  }
  return(typeII)
}

# 示例用法
num_particles <- 300
num_dimensions <- 4  # r1, n1, r, n
max_iterations <- 200000
p0 <- 0.2
p1 <- 0.4
alpha <- 0.05
beta <- 0.1

result <- pso(num_particles, num_dimensions, max_iterations, p0, p1, alpha, beta)
cat("Global best position (r1, n1, r, n):", result$global_best_position, "\n")
cat("Global best value (EN_p0):", result$global_best_value, "\n")
cat("PET_p0:", result$PET_p0, "\n")
cat("type_I_error:", calculate_type_I_error(result$global_best_position[1], result$global_best_position[3], p0, result$global_best_position[2], result$global_best_position[4]), "\n")
cat("type_II_error:", calculate_type_II_error(result$global_best_position[1], result$global_best_position[3], p1, result$global_best_position[2], result$global_best_position[4]), "\n")
