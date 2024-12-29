# Defining sample size calculation function
calculate_sample_size <- function(n1, n2, PET) {
  return(n1 + (1 - PET) * n2)
}

# Defining cumulative binomial distribution function
binomial_cdf <- function(r, p, n) {
  return(pbinom(r, n, p))
}

# Defining binomial distribution probability mass function
binomial_pmf <- function(x, p, n) {
  return(dbinom(x, n, p))
}

# Defining early termination probability function
PET_function <- function(r1, p, n1) {
  return(binomial_cdf(r1, p, n1))
}

# Defining drug rejection probability function
reject_probability <- function(r1, r, p, n1, n2) {
  if (r1 >= n1) {
    return(1)  # 若 r1 超过 n1，概率为 1
  }
  term1 <- binomial_cdf(r1, p, n1)
  term2 <- sum(sapply((r1 + 1):min(n1, r), function(x) {
    binomial_pmf(x, p, n1) * binomial_cdf(r - x, p, n2)
  }))
  prob <- term1 + term2
  return(min(prob, 1))  # 确保概率不超过 1
}


# Defining the fitness function
fitness_function <- function(params, p0, p1, alpha, beta) {
  r1 <- floor(params[1])
  n1 <- floor(params[2])
  r <- floor(params[3])
  n <- floor(params[4])
  n2 <- n - n1
  
  # 检查无效的参数组合
  if (n1 <= 0 || n2 <= 0 || r1 < 0 || r < 0 || r1 > r || r > n) {
    return(Inf)
  }
  
  # 计算早期终止概率
  PET_p0 <- PET_function(r1, p0, n1)
  PET_p1 <- PET_function(r1, p1, n1)
  
  # 计算拒绝药物的概率
  reject_p0 <- reject_probability(r1, r, p0, n1, n2)
  reject_p1 <- 1 - reject_probability(r1, r, p1, n1, n2)
  
  # 打印调试信息
  cat("Parameters: r1 =", r1, "n1 =", n1, "r =", r, "n =", n, "n2 =", n2, "\n")
  cat("PET_p0:", PET_p0, "PET_p1:", PET_p1, "reject_p0:", reject_p0, "reject_p1:", reject_p1, "\n")
  
  # 检查是否满足类型 I 和类型 II 错误概率的约束
  if (reject_p0 > alpha || reject_p1 > beta) {
    return(Inf)  # 如果不满足约束，返回一个大的值
  }
  
  # 计算预期样本大小
  EN_p0 <- calculate_sample_size(n1, n2, PET_p0)
  
  return(EN_p0)  # 返回预期样本大小作为适应度值
}


# Particle Swarm Optimization algorithm
pso <- function(num_particles, num_dimensions, max_iter, p0, p1, alpha, beta) {
  positions <- matrix(runif(num_particles * num_dimensions, 1, 70), num_particles, num_dimensions)
  velocities <- matrix(0, num_particles, num_dimensions)
  
  personal_best_positions <- positions
  personal_best_values <- sapply(1:num_particles, function(i) {
    fitness_function(positions[i, ], p0, p1, alpha, beta)
  })
  global_best_index <- which.min(personal_best_values)
  global_best_position <- personal_best_positions[global_best_index, ]
  global_best_value <- personal_best_values[global_best_index]
  
  for (iter in 1:max_iter) {
    for (i in 1:num_particles) {
      inertia_weight <- runif(1, 0.4, 0.9)
      cognitive_factor <- 2
      social_factor <- 2
      beta1 <- runif(num_dimensions)
      beta2 <- runif(num_dimensions)
      velocity_update <- inertia_weight * velocities[i, ] +
        cognitive_factor * beta1 * (personal_best_positions[i, ] - positions[i, ]) +
        social_factor * beta2 * (global_best_position - positions[i, ])
      positions[i, ] <- positions[i, ] + velocity_update
      velocities[i, ] <- velocity_update
      
      # 施加参数限制
      positions[i, ] <- pmax(pmin(positions[i, ], 50), 1)
      if (positions[i, 1] > positions[i, 3]) {  # 确保 r1 <= r
        positions[i, 1] <- positions[i, 3]
      }
      if (positions[i, 3] > positions[i, 4]) {  # 确保 r <= n
        positions[i, 3] <- positions[i, 4]
      }
      
      current_value <- fitness_function(positions[i, ], p0, p1, alpha, beta)
      if (current_value < personal_best_values[i]) {
        personal_best_positions[i, ] <- positions[i, ]
        personal_best_values[i] <- current_value
      }
    }
    
    new_global_best_index <- which.min(personal_best_values)
    new_global_best_value <- personal_best_values[new_global_best_index]
    if (new_global_best_value < global_best_value) {
      global_best_index <- new_global_best_index
      global_best_position <- personal_best_positions[global_best_index, ]
      global_best_value <- new_global_best_value
    }
  }
  
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
max_iterations <- 100
p0 <- 0.10
p1 <- 0.30
alpha <- 0.10
beta <- 0.10

result <- pso(num_particles, num_dimensions, max_iterations, p0, p1, alpha, beta)
cat("Global best position (r1, n1, r, n):", result$global_best_position, "\n")
cat("Global best value (EN_p0):", result$global_best_value, "\n")
cat("Expected sample size EN(p_0):", result$EN_p0, "\n")
cat("Probability of early termination PET(p_0):", result$PET_p0, "\n")

