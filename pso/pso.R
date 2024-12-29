# 定义 eta 函数
eta <- function(x, theta) {
  theta[3] * (exp(-theta[2] * x) - exp(-theta[1] * x))
}

# 定义 eta 的梯度
gradient_eta <- function(x, theta) {
  gradient <- numeric(3)
  gradient[1] <- x * theta[3] * exp(-theta[1] * x)
  gradient[2] <- -x * theta[3] * exp(-theta[2] * x)
  gradient[3] <- exp(-theta[2] * x) - exp(-theta[1] * x)
  return(gradient)
}

# 计算信息矩阵
information_matrix <- function(X, weights, theta) {
  num_parameters <- length(theta)
  M <- matrix(0, num_parameters, num_parameters)
  for (i in seq_along(X)) {
    grad <- gradient_eta(X[i], theta)
    M <- M + weights[i] * outer(grad, grad)
  }
  return(M)
}

# 定义适应度函数
fitness_function <- function(X, weights, theta) {
  num_points <- 3  # Number of x points
  X_points <- X[1:num_points]
  M <- information_matrix(X_points, weights, theta)
  det_M <- det(M)
  if (det_M > 0) {
    return(-log(det_M))
  } else {
    return(Inf)  # If the matrix is not positive definite, return a large value
  }
}

# 粒子群优化算法
pso <- function(num_particles, num_dimensions, max_iter, theta) {
  num_points <- 3
  positions <- matrix(runif(num_particles * num_dimensions, 0.1, 1), num_particles, num_dimensions)
  velocities <- matrix(0, num_particles, num_dimensions)
  
  personal_best_positions <- positions
  personal_best_values <- apply(positions, 1, function(position) {
    weights <- position[(num_dimensions - num_points + 1):num_dimensions]
    weights <- weights / sum(weights)  # Normalize weights
    fitness_function(position[1:num_points], weights, theta)
  })
  global_best_index <- which.min(personal_best_values)
  global_best_position <- personal_best_positions[global_best_index, ]
  global_best_value <- personal_best_values[global_best_index]
  
  lower_bound <- c(rep(0.1, num_points), rep(0, num_points))
  upper_bound <- c(rep(30, num_points), rep(1, num_points))
  
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
      
      positions[i, ] <- pmax(pmin(positions[i, ], upper_bound), lower_bound)
      
      weights <- positions[i, (num_dimensions - num_points + 1):num_dimensions]
      weights <- weights / sum(weights)  # Normalize weights
      positions[i, (num_dimensions - num_points + 1):num_dimensions] <- weights  # Update positions with normalized weights
      
      current_value <- fitness_function(positions[i, 1:num_points], weights, theta)
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
  
  return(list(global_best_position = global_best_position, global_best_value = global_best_value))
}

# 示例用法
num_particles <- 100
num_dimensions <- 6  # Three x values and three weights
max_iterations <- 100
theta <- c(0.05884, 4.298, 21.8)  # Example theta values

result <- pso(num_particles, num_dimensions, max_iterations, theta)
cat("Global best position (x1, x2, x3, w1, w2):", result$global_best_position[1:5], "\n")
cat("Global best value:", result$global_best_value, "\n")

