pso <- function(num_particles, num_dimensions, max_iter, p0, p1, alpha, beta) {
  n_fixed <- ceiling(((p0 + p1) / 2) * (1 - ((p0 + p1) / 2)) * ((qnorm(1 - alpha) + qnorm(1 - beta)) / (p1 - p0))^2)
  
  # Initialize positions for N, q1, q2
  positions <- matrix(runif(num_particles * 3), num_particles, 3)
  positions[, 1] <- runif(num_particles, 0, n_fixed)  # Initialize N in the range of [0, n_fixed]
  
  # Initialize w1 and w2 such that they sum to 1
  weights <- matrix(runif(num_particles * 2), num_particles, 2)
  weights <- weights / rowSums(weights)
  
  velocities <- matrix(runif(num_particles * num_dimensions, -1, 1), num_particles, num_dimensions)
  
  personal_best_positions <- cbind(positions, weights)
  personal_best_values <- apply(personal_best_positions, 1, function(pos) fitness_function(pos[1:3], pos[4:5], pos[2:3], p0, p1, alpha, beta))
  global_best_index <- which.min(personal_best_values)
  global_best_position <- personal_best_positions[global_best_index, ]
  global_best_value <- personal_best_values[global_best_index]
  
  for (iter in 1:max_iter) {
    for (i in 1:num_particles) {
      inertia_weight <- runif(1, 0.4, 0.9)
      cognitive_factor <- 2
      social_factor <- 2
      velocity_update <- inertia_weight * velocities[i, ] +
        cognitive_factor * runif(num_dimensions) * (personal_best_positions[i, ] - c(positions[i, ], weights[i, ])) +
        social_factor * runif(num_dimensions) * (global_best_position - c(positions[i, ], weights[i, ]))
      
      positions[i, ] <- positions[i, ] + velocity_update[1:3]
      weights[i, ] <- weights[i, ] + velocity_update[4:5]
      velocities[i, ] <- velocity_update
      
      positions[i, 1] <- pmax(pmin(positions[i, 1], n_fixed), 0)  # Ensure N is within [0, n_fixed]
      positions[i, 2:3] <- pmax(pmin(positions[i, 2:3], 1), 0)  # Ensure q1, q2 are within [0, 1]
      
      # Ensure weights sum to 1
      weights[i, 1] <- pmax(pmin(weights[i, 1], 1), 0)
      weights[i, 2] <- 1 - weights[i, 1]
      
      current_value <- fitness_function(positions[i, ], weights[i, ], positions[i, 2:3], p0, p1, alpha, beta)$value
      if (is.finite(current_value) && current_value < personal_best_values[i]) {
        personal_best_positions[i, ] <- c(positions[i, ], weights[i, ])
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
  
  global_best_position <- round(global_best_position, 2)
  best_fit <- fitness_function(global_best_position[1:3], global_best_position[4:5], global_best_position[2:3], p0, p1, alpha, beta)
  
  return(list(global_best_position = global_best_position,
              global_best_value = global_best_value,
              EN_p0 = best_fit$value,
              PET_p0 = best_fit$PET_p0,
              r1 = best_fit$r1,
              n1 = best_fit$n1,
              r2 = best_fit$r2,
              n2 = best_fit$n2))
}

fitness_function <- function(position, weights, scale_n, p0, p1, alpha, beta) {
  N  <- position[1]
  q1 <- scale_n[1]
  w1 <- weights[1]
  q2 <- scale_n[2]
  w2 <- weights[2]
  
  r1 <- floor(N * w1 * q1)
  n1 <- floor(N * w1)
  r2 <- floor(N * w2 * q2)
  n2 <- floor(N * w2)
  r  <- floor(N * q1 * w1 + N * q2 * w2)
  
  typeI <- pbinom(r1, n1, p0)
  typeII <- pbinom(r1, n1, p1)
  
  if (r1 + 1 <= min(r, n1)) {
    for (x in seq(r1 + 1, min(r, n1), 1)) {
      typeI <- typeI + dbinom(x, n1, p0) * pbinom(r - x, n2, p0)
      typeII <- typeII + dbinom(x, n1, p1) * pbinom(r - x, n2, p1)
    }
  }
  
  PET_p0 <- PET_function(r1, p0, n1)
  
  if ((1 - typeI) < alpha & typeII < beta) {
    value <- n1  + n2 * (1 - pbinom(r1, n1, p0))
  } else {
    value <- Inf
  }
  
  return(list(value = value, PET_p0 = PET_p0, r1 = r1, n1 = n1, r2 = r2, n2 = n2))
}

PET_function <- function(r1, p0, n1) {
  return(pbinom(r1, n1, p0))
}

calculate_type_I_error <- function(position, weights, scale_n, p0) {
  N  <- ceiling(position[1])
  q1 <- scale_n[1]
  w1 <- weights[1]
  q2 <- scale_n[2]
  w2 <- weights[2]
  
  n1 <- ceiling(N * w1)
  n2 <- ceiling(N * w2)
  r1 <- ceiling(n1 * q1)
  r2 <- ceiling(n2 * q2)
  
  typeI <- pbinom(r1, n1, p0)
  if (r1 + 1 <= min(r2, n1)) {
    for (obs in seq(r1 + 1, min(r2, n1), 1)) {
      typeI <- typeI + dbinom(obs, n1, p0) * pbinom(r2 - obs, n2, p0)
    }
  }
  return(1 - typeI)
}

calculate_type_II_error <- function(position, weights, scale_n, p1) {
  N  <- ceiling(position[1])
  q1 <- scale_n[1]
  w1 <- weights[1]
  q2 <- scale_n[2]
  w2 <- weights[2]
  
  n1 <- ceiling(N * w1)
  n2 <- ceiling(N * w2)
  r1 <- ceiling(n1 * q1)
  r2 <- ceiling(n2 * q2)
  
  typeII <- pbinom(r1, n1, p1)
  if (r1 + 1 <= min(r2, n1)) {
    for (obs in seq(r1 + 1, min(r2, n1), 1)) {
      typeII <- typeII + dbinom(obs, n1, p1) * pbinom(r2 - obs, n2, p1)
    }
  }
  return(typeII)
}

num_particles <- 100
num_dimensions <- 5
max_iterations <- 20000
p0 <- 0.10
p1 <- 0.30
alpha <- 0.1
beta <- 0.1

result <- pso(num_particles, num_dimensions, max_iterations, p0, p1, alpha, beta)
cat("Global best position (N, q1, q2, w1, w2):", result$global_best_position, "\n")
cat("Global best value (EN_p0):", result$global_best_value, "\n")
cat("PET_p0:", result$PET_p0, "\n")
cat("r1:", result$r1, "\n")
cat("n1:", result$n1, "\n")
cat("r2:", result$r2, "\n")
cat("n2:", result$n2, "\n")
cat("type_I_error:", calculate_type_I_error(result$global_best_position[1:3], result$global_best_position[4:5], result$global_best_position[2:3], p0), "\n")
cat("type_II_error:", calculate_type_II_error(result$global_best_position[1:3], result$global_best_position[4:5], result$global_best_position[2:3], p1), "\n")
