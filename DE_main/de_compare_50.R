setwd("~/Optimizing-Multi-Stage-Designs-for-Phase-II-Clinical-Trials-")
library(DEoptim)
source("run/util.R")
source("run/kStageP2A_Objective.R")
source("run/run_de_stage.R")

nStage <- 3
nMaxRange <- c(10, 70)
seeds <- 1:50

p_combinations <- data.frame(
  p0 = seq(0.05, 0.75, by = 0.05),
  p1 = seq(0.25, 0.95, by = 0.05)
)

optimal_vals <- c(13.78, 17.79, 21.11, 23.86, 25.55, 26.75, 27.84,
                  27.59, 27.22, 26.08, 24.36, 22.3, 19.23, 15.45, 10.94)
minimax_vals <- c(15.86, 19.10, 22.23, 26.41, 29.83, 33.71, 28.59,
                  31.56, 28.94, 29.00, 27.74, 22.97, 19.89, 17.85, 10.94)

iter_combinations <- list(
  single = list(maxIter = 600),
  two_stage = list(
    c(100, 500), c(100, 400), c(100, 300), c(100, 200), c(100, 100)
  )
)

all_results <- list()

for (i in 1:nrow(p_combinations)) {
  p0 <- p_combinations$p0[i]
  p1 <- p_combinations$p1[i]
  cat(sprintf("\n========== DE Start p0 = %.2f, p1 = %.2f ==========\n", p0, p1))
  cliRequirement <- list(p0 = p0, p1 = p1, alpha = 0.1, beta = 0.1)
  results <- list()
  
  for (seed in seeds) {
    cat(sprintf("Running DE: p0 = %.2f, p1 = %.2f, seed = %d\n", p0, p1, seed))
    result_seed <- list()
    
    maxIter_single <- iter_combinations$single$maxIter
    result_seed[[paste0("single_opt_", maxIter_single)]] <- run_de_single_stage(seed, 0, cliRequirement, nStage, maxIter_single, nMaxRange)
    result_seed[[paste0("single_min_", maxIter_single)]] <- run_de_single_stage(seed, 1, cliRequirement, nStage, maxIter_single, nMaxRange)
    
    for (combo in iter_combinations$two_stage) {
      maxIter1 <- combo[1]
      maxIter2 <- combo[2]
      result_seed[[paste0("opt_", maxIter1, "_", maxIter2)]] <-
        run_de_two_stage(seed, 1, 0, cliRequirement, nStage, maxIter1, maxIter2, nMaxRange, upper2_is_fixed = FALSE)
      result_seed[[paste0("min_", maxIter1, "_", maxIter2)]] <-
        run_de_two_stage(seed, 1, 1, cliRequirement, nStage, maxIter1, maxIter2, nMaxRange, upper2_is_fixed = TRUE)
    }
    results[[as.character(seed)]] <- result_seed
  }
  
  ref_vals <- list(
    opt = optimal_vals[i],
    min = minimax_vals[i],
    single_opt = optimal_vals[i],
    single_min = minimax_vals[i]
  )
  
  summary_df <- data.frame(Seed = seeds)
  for (seed in seeds) {
    for (key in names(results[[as.character(seed)]])) {
      prefix <- sub("_\\d+.*", "", key)
      target_val <- ref_vals[[prefix]]
      t_col <- paste0("Time_", key)
      s_col <- paste0("Success_", key)
      if (!(t_col %in% names(summary_df))) summary_df[[t_col]] <- NA_real_
      if (!(s_col %in% names(summary_df))) summary_df[[s_col]] <- NA_real_
      en_val <- round(results[[as.character(seed)]][[key]]$en, 2)
      summary_df[summary_df$Seed == seed, t_col] <- results[[as.character(seed)]][[key]]$cputime
      summary_df[summary_df$Seed == seed, s_col] <- en_val <= target_val + 0.01
    }
  }
  
  success_rates <- colMeans(summary_df[, grep("Success_", names(summary_df))])
  avg_times <- colMeans(summary_df[, grep("Time_", names(summary_df))])
  all_results[[paste0("p0_", p0, "_p1_", p1)]] <- list(
    success_rates = success_rates,
    avg_times = avg_times,
    raw_results = results
  )
}

saveRDS(all_results, file = "results/compare_results_de.rds")


# === 整理結果 ===
opt_summary_df <- data.frame()
min_summary_df <- data.frame()

for (key in names(all_results)) {
  parts <- strsplit(key, "_")[[1]]
  p0 <- as.numeric(parts[2])
  p1 <- as.numeric(parts[4])
  res <- all_results[[key]]
  
  row_opt <- data.frame(p0 = p0, p1 = p1)
  row_min <- data.frame(p0 = p0, p1 = p1)
  
  # OPT
  opt_time_keys <- grep("^Time_.*(opt|single_opt)", names(res$avg_times), value = TRUE)
  opt_succ_keys <- grep("^Success_.*(opt|single_opt)", names(res$success_rates), value = TRUE)
  for (k in opt_time_keys) row_opt[[k]] <- res$avg_times[[k]]
  for (k in opt_succ_keys) row_opt[[k]] <- res$success_rates[[k]]
  opt_summary_df <- rbind(opt_summary_df, row_opt)
  
  # MIN
  min_time_keys <- grep("^Time_.*(min|single_min)", names(res$avg_times), value = TRUE)
  min_succ_keys <- grep("^Success_.*(min|single_min)", names(res$success_rates), value = TRUE)
  for (k in min_time_keys) row_min[[k]] <- res$avg_times[[k]]
  for (k in min_succ_keys) row_min[[k]] <- res$success_rates[[k]]
  min_summary_df <- rbind(min_summary_df, row_min)
}

# 欄位名稱轉簡潔
names(opt_summary_df) <- names(opt_summary_df) |>
  gsub("^Time_single_opt_", "T_so_", x = _) |>
  gsub("^Time_opt_", "T_o_", x = _) |>
  gsub("^Success_single_opt_", "S_so_", x = _) |>
  gsub("^Success_opt_", "S_o_", x = _)

names(min_summary_df) <- names(min_summary_df) |>
  gsub("^Time_single_min_", "T_sm_", x = _) |>
  gsub("^Time_min_", "T_m_", x = _) |>
  gsub("^Success_single_min_", "S_sm_", x = _) |>
  gsub("^Success_min_", "S_m_", x = _)

# 輸出確認
print(opt_summary_df)
print(min_summary_df)