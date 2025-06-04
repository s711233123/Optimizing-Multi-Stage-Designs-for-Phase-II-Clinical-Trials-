run_ga_single_stage <- function(seed, q_val, cliRequirement,
                                nStage,
                                maxIter = 600,
                                popSize = 256,
                                nMaxRange = c(10, 70)) {
  cat(sprintf(" Stage: maxIter = %d\n", maxIter))  
  set.seed(seed)
  nMinEachInterim <- rep(1, nStage)
  lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  fitness_function <- function(x) {
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim,
                               cliRequirement = cliRequirement, q = q_val)
    if (is.nan(val) || is.infinite(val)) return(-1e8)
    return(-val)
  }
  
  # ⏱ 手動計時
  start_time <- Sys.time()
  
  ga_res <- ga(
    type = "real-valued",
    fitness = fitness_function,
    lower = lower,
    upper = upper,
    popSize = popSize,
    maxiter = maxIter,
    pmutation = 0.1,
    pcrossover = 0.9,
    seed = seed,
    monitor = FALSE
  )
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  sol <- ga_res@solution[1, ]
  design <- kStageFreqCrit(
    nPolarized = sol[2:nStage],
    rProportion = sol[(nStage + 1):(2 * nStage)],
    nMax = sol[1],
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement
  )
  
  return(list(
    en = design$en,
    rseq = design$rseq,
    nseq = design$nseq,
    pet_seq = design$pet_seq,
    cputime = elapsed_time  # ✅ 用秒數表示
  ))
}
#----------------------------------------------------------------------------------
run_ga_two_stage <- function(seed, q_stage1, q_stage2, cliRequirement,
                             nStage,
                             maxIter1 = 100, maxIter2 = 500,
                             popSize = 256,
                             nMaxRange = c(10, 70),
                             upper2_is_fixed = FALSE) {
  cat(sprintf("Stage1: maxIter = %d, Stage2: maxIter = %d, popSize = %d\n", maxIter1, maxIter2, popSize))
  set.seed(seed)
  nMinEachInterim <- rep(1, nStage)
  lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  # === Stage 1: 搜尋 nMax ===
  fitness1 <- function(x) {
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim, cliRequirement = cliRequirement, q = q_stage1)
    if (is.nan(val) || is.infinite(val)) return(-1e8)
    return(-val)
  }
  
  start1 <- Sys.time()
  ga1 <- ga(
    type = "real-valued",
    fitness = fitness1,
    lower = lower,
    upper = upper,
    popSize = popSize,
    maxiter = maxIter1,
    pmutation = 0.1,
    pcrossover = 0.9,
    seed = seed,
    monitor = FALSE
  )
  end1 <- Sys.time()
  time1 <- as.numeric(difftime(end1, start1, units = "secs"))
  
  fixed_nMax <- round(ga1@solution[1, 1])
  
  # === Stage 2: 精調 nPolarized / rProportion ===
  refined_lower <- c(fixed_nMax, rep(0.0 * pi, nStage - 1), rep(0, nStage))
  refined_upper <- c(if (upper2_is_fixed) fixed_nMax else fixed_nMax + 20,
                     rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  fitness2 <- function(x) {
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim, cliRequirement = cliRequirement, q = q_stage2)
    if (is.nan(val) || is.infinite(val)) return(-1e8)
    return(-val)
  }
  
  start2 <- Sys.time()
  ga2 <- ga(
    type = "real-valued",
    fitness = fitness2,
    lower = refined_lower,
    upper = refined_upper,
    popSize = popSize,
    maxiter = maxIter2,
    pmutation = 0.1,
    pcrossover = 0.9,
    seed = seed,
    monitor = FALSE
  )
  end2 <- Sys.time()
  time2 <- as.numeric(difftime(end2, start2, units = "secs"))
  
  # === 整合結果 ===
  sol <- ga2@solution[1, ]
  design <- kStageFreqCrit(
    nPolarized = sol[2:nStage],
    rProportion = sol[(nStage + 1):(2 * nStage)],
    nMax = sol[1],
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement
  )
  
  return(list(
    en = design$en,
    rseq = design$rseq,
    nseq = design$nseq,
    pet_seq = design$pet_seq,
    cputime = time1 + time2  # ⏱ 兩段加總秒數
  ))
}
