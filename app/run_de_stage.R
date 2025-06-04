run_de_single_stage <- function(seed, q_val, cliRequirement,
                                nStage,
                                maxIter = 600,
                                NP = 256,
                                nMaxRange = c(10, 70)) {
  cat(sprintf("Stage: maxIter = %d, NP = %d\n", maxIter, NP))  
  set.seed(seed)
  nMinEachInterim <- rep(1, nStage)
  lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  fitness_fn <- function(x) {
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim, cliRequirement = cliRequirement, q = q_val)
    if (is.nan(val) || is.infinite(val)) return(1e8)
    return(val)
  }
  
  time_used <- system.time({
    res <- DEoptim(
      fn = fitness_fn,
      lower = lower, upper = upper,
      control = DEoptim.control(
        NP = NP,
        F = 0.7,
        CR = 0.6,
        itermax = maxIter,
        trace = FALSE,
        strategy = 6,
        parallelType = 2,
        reltol = 1e-6
      )
    )
  })
  
  sol <- res$optim$bestmem
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
    cputime = time_used["elapsed"]
  ))
}

#----------------------------------------------------------------------------------
run_de_two_stage <- function(seed, q_stage1, q_stage2, cliRequirement,
                             nStage,
                             maxIter1 = 100, maxIter2 = 500,
                             NP = 256,
                             nMaxRange = c(10, 70),
                             upper2_is_fixed = FALSE) {
  cat(sprintf("Stage1: maxIter = %d, Stage2: maxIter = %d, NP = %d\n", maxIter1, maxIter2, NP))  
  set.seed(seed)
  nMinEachInterim <- rep(1, nStage)
  lower1 <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper1 <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  fitness1 <- function(x) {
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim, cliRequirement = cliRequirement, q = q_stage1)
    if (is.nan(val) || is.infinite(val)) return(1e8)
    return(val)
  }
  
  time1 <- system.time({
    res1 <- DEoptim(
      fn = fitness1,
      lower = lower1, upper = upper1,
      control = DEoptim.control(
        NP = NP,
        F = 0.7,
        CR = 0.6,
        itermax = maxIter1,
        trace = FALSE,
        strategy = 6,
        parallelType = 2,
        reltol = 1e-6
      )
    )
  })
  
  fixed_nMax <- round(res1$optim$bestmem[1])
  refined_lower <- c(fixed_nMax, rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper2_nMax <- if (upper2_is_fixed) fixed_nMax else nMaxRange[2]
  refined_upper <- c(upper2_nMax, rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  fitness2 <- function(x) {
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim, cliRequirement = cliRequirement, q = q_stage2)
    if (is.nan(val) || is.infinite(val)) return(1e8)
    return(val)
  }
  
  time2 <- system.time({
    res2 <- DEoptim(
      fn = fitness2,
      lower = refined_lower, upper = refined_upper,
      control = DEoptim.control(
        NP = NP,
        F = 0.7,
        CR = 0.6,
        itermax = maxIter2,
        trace = FALSE,
        strategy = 6,
        parallelType = 2,
        reltol = 1e-6
      )
    )
  })
  
  sol <- res2$optim$bestmem
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
    cputime = time1["elapsed"] + time2["elapsed"]
  ))
}
