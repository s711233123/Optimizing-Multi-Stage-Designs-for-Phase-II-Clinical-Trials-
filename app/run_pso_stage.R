run_pso_single_stage <- function(seed, q_val, cliRequirement,
                                 nStage,
                                 maxIter = 600,
                                 nSwarm = 256,
                                 nMaxRange = c(10, 70)) {
  cat(sprintf(" Stage: maxIter = %d\n", maxIter))  
  nMinEachInterim <- rep(1, nStage)
  lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "basic")
  
  res <- globpso(
    objFunc = kStageCompromiseObj,
    PSO_INFO = alg_setting,
    lower = lower, upper = upper,
    seed = seed, verbose = FALSE,
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement,
    q = q_val
  )
  
  design <- kStageFreqCrit(
    nPolarized = res$par[2:nStage],
    rProportion = res$par[(nStage + 1):(2 * nStage)],
    nMax = res$par[1],
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement
  )
  
  return(list(
    en = design$en,
    rseq = design$rseq,
    nseq = design$nseq,
    pet_seq = design$pet_seq,
    cputime = res$cputime
  ))
}

#----------------------------------------------------------------------------------
run_pso_two_stage <- function(seed, q_stage1, q_stage2, cliRequirement,
                              nStage,
                              maxIter1 = 100, maxIter2 = 500,
                              nSwarm = 256,
                              nMaxRange = c(10, 70),
                              upper2_is_fixed = FALSE) {
  cat(sprintf("Stage1: maxIter = %d, Stage2: maxIter = %d, nSwarm = %d\n", maxIter1, maxIter2, nSwarm))
  nMinEachInterim <- rep(1, nStage)
  lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  # Stage 1
  alg1 <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter1, psoType = "basic")
  res1 <- globpso(
    objFunc = kStageCompromiseObj,
    PSO_INFO = alg1,
    lower = lower, upper = upper,
    seed = seed, verbose = FALSE,
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement,
    q = q_stage1
  )
  fixed_nMax <- round(res1$par[1])
  
  # Stage 2
  refined_lower <- c(fixed_nMax, rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper2_nMax <- if (upper2_is_fixed) fixed_nMax else fixed_nMax + 20
  refined_upper <- c(upper2_nMax, rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  alg2 <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter2, psoType = "basic")
  res2 <- globpso(
    objFunc = kStageCompromiseObj,
    PSO_INFO = alg2,
    lower = refined_lower, upper = refined_upper,
    seed = seed, verbose = FALSE,
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement,
    q = q_stage2
  )
  
  design <- kStageFreqCrit(
    nPolarized = res2$par[2:nStage],
    rProportion = res2$par[(nStage + 1):(2 * nStage)],
    nMax = res2$par[1],
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement
  )
  
  return(list(
    en = design$en,
    rseq = design$rseq,
    nseq = design$nseq,
    pet_seq = design$pet_seq,
    cputime = res1$cputime + res2$cputime
  ))
}


