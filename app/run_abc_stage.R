abc_opt <- function(fn, lower, upper, foodNumber = 256, maxIter = 100, limit = 100) {
  D <- length(lower)
  if (any(is.na(lower)) || any(is.na(upper))) {
    stop("❌ abc_opt() 偵測到 NA in lower/upper — 請檢查 refined_upper 是否有 NA")
  }
  
  Foods <- matrix(runif(foodNumber * D, lower, upper), nrow = foodNumber, byrow = TRUE)
  Fitness <- apply(Foods, 1, fn)
  
  trial <- rep(0, foodNumber)
  bestIndex <- suppressWarnings(which.min(Fitness))
  if (length(bestIndex) == 0) {
    stop("❌ abc_opt() 找不到有效最小 Fitness，which.min(Fitness) 為空")
  }
  bestFood <- Foods[bestIndex, ]
  bestFitness <- Fitness[bestIndex]
  
  for (iter in 1:maxIter) {
    for (i in 1:foodNumber) {
      param2change <- sample(1:D, 1)
      neighbour <- sample(setdiff(1:foodNumber, i), 1)
      phi <- runif(1, -1, 1)
      
      vi <- Foods[i, ]
      vi[param2change] <- Foods[i, param2change] + phi * (Foods[i, param2change] - Foods[neighbour, param2change])
      vi <- pmin(pmax(vi, lower), upper)
      
      newFit <- fn(vi)
      if (newFit < Fitness[i]) {
        Foods[i, ] <- vi
        Fitness[i] <- newFit
        trial[i] <- 0
      } else {
        trial[i] <- trial[i] + 1
      }
    }
    
    for (i in 1:foodNumber) {
      prob <- Fitness[i] / sum(Fitness)
      if (runif(1) > prob) {
        vi <- runif(D, lower, upper)
        newFit <- fn(vi)
        if (newFit < Fitness[i]) {
          Foods[i, ] <- vi
          Fitness[i] <- newFit
          trial[i] <- 0
        }
      }
    }
    
    for (i in 1:foodNumber) {
      if (trial[i] > limit) {
        Foods[i, ] <- runif(D, lower, upper)
        Fitness[i] <- fn(Foods[i, ])
        trial[i] <- 0
      }
    }
    
    currentBest <- which.min(Fitness)
    if (Fitness[currentBest] < bestFitness) {
      bestFitness <- Fitness[currentBest]
      bestFood <- Foods[currentBest, ]
    }
  }
  return(list(par = bestFood, value = bestFitness))
}

#----------------------------------------------------------------------------------
run_abc_single_stage <- function(seed, q_val, cliRequirement,
                                 nStage,
                                 maxIter = 100,
                                 foodNumber = 256,
                                 nMaxRange = c(10, 70)) {
  cat(sprintf(" Stage: maxIter = %d\n", maxIter))  
  set.seed(seed)
  nMinEachInterim <- rep(1, nStage)
  lower <- c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage))
  upper <- c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage))
  
  fitness_fn <- function(x) {
    nMax <- x[1]
    if (sum(nMinEachInterim) > floor(nMax)) return(1e8)
    val <- kStageCompromiseObj(x, nMin = nMinEachInterim, cliRequirement = cliRequirement, q = q_val)
    if (is.nan(val) || is.infinite(val)) return(1e8)
    return(val)
  }
  
  start_time <- Sys.time()
  res <- abc_opt(fitness_fn, lower, upper, foodNumber, maxIter, limit = 100)
  end_time <- Sys.time()
  
  sol <- res$par
  design <- kStageFreqCrit(
    nPolarized = sol[2:nStage],
    rProportion = sol[(nStage + 1):(2 * nStage)],
    nMax = sol[1],
    nMin = nMinEachInterim,
    cliRequirement = cliRequirement
  )
  
  attr(design, "fixed_nMax") <- sol[1]
  return(list(
    en = design$en,
    rseq = design$rseq,
    nseq = design$nseq,
    pet_seq = design$pet_seq,
    cputime = as.numeric(difftime(end_time, start_time, units = "secs"))
  ))
}

#----------------------------------------------------------------------------------
run_abc_two_stage <- function(seed, q_stage1, q_stage2, cliRequirement,
                                    nStage,
                                    maxIter1 = 100, maxIter2 = 500,
                                    foodNumber = 256,
                                    nMaxRange = c(10, 70),
                                    upper2_is_fixed = FALSE) {
  cat(sprintf("Stage1: maxIter = %d, Stage2: maxIter = %d, foodNumber = %d\n", maxIter1, maxIter2, foodNumber)) 
  set.seed(seed)
  nMinEachInterim <- rep(1, nStage)
  
  # === 第一階段搜尋 nMax ===
  attempt <- 1
  repeat {
    res1 <- tryCatch({
      abc_opt(
        fn = function(x) {
          nMax <- x[1]
          if (sum(nMinEachInterim) > floor(nMax)) return(1e8)
          val <- kStageCompromiseObj(x, nMin = nMinEachInterim,
                                     cliRequirement = cliRequirement, q = q_stage1)
          if (is.nan(val) || is.infinite(val)) return(1e8)
          return(val)
        },
        lower = c(nMaxRange[1], rep(0.0 * pi, nStage - 1), rep(0, nStage)),
        upper = c(nMaxRange[2], rep(0.5 * pi, nStage - 1), rep(1, nStage)),
        foodNumber = foodNumber,
        maxIter = maxIter1,
        limit = 100
      )
    }, error = function(e) return(NULL))
    
    if (!is.null(res1) && !is.null(res1$par) && !is.na(res1$par[1])) break
    attempt <- attempt + 1
    if (attempt > 10) return(NULL)
  }
  
  fixed_nMax <- round(res1$par[1])
  
  # === 第二階段固定 nMax 後搜尋設計 ===
  res2 <- abc_opt(
    fn = function(x) {
      nMax <- x[1]
      if (sum(nMinEachInterim) > floor(nMax)) return(1e8)
      val <- kStageCompromiseObj(x, nMin = nMinEachInterim,
                                 cliRequirement = cliRequirement, q = q_stage2)
      if (is.nan(val) || is.infinite(val)) return(1e8)
      return(val)
    },
    lower = c(fixed_nMax, rep(0.0 * pi, nStage - 1), rep(0, nStage)),
    upper = c(if (upper2_is_fixed) fixed_nMax else nMaxRange[2],
              rep(0.5 * pi, nStage - 1), rep(1, nStage)),
    foodNumber = foodNumber,
    maxIter = maxIter2,
    limit = 100
  )
  
  # === 建立設計結果 ===
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
    pet_seq = design$pet_seq
  ))
}
