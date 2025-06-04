library(Rcpp)
library(RcppArmadillo)
sourceCpp('run/kstagecli/kStageP2A.cpp')
### ---------------------------------------------------------------------------
### The objective function of the K-stage compromise design
### ---------------------------------------------------------------------------
kStageCompromiseObj <- function(particle, nMin, cliRequirement = NULL, q = 0.5) {
  if (is.null(cliRequirement)) {
    p0 <- 0.1; p1 <- 0.3; t1eThres <- 0.1; t2eThres <- 0.1
  } else {
    p0 <- cliRequirement$p0
    p1 <- cliRequirement$p1
    t1eThres <- cliRequirement$alpha
    t2eThres <- cliRequirement$beta
  }
  
  nStage <- as.integer(length(particle)/2)
  nMax <- particle[1]
  nPolarized <- particle[2:nStage]
  rProportion <- particle[(nStage + 1):length(particle)]
  
  result <- kStageFreqCrit(nPolarized, rProportion, nMax, nMin, cliRequirement)
  if ((result$t1e <= t1eThres) & (result$t2e <= t2eThres)) {
    loss <- q * (result$en/nMax + nMax) + (1 - q) * result$en
    return(loss)
  } else {
    return(1e8)
  }
}


### ---------------------------------------------------------------------------
### The function computing the expected sample size under null hypothesis for 
###  the input PSO particle 
### ---------------------------------------------------------------------------
kStageFreqCrit <- function(nPolarized, rProportion, nMax, nMin, cliRequirement = NULL) {
  
  if (is.null(cliRequirement)) {
    p0 <- 0.1
    p1 <- 0.3
    t1eThres <- 0.1
    t2eThres <- 0.1
  } else {
    p0 <- cliRequirement$p0
    p1 <- cliRequirement$p1
    t1eThres <- cliRequirement$alpha
    t2eThres <- cliRequirement$beta
  }
  
  nStage <- length(nPolarized) + 1
  nseq <- get_cohort(nPolarized, nStage, nMax, nMin)
  rseq <- get_cutoff(rProportion, nseq)
  #
  return(kStageP2A_Cpp(p0, p1, nseq, rseq))
}