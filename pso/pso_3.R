---
  title: "Metaheuristic Optimization Algorithm"
author: "Ping-Yang Chen"
output: html_document
date: "2024-05-24"
runtime: shiny
---
  
  ```{r setup, include=FALSE}
library(shiny)
knitr::opts_chunk$set(echo = TRUE)
```


EGGHOLDER FUNCTION
For function details and reference information, see: http://www.sfu.ca/~ssurjano/
  
  ```{r egg, eval=TRUE, echo=TRUE}
egg <- function(xx) {
  ###########################
  # INPUT:
  #   xx = c(x1, x2)
  ###########################
  # Force save xx to global R variable 'pars'
  i <<- i + 1
  pars[[i]] <<- xx
  #
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- -(x2 + 47) * sin(sqrt(abs(x2 + x1/2 + 47)))
  term2 <- -x1 * sin(sqrt(abs(x1 - (x2 + 47))))
  
  y <- term1 + term2
  return(y)
}
```


In R, the `optim()` function provides a few built-in optimization algorithms. Let's try implementing L-BFGS-B algorithm, which is one of the popularly used quasi-Newton method. 

```{r qn, fig.align='center', fig.width=6, fig.height=6}
i <- 0
pars <- list() # Store X during L-BFGS-B update
res_qn <- optim(runif(2)*1024-512, egg, method = "L-BFGS-B")
pars_qn <- pars
```



Here is the quick development of the PSO algorithm. 

```{r psoalg, fig.align='center', fig.width=6, fig.height=6}
pso <- function(fun, lower, upper, nswarm = 40, niter = 100) {
  # Tuning Parameters (fixed in this version)
  c1 <- 2
  c2 <- 2
  vk <- 4
  w <- 0.9
  # Set matrices of the search domain 
  stopifnot(length(lower) == length(upper))
  dswarm <- length(lower)
  lowmat <- matrix(lower, nrow = nswarm, ncol = dswarm, byrow = TRUE)
  uppmat <- matrix(upper, nrow = nswarm, ncol = dswarm, byrow = TRUE)
  # Initialize Velocity
  vel <- matrix(0, nrow = nswarm, ncol = dswarm)
  vmax <- (uppmat - lowmat)/vk
  # Initialize Swarm and particles' function values
swarm <- matrix(runif(nswarm*dswarm), nrow = nswarm, ncol = dswarm)*(uppmat - lowmat) + lowmat
fswarm <- numeric(nSwarm)
for (i in 1:nswarm) {
  fswarm[i] <- fun(swarm[i,])
}
# Initialize Local Bests and Global Best
lbest <- swarm
flbest <- fswarm
gbest <- swarm[which.min(flbest),]
fgbest <- min(fswarm)
# Start PSO run
for (iter in 1:niter) {
  # Update Velocity
  gbmat <- matrix(gbest, nrow = nswarm, ncol = dswarm, byrow = TRUE)
  r1 <- matrix(runif(nswarm*dswarm), nrow = nswarm, ncol = dswarm)
  r2 <- matrix(runif(nswarm*dswarm), nrow = nswarm, ncol = dswarm)
  vel <- w*vel + c1*r1*(lbest - swarm) + c2*r2*(gbmat - swarm)
  # Refine Velocity not exceeding Vmax
  vlow <- which(vel < -vmax)
  vel[vlow] <- -vmax[vlow]
  vupp <- which(vel > vmax)
  vel[vupp] <- vmax[vupp]
  # Update Swarm
  swarm <- vel + swarm
  # Refine Swarm not exceeding boundary
  swlow <- which(swarm < lowmat)
  swarm[swlow] <- lowmat[swlow]
  swupp <- which(swarm > uppmat)
  swarm[swupp] <- uppmat[swupp]
  # Update particles' function values, local and global bests
  fswarm <- numeric(nSwarm)
  for (i in 1:nswarm) {
    fswarm[i] <- fun(swarm[i,])
    if (fswarm[i] < flbest[i]) {
      flbest[i] <- fswarm[i]
      lbest[i,] <- swarm[i,]
    }
  }
  gbest <- swarm[which.min(flbest),]
  fgbest <- min(fswarm)
}
return(list( par = gbest, val = fgbest ))
}

```


Let's run PSO with 64 particles for 100 iterations.
```{r pso, fig.align='center', fig.width=6, fig.height=6}
nSwarm <- 64
nIter <- 100
#
i <- 0
pars <- list() # Store X during PSO update
res_pso <- pso(egg, lower = c(-512, -512), upper = c(512, 512), nswarm = nSwarm, niter = nIter)
pars_pso <- pars
```


```{r psocollect, fig.align='center', fig.width=6, fig.height=6}
# Refine the PSO storage using list of length niter + 1
ct <- 1
pso_x <- lapply(1:(nIter+1), function(k) matrix(0, nSwarm, 2))
for (i in 1:(nIter+1)) {
  for (j in 1:nSwarm) {
    pso_x[[i]][j,] <- pars_pso[[ct]]
    ct <- ct + 1
  }
}
```


Observe the updating procedure of the L-BFGS-B algorithm in `optim()` and the PSO algorithm.


```{r plotoptim, eval = TRUE, echo = FALSE, fig.align='center', fig.width=6, fig.height=6}
drawOptim <- function(x, bg) {
  image(bg$xg, bg$xg, bg$eggmat, col = hcl.colors(20, "TealGrn", rev = TRUE),
        xlim = 520*c(-1, 1), ylim = 520*c(-1, 1), xlab = "x1", ylab = "x2")
  points(512, 404.2319, pch = 18, cex = 1.8, col = "red")
  #
  points(x[,1], x[,2], pch = 16, cex = 1, col = "black")
}
```


```{r qnshiny, eval = TRUE, echo = FALSE, fig.align='center', fig.width=6, fig.height=10}
fluidPage(
  fluidRow(
    column(6, 
           sliderInput("iter1", tags$h3("L-BFGS-B Iteration"), 1, length(pars_qn), 
                       value = 1, step = 1, width = 480),
           imageOutput("fig1")
           ),
    column(6, 
           sliderInput("iter2", tags$h3("PSO Iteration"), 0, 100, 
                       value = 0, step = 1, width = 480),
           imageOutput("fig2")
           )
  ))

eggmat <- matrix(0, 200, 200)
xg <- seq(-512, 512, length = 200)
for (ix in 1:200) {
  for (iy in 1:200) {
    eggmat[ix,iy] <- egg( c(xg[ix], xg[iy]) )
  }
}

rlist <- reactiveValues(bg = list(xg = xg, eggmat = eggmat))

output$fig1 <- renderImage({
 outfile <- tempfile(fileext='ex0.png')
 png(outfile, width=480, height=480)
 drawOptim(matrix(pars_qn[[input$iter1]], 1, 2), rlist$bg)
 dev.off()
 list(src = outfile, contentType = 'image/png', width = 480, height = 480, alt = "")
}, deleteFile = TRUE)
output$fig2 <- renderImage({
 outfile <- tempfile(fileext='ex1.png')
 png(outfile, width=480, height=480)
 drawOptim(pso_x[[input$iter2+1]], rlist$bg)
 dev.off()
 list(src = outfile, contentType = 'image/png', width = 480, height = 480, alt = "")
}, deleteFile = TRUE)
```




<!--
## D-optimal Design Simulation

```{r, fig.align='center', fig.width=6, fig.height=6}
set.seed(1)
n <- 20
s <- 0.5
a <- .5; b <- 5
xr <- runif(n, 0, 1)

ds <- c(0.25, 0.5, 1)
for (di in 1:3) {
  #di <- 3
  x <- xr*ds[di]
  e <- rnorm(n, 0, s)
  y <- a + b*x + e
  mdl <- lm(y ~ x)
  ttr <- summary(mdl)$coefficients
  
  rownames(ttr) <- c("a", "b")
  print(round(ttr, 4))
  
  cat("\n")
  plot(x, y, xlim = c(0, max(ds)), ylim = c(0, a+b*max(ds)), pch = 20)
  abline(a = a, b = b, col = 'red', lty = 2)
  title(paste0("Reg. Model Fitted w. Sampled Data from [0, ", ds[di], "]"))
}

n <- 30
xg <- list(c(0, 0.25, 0.5, 0.75, 1),
            c(0, 0.5, 1),
            c(0, 1))
for (i in 1:length(xg)) {
  x <- rep(xg[[i]], each = 30/length(xg[[i]]))  
  e <- rnorm(n, 0, s)
  y <- a + b*x + e
  mdl <- lm(y ~ x)
  ttr <- summary(mdl)$coefficients
  
  rownames(ttr) <- c("a", "b")
  print(round(ttr, 4))
  
  cat("\n")
  plot(x, y, xlim = c(0, max(ds)), ylim = c(0, a+b*max(ds)), pch = 20)
  abline(a = a, b = b, col = 'red', lty = 2)
  #title(paste0("Reg. Model Fitted w. ", length(xg[[i]]), "grids"))
}
```

## T-optimal Design Illustration

```{r, fig.align='center', fig.width=6, fig.height=6}
f0func <- function(x) { 0.5*x^2 + x*2 + 1 }
f1func <- function(x, p) { p[1] + p[2]*x }
obj <- function(x, design) {
  val <- 0
  for (i in 1:length(design)) {
    aa <- f0func(design[i])
    bb <- f1func(design[i], x)
    val <- val + (aa - bb)^2
  }
  return(val)
}

#
plotdisc <- function(xx, p2, val, frange, cSet) {
  lwdSet <- 2
  n <- 200
  xseq <- seq(-2, 2, length = n)
  f0 <- f0func(xseq)
  f1 <- p2[1] + p2[2]*xseq
  #
  plot(xseq, f0, type = "l", lwd = lwdSet, col = cSet[1], 
     ylim = frange, xlim = c(min(xseq), max(xseq)*1.5), xlab = "X", ylab = "eta")
  text(xseq[n], f0[n], "eta_t", col = cSet[1], pos = 4)
  #
  points(xseq, f1, type = "l", lwd = lwdSet, col = cSet[2])  
  text(xseq[n], f1[n], sprintf("eta_2"), col = cSet[2], pos = 4)
  #
  for (i in 1:length(xx)) {
    arrows(xx[i], frange[1] - 100, xx[i], max(f0func(xx[i]), f1func(xx[i], p2)), lty = 2, code = 0, col = "firebrick")
    text(xx[i], frange[1], sprintf("x%d", i), col = "firebrick", bg  = "white", cex = 1.1)
    arrows(xx[i], f0func(xx[i]), xx[i], f1func(xx[i], p2), code = 0, lwd = 4, lend = 2, col = "red")
  }
  text(xseq[1], frange[2], sprintf("Delta = %.4f", val), pos = 4, cex = 1.5)
}


#
n <- 200
xseq <- seq(-2, 2, length = n)
f0 <- f0func(xseq)
f1set <- matrix(0, nrow = length(xseq), ncol = 0)
####
xx <- c(-1, 0, 1)
p21 <- c(2, 3)
val21 <- obj(p21, xx)
f1set <- cbind(f1set, f1func(xseq, p21))
####
xx <- c(-1, 0, 1)
out <- optim(c(0, 0), obj, method = "L-BFGS-B", design = xx)
f1set <- cbind(f1set, f1func(xseq, out$par))
####
xx2 <- c(-2, 0, 2)
out2 <- optim(c(0, 0), obj, method = "L-BFGS-B", design = xx2)
f1set <- cbind(f1set, f1func(xseq, out2$par))
####

frange <- range(c(range(f1set), range(f0)))
cSet <- c("dodgerblue2", "darkorange1", "forestgreen", "firebrick")


plotdisc(xx = xx, p2 = p21, val = val21, frange = frange, cSet = c("dodgerblue2", "#666666"))
plotdisc(xx = xx, p2 = out$par, val = out$value, frange = frange, cSet = c("dodgerblue2", "darkorange1"))

plotdisc(xx = xx2, p2 = out2$par, val = out2$value, frange = frange, cSet = c("dodgerblue2", "darkorange1"))


```

-->
