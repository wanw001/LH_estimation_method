
setwd("C:/Users/wanwank/OneDrive - University of Tasmania/PhD_Wanwan-Kurniawan/Code/Ch1/Github/LH_estimation_method")

library(rstan)

## Import the length composition sample data
#Comp

## Format the data for rstan
LOW <- Comp$Ll
MID <- Comp$L
UP <- Comp$Lu
nl <- length(LOW)
Lmin <- LOW[1]
Lmax <- UP[nl]
N <- sum(Comp$Nl.obs)
Nl <- unlist(Comp$Nl.obs)
Lopt <- UP[which(Nl==max(Nl))][1]
X <- 200
sa.e <- 0.01
Md <- -log(sa.e)/(X-1)
CV <- 0.1
Stan.file <- c("lbspr-dirichlet.stan")

stan_dat <- list(
  nl = nl,           
  LOW = LOW, 
  MID = MID,
  UP = UP,
  Lmin = Lmin,
  Lmax = Lmax,
  Lopt = Lopt,
  N = N,
  Nl = Nl,
  X = X,
  sa = sa.e,
  Md = Md,
  CV = CV)

chains <- 1
warmup <- 500
iter <- 3000

## Fitting using rstan
fit <- stan(
  file = stf,         # Stan program
  data = stan_dat,    # named list of data
  chains = chains,    # number of Markov chains
  warmup = warmup,    # number of warmup iterations per chain
  iter = iter,        # total number of iterations per chain
)

## Check traceplot and show summary of results
traceplot(fit, inc_warmup = TRUE)
fit

## Store the posterior estimates
post_est <- as.matrix(fit)

