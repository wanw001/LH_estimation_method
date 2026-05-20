#### EM using LH estimation method ####

setwd("C:/Users/wanwank/OneDrive - University of Tasmania/PhD_Wanwan-Kurniawan/Code/Ch1/Github/LH_estimation_method")

library(rstan)

### Import the length composition sample data
Comp <- read.csv("Comp.csv")

### Format the data for rstan

LOW <- Comp$Ll #vector of lower bounds of size bins
MID <- Comp$L #vector of mid size of bins
UP <- Comp$Lu #vector of upper bounds of size bins
nl <- length(LOW) #number of size bins
Lmin <- LOW[1] #lower bound of the lowest bin
Lmax <- UP[nl] #upper bound of the highest bin
N <- sum(Comp$Nl.obs) #total abundance
Nl <- unlist(Comp$Nl.obs) #vector of observed number at bin
Lopt <- UP[which(Nl==max(Nl))][1] #mid value of bin with highest abundance
X <- 200 #number of relative age classes
sa.e <- 0.01 #fraction of individuals surviving at maximum age
Md <- -log(sa.e)/(X-1) #natural mortality rate per relative age
CV <- 0.1 #coefficient of variation of length at age

stan_file <- "lh_est_code.stan"

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

### Fitting using rstan

fit <- stan(
  file = stan_file,         # Stan program
  data = stan_dat,    # named list of data
  chains = chains,    # number of Markov chains
  warmup = warmup,    # number of warmup iterations per chain
  iter = iter,        # total number of iterations per chain
)

### Check traceplot and show summary of results
traceplot(fit, inc_warmup = TRUE)
fit

### Store the posterior estimates
post_par <- as.matrix(fit)

### Compute median and 90% CrI
MED <- CrI <- NULL
np <- 5
for (p in 1:5){
  name <- colnames(as.data.frame(post_par))[p]
  MED <- apply(post_par[, 1:np], 2, median)
  CrI <- apply(post_par[, 1:np], 2, quantile, probs = c(0.05, 0.95))
}
