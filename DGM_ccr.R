#### DGM with Continuous and constant recruitment ####

setwd("C:/Users/wanwank/OneDrive - University of Tasmania/PhD_Wanwan-Kurniawan/Code/Ch1/Github/LH_estimation_method")

library(ggplot2)
library(arrow)

### Generate length distribution in fine binning  
sizeProp.ccr <- function(linf, mk, l50, l95r, l0r, sa, cv, X, w){
  #linf = Linf
  #mk = M/K
  #l50 = L50
  #l95r = L95/L50
  #l0r = length at zero age
  #sa = fraction of individuals surviving at maximum age
  #cv = CV of length at age
  #X = number of relative age class
  #w = fine bin width representing population length
  
  ## Construct relative abundance at age (relative with respect to recruitment)
  m <- 0.3
  k <- m/mk
  amax <- -log(sa)/m
  md <- -log(sa)/(X-1)
  a0 <- 1/k*log(1-l0r)
  Ax <- (0:(X-1))/(X-1)*amax
  Lx <- linf*(1-exp(-k*(Ax-a0)))

  ## Determine size bins
  lmin <- 0
  lmax0 <- round(qnorm(0.999, max(Lx), cv*max(Lx)))
  nl <- ceiling((lmax0-lmin)/w)
  lmax <- lmin + nl*w
  Ll <- seq(lmin, lmax-w, by=w)
  Lu <- seq(lmin+w, lmax, by=w)
  L <- (Ll+Lu)/2
  
  R <- 1
  Surv <- 1
  for (x in 2:X){
    Surv[x] <- Surv[x-1]*exp(-md)
    #if (x==X){Surv[x] <- Surv[x]/(1-exp(-M))}
  }
  Nx <- R*Surv
  
  ## Converting relative abundance at age to proportion at size bin
  PP <- matrix(rep(NA, nl*X), nrow=nl, ncol=X)
  for (x in 1:X){
    la <- Lx[x]
    sd <- cv*la
    lu <- lmin + w
    PP[1,x] <- pnorm((lu-la)/sd)
    for (i in 2:(nl-1)){
      ll <- lmin + (i-1)*w
      lu <- lmin + i*w
      PP[i,x] <- pnorm((lu-la)/sd) - pnorm((ll-la)/sd)
    }
    ll <- lmax - w
    PP[nl,x] <- 1-pnorm((ll-la)/sd)
    
    PP[,x] <- PP[,x]/sum(PP[,x])
  }
  
  Nl0 <- unlist((matrix(Nx, nrow=1)%*%t(PP))[1,])
  pl <- Nl0/sum(Nl0)
  
  ## Computing probability of observed proportion of fish at size bins
  l95 <- l95r*l50
  alp <- log(19)/(l95-l50)
  Sl <- 1/(1 + exp(-alp*(L-l50))) #Observation selectivity
  Nl.o0 <- Sl*Nl0
  pl.o <- Nl.o0/sum(Nl.o0)
  
  df <- data.frame(L=L, Ll=Ll, Lu=Lu, pl=pl, pl.obs=pl.o)
  return(df)
}

### Generate length composition sample for a given sample size
sizeComp <- function(r, Prop, ns){
  
  #r = index of replicates
  #Prop = a dataframe or vector containing vector of length binning and proportion at length
  #ns = sample size
  
  set.seed(r)
  
  L <- Prop$L; Ll=Prop$Ll; Lu=Prop$Lu; nl <- dim(Prop)[1]
  pl <- Prop$pl; pl.o <- Prop$pl.obs
  
  ## Generating actual (random) abundance of fish at size bin from the proportion at size
  # #Directly draw the number at size from multinomial distribution
  # Nl <- Nl.o <- rep(0, nl)
  # A <- rmultinom(1, ns, pl); B <- rmultinom(1, ns, pl.o)
  # for (j in 1:nl){
  #   Nl[j] <- A[j,1]; Nl.o[j] <- B[j,1]
  # }
  
  #Draw size using function sample
  B1 <- sample(L, ns, replace=TRUE, pl.o); Nl.o <- rep(0, nl)
  for (i in 1:nl){
    B2 <- which(B1==L[i]); Nl.o[i] <- Nl.o[i] + length(B2)
  }
  
  df2 <- data.frame(L=L, Nl.obs=Nl.o, nl=nl, N.obs=sum(Nl.o), Ll=Ll, Lu=Lu, pl=pl, pl.obs=pl.o)
  
  return (df2)
}

### Rebin the size composition sample into observational binning
sizeComp.rebin <- function(B, Comp){
  
  #B = vector of bounds of observational bins from lower bound of lowest bin to upper bound of highest bin
  #Comp = data frame or vector containing vector of length binning and number at bin for population length binning
  
  nl <- length(B) - 1
  Ll <- B[1:nl]; Lu <- B[2:(nl+1)]; L <- (Ll + Lu)/2
  
  Nl.o <- NULL 
  for (i in 1:nl){
    q4 <- 0 
    J <- which(Comp$Ll>=Ll[i] & Comp$Ll<=Lu[i])
    lJ <- length(J)
    if (lJ>0){
      if (Ll[i]>Comp$Ll[1] & Ll[i]<Comp$Ll[J[1]]){
        pr1 <- (Comp$Ll[J[1]] - Ll[i]) / (Comp$Ll[J[1]] - Comp$Ll[J[1]-1])
        q4 <- q4 + pr1*Comp$Nl.obs[J[1]-1]
      }
      if (lJ>1){
        for (j in J[-lJ]){
          q4 <- q4 + Comp$Nl.obs[j] 
        }
      }
      if (Lu[i]<Comp$Lu[dim(Comp)[1]]){
        pr2 <- (Lu[i] - Comp$Ll[J[lJ]]) / (Comp$Lu[J[lJ]] - Comp$Ll[J[lJ]])
        q4 <- q4 + pr2*Comp$Nl.obs[J[lJ]]
      }else{
        q4 <- q4 + Comp$Nl.obs[J[lJ]]
      }
    }else if(Ll[i]>Comp$Ll[1] & Lu[i]<Comp$Lu[dim(Comp)[1]]){
      j <- max(which(Comp$Ll<Ll[i]))
      pr1 <- (Lu[i]-Ll[i]) / (Comp$Lu[j] - Comp$Ll[j])
      q4 <- q4 + pr1*Comp$Nl.obs[j]
    }else if (Ll[i]<Comp$Lu[dim(Comp)[1]] & Lu[i]>Comp$Lu[dim(Comp)[1]]){
      pr1 <- (Comp$Lu[dim(Comp)[1]] - Ll[nl]) / (Comp$Lu[dim(Comp)[1]] - Comp$Lu[dim(Comp)[1]-1])
      q4 <- q4 + pr1*Comp$Nl.obs[dim(Comp)[1]]
    }
    Nl.o <- c(Nl.o, round(q4))
  }
  df <- data.frame(L=L, Nl.obs=Nl.o, nl=rep(nl, nl), N.obs=rep(sum(Nl.o), nl), Ll=Ll, Lu=Lu)
  
  return(df)
}

