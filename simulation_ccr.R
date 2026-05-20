#### Simulating the length data under continulous recruitment ####

setwd("C:/Users/wanwank/OneDrive - University of Tasmania/PhD_Wanwan-Kurniawan/Code/Ch1/Github/LH_estimation_method")

library(ggplot2)
library(arrow)

source("DGM_ccr.R")

#### 2. New way of sensitivity analysis 

#### Parameters
LINF0 <- c(10, 20, 30, 40, 55, 70, 90) 
MK <- c(0.1, 0.5, 0.8, 1, 1.5, 2.2, 2.3, 3.0, 3.5)
L50 <- c(2, 5, 10, 20)
L95r <- c(1.2, 1.5, 2)
Sel <- list(c(2,2), c(1,1), c(1,2), c(1,3), 
            c(2,1), c(2,3), 
            c(3,1), c(3,2), c(3,3),
            c(4,1), c(4,2), c(4,3))
ISel <- 1:(length(Sel))
L0r <- c(0, 0.2)
SA <- c(0.001, 0.02) #Proxy for amax: 0.01-0.001, amax corresponds to 0.01-0.001 of surviving cohort
CV <- c(0.05, 0.2)
PR <- c(1, 0) #Recruitment periodicity: 1 continuous, 2. annual
pr <- PR[1]
#nl <- 100  #Original number of size bins (unbinned sizes)
w <- 0.5
X <-  1000 #Number of (relative) age class
NS <- c(200, 500, 1000, 5000, 20000, 100000) #Sample size
Red <- 30 #Replicates for discrete RV
Rec <- 100 #Replicates for continuous RV
WRB <- c(1, 2, 4, 8, 12, 0) #Bin width (Binned)
B.rls <- c(1.25, 3.75, 6.25, 8.75, 11.25, 13.75, 17.5, 22.5, 27.5, 32.5, 37.5, 45, 
           56.25, 68.75, 81.25, 118.75, 131.25)


## Sensitivity analysis for a0 - with replicates
pr <- 1
cv <- 0.1
L0r <- seq(0, 0.2, by=0.05)
sa <- 0.01
WRB2 <- WRB[c(2, 6)]
NS2 <- NS[4:6]
linf <- LINF0[4]
l50 <- L50[2]
l95r <- L95r[2]
Re <- 10

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (l0r in L0r){
  for (mk in MK){
    scp <- scp + 1
    Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, w=w)
    nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
    Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                     s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                     MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                     pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
    PROP <- rbind(PROP, Q1)
    Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                       L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
    PAR.p <- rbind(PAR.p, Par1)
    for (r in 1:Re){
      for (ns in NS2){
        if (r==1){scu <- scu + 1}else{
          Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
          scu <- Comp.ub2$Scen.u[1]
        }
        Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
        Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                         CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                         L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                         Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                         Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                         pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
        COMP.ub <- rbind(COMP.ub, Q2)
        Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                           nS=ns, Scen.p=scu)
        PAR.ub <- rbind(PAR.ub, Par2)
        for (wr in WRB2){
          if (wr>0){
            lmin <- Comp.ub$Ll[1]
            lmax <- Comp.ub$Lu[nl]
            nlr <- ceiling((lmax - lmin)/wr)
            B <- lmin + (0:nlr)*wr
          }else{
            ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
            B <- B.rls[1:ima]
          }
          if (r==1){scr <- scr + 1}else{
            Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
            scr <- Comp.rb2$Scen.r[1]
          }
          Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
          nlr <- dim(Comp.rb)[1]
          Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                           CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                           L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                           wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                           MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                           Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                           nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
          COMP.rb <- rbind(COMP.rb, Q3)
          Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                             nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                             wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
          PAR.rb <- rbind(PAR.rb, Par3)
        }
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-a0b-order_rev01.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-a0b-order_rev01.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-a0b-order_rev01.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-a0b-order_rev01.csv", row.names=FALSE)

## Sensitivity analysis for amax
pr <- 1
cv <- 0.1
l0r <- 0
SA <- seq(0.001, 0.05, by=0.001)
WRB2 <- WRB[c(2, 6)]
NS2 <- NS[4:5]
linf <- LINF[1]
l50 <- L50[2]
l95r <- L95r[2]
Re <- 10

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (sa in SA){
  for (mk in MK){
    scp <- scp + 1
    Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, nl=nl)
    nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
    Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                     s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                     MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                     pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
    PROP <- rbind(PROP, Q1)
    Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                       L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
    PAR.p <- rbind(PAR.p, Par1)
    Re <- 1
    for (r in 1:Re){
      for (ns in NS2){
        if (r==1){scu <- scu + 1}else{
          Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
          scu <- Comp.ub2$Scen.u[1]
        }
        Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
        Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                         CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                         L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                         Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                         Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                         pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
        COMP.ub <- rbind(COMP.ub, Q2)
        Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                           nS=ns, Scen.p=scu)
        PAR.ub <- rbind(PAR.ub, Par2)
        for (wr in WRB2){
          if (wr>0){
            lmin <- Comp.ub$Ll[1]
            lmax <- Comp.ub$Lu[nl]
            nlr <- ceiling((lmax - lmin)/wr)
            B <- lmin + (0:nlr)*wr
          }else{
            ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
            B <- B.rls[1:ima]
          }
          if (r==1){scr <- scr + 1}else{
            Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
            scr <- Comp.rb2$Scen.r[1]
          }
          Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
          nlr <- dim(Comp.rb)[1]
          Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                           CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                           L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                           wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                           MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                           Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                           nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
          COMP.rb <- rbind(COMP.rb, Q3)
          Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                             nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                             wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
          PAR.rb <- rbind(PAR.rb, Par3)
        }
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-amax-order.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-amax-order.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-amax-order.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-amax-order.csv", row.names=FALSE)

## Sensitivity analysis for amax - with replicates
pr <- 1
cv <- 0.1
l0r <- 0
SA <- c(0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05)
WRB2 <- WRB[c(2, 6)]
NS2 <- NS[4:6]
linf <- LINF[1]
l50 <- L50[2]
l95r <- L95r[2]
Re <- 10

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (sa in SA){
  for (mk in MK){
    scp <- scp + 1
    Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, nl=nl)
    nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
    Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                     s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                     MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                     pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
    PROP <- rbind(PROP, Q1)
    Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                       L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
    PAR.p <- rbind(PAR.p, Par1)
    for (r in 1:Re){
      for (ns in NS2){
        if (r==1){scu <- scu + 1}else{
          Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
          scu <- Comp.ub2$Scen.u[1]
        }
        Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
        Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                         CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                         L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                         Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                         Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                         pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
        COMP.ub <- rbind(COMP.ub, Q2)
        Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                           nS=ns, Scen.p=scu)
        PAR.ub <- rbind(PAR.ub, Par2)
        for (wr in WRB2){
          if (wr>0){
            lmin <- Comp.ub$Ll[1]
            lmax <- Comp.ub$Lu[nl]
            nlr <- ceiling((lmax - lmin)/wr)
            B <- lmin + (0:nlr)*wr
          }else{
            ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
            B <- B.rls[1:ima]
          }
          if (r==1){scr <- scr + 1}else{
            Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
            scr <- Comp.rb2$Scen.r[1]
          }
          Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
          nlr <- dim(Comp.rb)[1]
          Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                           CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                           L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                           wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                           MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                           Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                           nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
          COMP.rb <- rbind(COMP.rb, Q3)
          Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                             nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                             wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
          PAR.rb <- rbind(PAR.rb, Par3)
        }
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-amaxb-order.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-amaxb-order.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-amaxb-order.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-amaxb-order.csv", row.names=FALSE)

## Sensitivity analysis for CV - with replicates
pr <- 1
CV <- c(0.05, 0.1, 0.15, 0.2)
l0r <- 0
sa <- 0.01
WRB2 <- WRB[c(2, 6)]
NS2 <- NS[4:6]
linf <- LINF[1]
l50 <- L50[2]
l95r <- L95r[2]
Re <- 10

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (cv in CV){
  for (mk in MK){
    scp <- scp + 1
    Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, nl=nl)
    nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
    Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                     s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                     MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                     pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
    PROP <- rbind(PROP, Q1)
    Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                       L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
    PAR.p <- rbind(PAR.p, Par1)
    for (r in 1:Re){
      for (ns in NS2){
        if (r==1){scu <- scu + 1}else{
          Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
          scu <- Comp.ub2$Scen.u[1]
        }
        Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
        Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                         CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                         L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                         Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                         Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                         pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
        COMP.ub <- rbind(COMP.ub, Q2)
        Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                           nS=ns, Scen.p=scu)
        PAR.ub <- rbind(PAR.ub, Par2)
        for (wr in WRB2){
          if (wr>0){
            lmin <- Comp.ub$Ll[1]
            lmax <- Comp.ub$Lu[nl]
            nlr <- ceiling((lmax - lmin)/wr)
            B <- lmin + (0:nlr)*wr
          }else{
            ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
            B <- B.rls[1:ima]
          }
          if (r==1){scr <- scr + 1}else{
            Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
            scr <- Comp.rb2$Scen.r[1]
          }
          Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
          nlr <- dim(Comp.rb)[1]
          Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                           CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                           L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                           wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                           MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                           Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                           nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
          COMP.rb <- rbind(COMP.rb, Q3)
          Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                             nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                             wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
          PAR.rb <- rbind(PAR.rb, Par3)
        }
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-cvb-order.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-cvb-order.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-cvb-order.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-cvb-order.csv", row.names=FALSE)


## Sensitivity analysis for Linf - with replicates
pr <- 1
cv <- 0.1
l0r <- 0
sa <- 0.01
WRB2 <- WRB[c(2, 6)]
NS2 <- NS[3:6]
LINF <- NULL
for (linf0 in LINF0){
  low <- 0.9*linf0; up <- 1.1*linf0; inc <- (up - low)/4
  linf <- seq(low, up, inc)
  LINF <- c(LINF, linf)
}
# LINF <- c(9.5, 10, 11,
#           19, 20, 22,
#           32, 34, 36,
#           58, 61, 65,
#           89, 92, 95)
l50 <- L50[2]
l95r <- L95r[2]
Re <- 10

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (linf in LINF){
  for (mk in MK){
    scp <- scp + 1
    Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, w=w)
    nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
    Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                     s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                     MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                     pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
    PROP <- rbind(PROP, Q1)
    Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                       L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
    PAR.p <- rbind(PAR.p, Par1)
    for (r in 1:Re){
      for (ns in NS2){
        if (r==1){scu <- scu + 1}else{
          Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
          scu <- Comp.ub2$Scen.u[1]
        }
        Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
        Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                         CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                         L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                         Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                         Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                         pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
        COMP.ub <- rbind(COMP.ub, Q2)
        Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                           nS=ns, Scen.p=scu)
        PAR.ub <- rbind(PAR.ub, Par2)
        for (wr in WRB2){
          if (wr>0){
            lmin <- Comp.ub$Ll[1]
            lmax <- Comp.ub$Lu[nl]
            nlr <- ceiling((lmax - lmin)/wr)
            B <- lmin + (0:nlr)*wr
          }else{
            ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
            B <- B.rls[1:ima]
          }
          if (r==1){scr <- scr + 1}else{
            Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
            scr <- Comp.rb2$Scen.r[1]
          }
          Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
          nlr <- dim(Comp.rb)[1]
          Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                           CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                           L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                           wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                           MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                           Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                           nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
          COMP.rb <- rbind(COMP.rb, Q3)
          Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                             nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                             wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
          PAR.rb <- rbind(PAR.rb, Par3)
        }
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-linfb-order_rev01.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-linfb-order_rev01.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-linfb-order_rev01.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-linfb-order_rev01.csv", row.names=FALSE)

## Sensitivity analysis for sample size - with replicates
l0r <- 0
sa <- 0.01
cv <- 0.1
WRB2 <- WRB[c(2, 3, 6)]
NS2 <- NS
linf <- LINF0[4]
l50 <- L50[2]
l95r <- L95r[2]

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (mk in MK){
  scp <- scp + 1
  Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, nl=nl)
  nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
  Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                   s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                   MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                   pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
  PROP <- rbind(PROP, Q1)
  Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                     L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
  PAR.p <- rbind(PAR.p, Par1)
  for (r in 1:Red){
    for (ns in NS2){
      if (r==1){scu <- scu + 1}else{
        Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
        scu <- Comp.ub2$Scen.u[1]
      }
      Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
      Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                       CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                       L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                       Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                       Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                       pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
      COMP.ub <- rbind(COMP.ub, Q2)
      Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                         L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                         nS=ns, Scen.p=scu)
      PAR.ub <- rbind(PAR.ub, Par2)
      for (wr in WRB2){
        if (wr>0){
          lmin <- Comp.ub$Ll[1]
          lmax <- Comp.ub$Lu[nl]
          nlr <- ceiling((lmax - lmin)/wr)
          B <- lmin + (0:nlr)*wr
        }else{
          ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
          B <- B.rls[1:ima]
        }
        if (r==1){scr <- scr + 1}else{
          Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
          scr <- Comp.rb2$Scen.r[1]
        }
        Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
        nlr <- dim(Comp.rb)[1]
        Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                         CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                         L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                         wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                         MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                         Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                         nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
        COMP.rb <- rbind(COMP.rb, Q3)
        Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                           wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
        PAR.rb <- rbind(PAR.rb, Par3)
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-ns-order.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-ns-order.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-ns-order.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-ns-order.csv", row.names=FALSE)

## Sensitivity analysis for bin width - with replicates
l0r <- 0
sa <- 0.01
cv <- 0.1
WRB2 <- WRB
NS2 <- NS[4:6]
linf <- LINF0[4]
l50 <- L50[2]
l95r <- L95r[2]

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (mk in MK){
  scp <- scp + 1
  Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, nl=nl)
  nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
  Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                   s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                   MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                   pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
  PROP <- rbind(PROP, Q1)
  Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                     L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
  PAR.p <- rbind(PAR.p, Par1)
  for (r in 1:Red){
    for (ns in NS2){
      if (r==1){scu <- scu + 1}else{
        Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
        scu <- Comp.ub2$Scen.u[1]
      }
      Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
      Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                       CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                       L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                       Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                       Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                       pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
      COMP.ub <- rbind(COMP.ub, Q2)
      Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                         L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                         nS=ns, Scen.p=scu)
      PAR.ub <- rbind(PAR.ub, Par2)
      for (wr in WRB2){
        if (wr>0){
          lmin <- Comp.ub$Ll[1]
          lmax <- Comp.ub$Lu[nl]
          nlr <- ceiling((lmax - lmin)/wr)
          B <- lmin + (0:nlr)*wr
        }else{
          ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
          B <- B.rls[1:ima]
        }
        if (r==1){scr <- scr + 1}else{
          Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
          scr <- Comp.rb2$Scen.r[1]
        }
        Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
        nlr <- dim(Comp.rb)[1]
        Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                         CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                         L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                         wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                         MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                         Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                         nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
        COMP.rb <- rbind(COMP.rb, Q3)
        Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                           wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
        PAR.rb <- rbind(PAR.rb, Par3)
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-wr-order.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-wr-order.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-wr-order.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-wr-order.csv", row.names=FALSE)

## Sensitivity analysis for selectivity parameters - with replicates
cv <- 0.1
l0r <- 0
sa <- 0.01
WRB2 <- WRB[c(2, 3, 6)]
NS2 <- NS[4:6]
linf <- LINF0[4]
Sel <- list(c(1,1), c(1,2), c(1,3),
            c(2,1), c(2,2), c(2,3),
            c(3,1), c(3,2), c(3,3),
            c(4,1), c(4,2), c(4,3))
#ISel <- 1:(length(Sel))

PROP <- COMP.ub <- COMP.rb <- NULL
PAR.p <- PAR.ub <- PAR.rb <- NULL
scp <- scu <- scr <- 0
for (s in Sel){
  l50 <- L50[s[1]]; l95r <- L95r[s[2]]
  for (mk in MK){
    scp <- scp + 1
    Prop <- sizeProp.ccr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, cv=cv, X=X, nl=nl)
    nl <- dim(Prop)[1]; w <- Prop$Lu[1] - Prop$Ll[1]
    Q1 <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), CV=rep(cv, nl), L0r=rep(l0r, nl),
                     s.amax=rep(sa, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), Linf=rep(linf, nl),
                     MK=rep(mk, nl),  L=Prop$L, Ll=Prop$Ll, Lu=Prop$Lu, pl=Prop$pl,
                     pl.obs=Prop$pl.obs, nl=rep(nl, nl), wl=rep(w, nl), X=rep(X,nl))
    PROP <- rbind(PROP, Q1)
    Par1 <- data.frame(Scen.p=scp, perR=pr, CV=cv, L0r=l0r, s.amax=sa, L50=l50,
                       L95r=l95r, Linf=linf, MK=mk, nl=nl, wl=w, X=X)
    PAR.p <- rbind(PAR.p, Par1)
    for (r in 1:Red){
      for (ns in NS2){
        if (r==1){scu <- scu + 1}else{
          Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
          scu <- Comp.ub2$Scen.u[1]
        }
        Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
        Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                         CV=rep(cv, nl), L0r=rep(l0r, nl), s.amax=rep(sa, nl), L50=rep(l50, nl),
                         L95r=rep(l95r, nl), nl=rep(nl, nl), wl=rep(w, nl), X=rep(X, nl),
                         Linf=rep(linf, nl), MK=rep(mk, nl), L=Prop$L, Ll=Prop$Ll,
                         Lu=Prop$Lu, nS=rep(ns, nl), Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl,
                         pl.obs=Prop$pl.obs, Scen.p=rep(scp, nl))
        COMP.ub <- rbind(COMP.ub, Q2)
        Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                           L50=l50, L95r=l95r, nl=nl, wl=w, X=X, Linf=linf, MK=mk,
                           nS=ns, Scen.p=scu)
        PAR.ub <- rbind(PAR.ub, Par2)
        for (wr in WRB2){
          if (wr>0){
            lmin <- Comp.ub$Ll[1]
            lmax <- Comp.ub$Lu[nl]
            nlr <- ceiling((lmax - lmin)/wr)
            B <- lmin + (0:nlr)*wr
          }else{
            ima <- min(which(B.rls>=Comp.ub$Lu[nl]))
            B <- B.rls[1:ima]
          }
          if (r==1){scr <- scr + 1}else{
            Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr  & COMP.rb$Rep==1,]
            scr <- Comp.rb2$Scen.r[1]
          }
          Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
          nlr <- dim(Comp.rb)[1]
          Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                           CV=rep(cv, nlr), L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                           L50=rep(l50, nlr), L95r=rep(l95r, nlr), nl.rb=rep(nlr, nlr),
                           wl.rb=rep(wr, nlr), nS=rep(ns, nlr), Linf=rep(linf, nlr),
                           MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                           Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr), X=rep(X, nlr),
                           nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
          COMP.rb <- rbind(COMP.rb, Q3)
          Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, CV=cv, L0r=l0r, s.amax=sa,
                             nl.rb=nlr, wl.rb=wr, nS=ns, Scen.ub=scu, X=X, nl.ub=nl,
                             wl.ub=w, Linf=linf, MK=mk, L50=l50, L95r=l95r)
          PAR.rb <- rbind(PAR.rb, Par3)
        }
      }
    }
  }
}

#Reorder based on scenarios and replicates
PAR.ub2 <- PAR.ub[order(PAR.ub$Scen.ub, PAR.ub$Rep),]
PAR.rb2 <- PAR.rb[order(PAR.rb$Scen.rb, PAR.rb$Rep),]
COMP.ub2 <- COMP.ub[order(COMP.ub$Scen.ub, COMP.ub$Rep),]
COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]

# write.csv(PAR.ub2, "parub-cr-sel-order.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-cr-sel-order.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-cr-sel-order.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-cr-sel-order.csv", row.names=FALSE)
