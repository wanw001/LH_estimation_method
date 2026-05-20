#### Simulating length data under pulsed recruitment. 

setwd("C:/Users/wanwank/OneDrive - University of Tasmania/PhD_Wanwan-Kurniawan/Code/Ch1/Github/LH_estimation_method")

library(ggplot2)
library(arrow)

source("DGM_acr.R")

#### Parameters
LINF0 <- c(10, 20, 30, 40, 55, 70, 90) #Already consider the maximum bin in RLS data
MK <- c(0.1, 0.5, 0.8, 1, 1.5, 2.2, 2.3, 3, 3.5)
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
pr <- PR[2]
#nl <- 100  #Original number of size bins (unbinned sizes)
w <- 0.5
NS <- c(200, 500, 1000, 5000, 20000, 100000) #Sample size
Red <- 30 #Replicates of random data generation
WRB <- c(1, 2, 4, 8, 12, 0) #Bin width (Binned)
B.rls <- c(1.25, 3.75, 6.25, 8.75, 11.25, 13.75, 17.5, 22.5, 27.5, 32.5, 37.5, 45, 
           56.25, 68.75, 81.25, 118.75, 131.25)
mr <- 2:4 #Month of recruitment
MO <- c(3, 9) #Month of observation for 1 time sampling
NM <- c(1, 2, 3, 4, 6, 12) #Number of sampling

## Sensitivity analysis for frequency (and month) of observation
linf <- 40
l50 <- 5
l95r <- 1.5
l0r <- 0
sa <- 0.01
cv <- 0.1

# Generate the list of parameter and proportion for each month
scp <- 0
PAR <- PROP <- NULL
for (mk in MK){
  for (mo in 1:12){
    scp <- scp + 1
    Prop <- sizeProp.acr(linf=linf, mk=mk, l50=l50, l95r=l95r, l0r=l0r, sa=sa, 
                         cv=cv, w=w, mr=mr, mo=mo)
    #w <- Prop$Lu[1]-Prop$Ll[1]; 
    nl <- nrow(Prop)
    Q <- data.frame(Scen.p=rep(scp, nl), perR=rep(pr, nl), mo=rep(mo, nl), CV=rep(cv, nl),
                    L0r=rep(l0r, nl), s.amax=rep(sa, nl), Linf=rep(linf, nl), MK=rep(mk, nl),
                    L50=rep(l50, nl), L95r=rep(l95r, nl), L=Prop$L, Ll=Prop$Ll, 
                    Lu=Prop$Lu, pl=Prop$pl, pl.obs=Prop$pl.obs, nl=rep(nl, nl), 
                    wl=rep(w, nl))
    PROP <- rbind(PROP, Q)
    Par <- data.frame(Scen.p=scp, perR=pr, mo=mo, CV=cv, L0r=l0r, s.amax=sa, 
                      Linf=linf, MK=mk, L50=l50, L95r=l95r, nl=nl, wl=w)
    PAR <- rbind(PAR, Par)
  }
}

# write.csv(PAR, "parp-ar-freq-month-order_rev02.csv", row.names=FALSE)
# write.csv(PROP, "prop-ar-freq-month-order_rev02.csv", row.names=FALSE)

#### Generate the random abundance of fish at size (unbinned and binned, 30 replicates)
WRB2 <- WRB[c(2, 3, 6)]
NS2 <- NS[3:6]
Re <- 30
PROP <- read.csv("prop-ar-freq-month-order_rev02.csv")
COMP.ub <- COMP.rb <- NULL
PAR.ub <- PAR.rb <- NULL
scu <- scr <- 0
for (nm in NM){
  if (nm==1){
    for (mo in MO){
      for (scp in unique(PROP[PROP$mo==mo,]$Scen.p)){
        Prop <- PROP[PROP$Scen.p==scp,]; #Prop <- Prop[order(Prop$L),]
        pr <- Prop$perR[1]; cv <- Prop$CV[1]; l0r <- Prop$L0r[1]; 
        sa <- Prop$s.amax[1]; linf <- Prop$Linf[1]; mk <- Prop$MK[1]; l50 <- Prop$L50[1]
        l95r <- Prop$L95r[1]; nl <- Prop$nl[1]; w <- Prop$wl[1]
        for (ns in NS2){
          scu <- scu + 1
          for (r in 1:Re){
            # if (r==1){scu <- scu +1}else{
            #   Comp.ub2 <- COMP.ub[COMP.ub$Scen.p==scp & COMP.ub$nm==nm & COMP.ub$mo==mo 
            #                       & COMP.ub$nS==ns  & COMP.ub$Rep==1,]
            #   scu <- Comp.ub2$Scen.u[1]
            # }
            Comp.ub <- sizeComp(r=r, Prop=Prop, ns=ns)
            Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                             nm=rep(nm, nl), mo=rep(mo, nl), CV=rep(cv, nl),
                             L0r=rep(l0r, nl), s.amax=rep(sa, nl), nl=rep(nl, nl),
                             wl=rep(w, nl), nS=rep(ns, nl), L50=rep(l50, nl), L95r=rep(l95r, nl), 
                             Linf=rep(linf, nl), MK=rep(mk, nl), L=Comp.ub$L, Ll=Comp.ub$Ll, 
                             Lu=Comp.ub$Lu, Nl.obs=Comp.ub$Nl.obs, pl=Prop$pl, pl.obs=Prop$pl.obs, 
                             Scen.p=rep(scp, nl))
            COMP.ub <- rbind(COMP.ub, Q2)
            Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, nm=nm, mo=mo, CV=cv, L0r=l0r,
                               s.amax=sa, nl=nl, wl=w, nS=ns, L50=l50, L95r=l95r, Linf=linf, 
                               MK=mk, Scen.p=scp)
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
              nlr <- nrow(Comp.rb)
              Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                               nm=rep(nm, nlr), mo=rep(mo, nlr), CV=rep(cv, nlr), 
                               L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                               nl.rb=rep(nlr, nlr), wl.rb=rep(wr, nlr), nS=rep(ns, nlr), 
                               L50=rep(l50, nlr), L95r=rep(l95r, nlr), Linf=rep(linf, nlr),
                               MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                               Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr),
                               nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
              COMP.rb <- rbind(COMP.rb, Q3)
              Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, nm=nm, mo=mo, CV=cv, 
                                 L0r=l0r, s.amax=sa, nl.rb=nlr, wl.rb=wr, nS=ns, 
                                 L50=l50, L95r=l95r, Linf=linf, MK=mk, Scen.ub=scu, 
                                 nl.ub=nl, wl.ub=w)
              PAR.rb <- rbind(PAR.rb, Par3)
            }
          }
        }
      }
    }
  }else{
    mo <- 0
    MO2 <- (1:nm)*12/nm - 2
    for (i in 1:length(MO2)){
      if (MO2[i]<=0){MO2[i] <- MO2[i]+12}
    }
    for (mk in MK){
      for (ns in NS2){
        scu <- scu + 1
        for (r in 1:Re){
          nl1 <- 1
          for (mo2 in MO2){
            nl0 <- nrow(PROP[PROP$mo==mo2 & PROP$MK==mk,])
            if (nl0>nl1){mo1 <- mo2; nl1 <- nl0} 
          }
          #mo2 <- MO2[1]
          Prop <- PROP[PROP$mo==mo1 & PROP$MK==mk,]; #Prop <- Prop[order(Prop$L),]
          scp <- Prop$Scen.p[1]; pr <- Prop$perR[1]; cv <- Prop$CV[1]; l0r <- Prop$L0r[1];
          sa <- Prop$s.amax[1]; linf <- Prop$Linf[1]; l50 <- Prop$L50[1];
          l95r <- Prop$L95r[1]; nl <- Prop$nl[1]; w <- Prop$wl[1]
          Comp.ub <- sizeComp(r=r, Prop=Prop, ns=round(ns/nm))
          for (mo2 in MO2[!MO2 %in% mo1]){
            Propi <- PROP[PROP$mo==mo2 & PROP$MK==mk,]
            if (mo2==MO2[length(MO2)]){nsi <- ns - sum(Comp.ub$Nl.obs)}else{nsi=round(ns/nm)}
            Comp.ubi <- sizeComp(r=r, Prop=Propi, ns=nsi)
            nl2 <- nrow(Comp.ubi)
            if (nl2<nl1){
              Row <- Comp.ub[(nl2+1):nl1,]; Row$Nl.obs <- 0
              Comp.ubi <- rbind(Comp.ubi, Row) 
            }
            Comp.ub$Nl.obs <- Comp.ub$Nl.obs + Comp.ubi$Nl.obs
          }
          ns2 <- sum(Comp.ub$Nl.obs)
          Q2 <- data.frame(Scen.ub=rep(scu, nl), Rep=rep(r, nl), perR=rep(pr, nl),
                           nm=rep(nm, nl), mo=rep(mo, nl), CV=rep(cv, nl),
                           L0r=rep(l0r, nl), s.amax=rep(sa, nl), nl=rep(nl, nl),
                           wl=rep(w, nl), nS=rep(ns2, nl), L50=rep(l50, nl), L95r=rep(l95r, nl),
                           Linf=rep(linf, nl), MK=rep(mk, nl), L=Comp.ub$L, Ll=Comp.ub$Ll,
                           Lu=Comp.ub$Lu, Nl.obs=Comp.ub$Nl.obs, pl=rep(0, nl), pl.obs=rep(0, nl),
                           Scen.p=rep(scp, nl))
          COMP.ub <- rbind(COMP.ub, Q2)
          Par2 <- data.frame(Scen.ub=scu, Rep=r, perR=pr, nm=nm, mo=mo, CV=cv, L0r=l0r,
                             s.amax=sa, nl=nl, wl=w, nS=ns2, L50=l50, L95r=l95r, Linf=linf,
                             MK=mk, Scen.p=scp)
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
              Comp.rb2 <- COMP.rb[COMP.rb$Scen.u==scu & COMP.rb$wl.rb==wr & COMP.rb$Rep==1,]
              scr <- Comp.rb2$Scen.r[1]
            }
            Comp.rb <- sizeComp.rebin(B=B, Comp=Comp.ub)
            nlr <- nrow(Comp.rb)
            Q3 <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                             nm=rep(nm, nlr), mo=rep(mo, nlr), CV=rep(cv, nlr), 
                             L0r=rep(l0r, nlr), s.amax=rep(sa, nlr),
                             nl.rb=rep(nlr, nlr), wl.rb=rep(wr, nlr), nS=rep(ns, nlr), 
                             L50=rep(l50, nlr), L95r=rep(l95r, nlr), Linf=rep(linf, nlr),
                             MK=rep(mk, nlr),  L=Comp.rb$L, Ll=Comp.rb$Ll, Lu=Comp.rb$Lu,
                             Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr),
                             nl.ub=rep(nl, nlr), wl.ub=rep(w, nlr))
            COMP.rb <- rbind(COMP.rb, Q3)
            Par3 <- data.frame(Scen.rb=scr, Rep=r, perR=pr, nm=nm, mo=mo, CV=cv, 
                               L0r=l0r, s.amax=sa, nl.rb=nlr, wl.rb=wr, nS=ns, 
                               L50=l50, L95r=l95r, Linf=linf, MK=mk, Scen.ub=scu, 
                               nl.ub=nl, wl.ub=w)
            PAR.rb <- rbind(PAR.rb, Par3)
          }
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

# write.csv(PAR.ub2, "parub-ar-freq-month-order_rev02.csv", row.names=FALSE)
# write.csv(PAR.rb2, "parrb-ar-freq-month-order_rev02.csv", row.names=FALSE)
# write.csv(COMP.ub2, "compub-ar-freq-month-order_rev02.csv", row.names=FALSE)
# write.csv(COMP.rb2, "comprb-ar-freq-month-order_rev02.csv", row.names=FALSE)

### Rebin the size frequency data
COMP <- data.frame(read_parquet("comp-unbin-acr-Linf45.parquet"))
Scen.ub <- sort(unique(COMP$Scen.ub))
Rep <- sort(unique(COMP$Rep))
# Scen.ub <- c(1,3)
# Rep <- c(1, 7)
COMP.rb <- NULL
QUANT.rb <- NULL
scr <- 0
for (scu in Scen.ub){
  for (r in Rep){
    Comp <- COMP[COMP$Scen.ub==scu & COMP$Rep==r,]; Comp <- Comp[order(Comp$L),]
    pr <- Comp$perR[1]; nm <- Comp$nm[1]; mo <- Comp$mo[1]; cv <- Comp$CV[1]; l0r <- Comp$L0r[1];
    sa <- Comp$s.amax[1]; wu <- Comp$wl[1]; linf <- Comp$Linf[1]; mk <- Comp$MK[1];
    l50 <- Comp$L50[1]; l95r <- Comp$L95r[1]; ns <- Comp$nS[1]; nlu <- Comp$nl[1]
    if (linf==LINF[2]){
      for (wr in WRB){
        if (wr>0){
          lmin <- Comp$Ll[1]
          lmax <- Comp$Lu[nlu]
          nlr <- ceiling((lmax - lmin)/wr)
          B <- lmin + (0:nlr)*wr
        }else{
          ima <- min(which(B.rls>=Comp$Lu[nlu]))
          B <- B.rls[1:ima]
        }
        Comp.rb <- sizeComp.rebin(B=B, Comp=Comp)
        nlr <- dim(Comp.rb)[1]
        if (r==1){scr <- scr + 1}else{
          Comp2 <- COMP.rb[COMP.rb$Scen.ub==scu & COMP.rb$Rep==1 & COMP.rb$wl.rb==wr,]
          scr <- Comp2$Scen.rb[1]
        }
        Q <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                        nm=rep(nm, nlr), mo=rep(mo, nlr), CV=rep(cv, nlr), L0r=rep(l0r, nlr),
                        s.amax=rep(sa, nlr), wl.rb=rep(wr, nlr), nl.rb=rep(nlr, nlr),
                        nS=rep(ns, nlr), Linf=rep(linf, nlr), MK=rep(mk, nlr),
                        L50=rep(l50, nlr), L95r=rep(l95r, nlr), L=Comp.rb$L, Ll=Comp.rb$Ll,
                        Lu=Comp.rb$Lu, Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr),
                        wl.ub=rep(wu, nlr), nl.ub=rep(nlu, nlr))
        COMP.rb <- rbind(COMP.rb, Q)
        Qu <- data.frame(Scen.rb=scr, Rep=r, perR=pr, nm=nm, mo=mo, CV=cv, L0r=l0r,
                         s.amax=sa, wl.rb=wr, nl.rb=nlr, nS=ns, Linf=linf, MK=mk,
                         L50=l50, L95r=l95r, Scen.ub=scu, wl.ub=wu, nl.ub=nlu)
        QUANT.rb <- rbind(QUANT.rb, Qu)
      }
    }else{
      wr <- WRB[length(WRB)]
      ima <- min(which(B.rls>=Comp$Lu[nlu]))
      B <- B.rls[1:ima]
      Comp.rb <- sizeComp.rebin(B=B, Comp=Comp)
      nlr <- dim(Comp.rb)[1]
      if (r==1){scr <- scr + 1}else{
        Comp2 <- COMP.rb[COMP.rb$Scen.ub==scu & COMP.rb$Rep==1 & COMP.rb$wl.rb==wr,]
        scr <- Comp2$Scen.rb[1]
      }
      Q <- data.frame(Scen.rb=rep(scr, nlr), Rep=rep(r, nlr), perR=rep(pr, nlr),
                      nm=rep(nm, nlr), mo=rep(mo, nlr), CV=rep(cv, nlr), L0r=rep(l0r, nlr),
                      s.amax=rep(sa, nlr), wl.rb=rep(wr, nlr), nl.rb=rep(nlr, nlr),
                      nS=rep(ns, nlr), Linf=rep(linf, nlr), MK=rep(mk, nlr),
                      L50=rep(l50, nlr), L95r=rep(l95r, nlr), L=Comp.rb$L, Ll=Comp.rb$Ll,
                      Lu=Comp.rb$Lu, Nl.obs=Comp.rb$Nl.obs, Scen.ub=rep(scu, nlr),
                      wl.ub=rep(wu, nlr), nl.ub=rep(nlu, nlr))
      COMP.rb <- rbind(COMP.rb, Q)
      Qu <- data.frame(Scen.rb=scr, Rep=r, perR=pr, nm=nm, mo=mo, CV=cv, L0r=l0r,
                       s.amax=sa, wl.rb=wr, nl.rb=nlr, nS=ns, Linf=linf, MK=mk,
                       L50=l50, L95r=l95r, Scen.ub=scu, wl.ub=wu, nl.ub=nlu)
      QUANT.rb <- rbind(QUANT.rb, Qu)
    }
  }
}

dim(COMP.rb)[1]
dim(QUANT.rb)[1]

COMP.rb2 <- COMP.rb[order(COMP.rb$Scen.rb, COMP.rb$Rep),]
QUANT.rb2 <- QUANT.rb[order(QUANT.rb$Scen.rb, QUANT.rb$Rep),]
# write_parquet(COMP.rb, "comp-rebin-acr-Linf40.parquet")
# write.csv(COMP.rb2, "comp-rebin-acr-Linf45.csv", row.names=FALSE)
# write.csv(QUANT.rb2, "param-rebin-acr-Linf45.csv", row.names=FALSE)


