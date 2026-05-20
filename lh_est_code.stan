// Estimate M/K, Linf, and L50 and L95

//Input data
data {
  int<lower=1>        nl; // number of size bins
  vector[nl]         LOW; // vector of lower bounds of size bins
  vector[nl]         MID; // vector of mid value of size bins
  vector[nl]          UP; // vector of upper bounds of size bins
  real<lower=0>     Lmin; // lower bound of minimum size bin
  real<lower=0>     Lopt; // mid value of bin with maximum abundance
  real<lower=0>     Lmax; // upper bound of maximum size bin
  int<lower=0>         N; // total abundance
  vector[nl]          Nl; // Vector of abundance at size bin        
  int<lower=1>         X; // number of relative ages
  real<lower=0>       sa; // proportion of cohort surviving to maximum age
  real<lower=0>       Md; // natural mortality rate per relative time
  real<lower=0>       CV; // CV of size at age
}

//Model parameters
parameters {
  real<lower=0.01, upper=6>                   MK; // M/K ratio
  real<lower=(Lmin+Lmax)/2, upper=1.5*Lmax> Linf; // Asymptotic size
  real<lower=0.1, upper=Lmax>                L50; // length at 50% observed
  real<lower=1.0001, upper=4>               L95r; // length at 50% observed relative to L50
  real<lower=0.001>                       theta0; // precision parameter
}

//Model likelihood
model {
  real               alp;
  real                La;
  real               sdL;
  matrix[nl,X]        PP; // Age-size transition matrix
  vector[X]           Nx;
  vector[nl]      Nreal0;
  vector[nl]       Nreal;
  vector[nl]          Sl;
  vector[nl]       Nobs0;
  vector[nl]        Nobs;
  vector[nl]       Nobs2;
  real           theta_i;
  real                 R;
  real               eps; // miss-classification rate
  

  //Set priors for parameters
  theta0 ~ exponential(0.005);
  MK     ~ normal(1.5, 3);
  Linf   ~ normal(0.9*Lmax, 0.3*Lmax);
  L50    ~ normal((Lmin+Lopt)/2, (Lopt-Lmin)/2);
  L95r   ~ normal(1.8, 1.5);
  
  //Age-size transition matrix in population
  for (x in 1:X){
    La = Linf*(1-sa^(1/MK*(x-1)/(X-1)));
    sdL = CV*La;
    if (sdL==0){
      sdL = 0.001;
      }
    PP[1,x] = normal_cdf(UP[1], La, sdL);
    for (i in 2:(nl-1)){
      PP[i,x] = normal_cdf(UP[i], La, sdL) - normal_cdf(LOW[i], La, sdL);
    }
    PP[nl,x] = 1 - normal_cdf(LOW[nl], La, sdL);
    //PP[,x] = PP[,x]/sum(PP[,x]);
  }
  
  //Age composition
  R = 1;
  Nx[1] = R;
  for (x in 2:X){
    Nx[x] = Nx[x-1]*exp(-Md);
  }
  
  //Population size distribution
  Nreal0 = (Nx' * PP')';
  Nreal = Nreal0/sum(Nreal0);
  
  //Observed size distribution
  alp = log(19)/(L95r*L50-L50);
  for (i in 1:nl){
    Sl[i] = 1/(1 + exp(-alp*(MID[i]-L50))); //Observation selectivity
    Nobs0[i] = Sl[i]*Nreal[i];
  }
  Nobs = Nobs0/sum(Nobs0);

  // // add missclassification
  // //Nobs2 = Nobs;
  // eps = 0.01;
  // for (i in 1:nl) {
  //   Nobs2[i] = eps + (1.0 - nl*eps)*Nobs[i];    
  // }
  Nobs2 = Nobs;
  
  target += lgamma(theta0) + lgamma(N+1) - lgamma(N+theta0); 
  
  //Likelihood according DM distribution
  for (i in 1:nl) { 
    theta_i = theta0*Nobs2[i];
    target += lgamma(Nl[i] + theta_i) - lgamma(theta_i) - lgamma(Nl[i]+1);
  }
}
