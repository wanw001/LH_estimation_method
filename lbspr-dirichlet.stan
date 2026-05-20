// Estimating M/K and Linf from length frequency data from underwater survey 

// Input data
data {
  // known fish life-history traits
  int<lower=1>        nl; // number of size classes
  vector[nl]         LOW; // lower bounds of size classes
  vector[nl]         MID; // mid size in each class
  vector[nl]          UP; // upper bounds of size classes
  real<lower=0>     Lmin;
  real<lower=0>     Lopt;
  real<lower=0>     Lmax; // maxmimum size
  int<lower=0>         N; // Total abundance of fish
  vector[nl]          Nl; // Vector of abundance at size        
  int<lower=1>         X; // maximum relative age +1
  real<lower=0>       sa; // proportion of cohort surviving to maximum age
  real<lower=0>       Md; // Natural mortality rate per relative time
  real<lower=0>       CV; // CV of size at age
}

// Model parameters
parameters {
  // estimated fish life-history parameters
  real<lower=0.1, upper=6>                    MK; // M/K ratio
  real<lower=(Lmin+Lmax)/2, upper=1.5*Lmax> Linf; // Infinite fish size
  real<lower=0.1, upper=Lmax>                L50; // length when 50% observe
  real<lower=1.0001, upper=4>               L95r; // steepness of observability w.r.t. length
  // variance parameter for observed abundances
  real<lower=0.001>       theta0;
}

// Model likelihood
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
  real            SNobs0;
  real           theta_k;
  real                 R;
  real               eps; // miss-classification rate
  

  // set priors for parameters
  theta0 ~ exponential(0.005);
  MK     ~ normal(1.5, 3);
  Linf   ~ normal(0.9*Lmax, 0.3*Lmax);
  L50    ~ normal((Lmin+Lopt)/2, (Lopt-Lmin)/2);
  L95r   ~ normal(1.8, 1.5);
  
  // Age-size transition matrix in population
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
  
  //Size composition (real)
  Nreal0 = (Nx' * PP')';
  Nreal = Nreal0/sum(Nreal0);
  
  //Size composition (observed)
  alp = log(19)/(L95r*L50-L50);
  for (i in 1:nl){
    Sl[i] = 1/(1 + exp(-alp*(MID[i]-L50))); //Observation selectivity
    Nobs0[i] = Sl[i]*Nreal[i];
  }
  SNobs0 = 0;
  for (i in 1:nl){
    SNobs0 = SNobs0 + Nobs0[i];
  }
  for (i in 1:nl){
    Nobs[i] = Nobs0[i]/SNobs0;
  }

  // //add missclassification
  // eps = 0.0001;
  // for (i in 1:nl) {
  //   Nobs2[i] = eps + (1.0 - nl*eps)*Nobs[i];    
  // }
  Nobs2 = Nobs;
  
  target += lgamma(theta0) + lgamma(N+1) - lgamma(N+theta0); 
  
  //Likelihood using DM distribution
  for (i in 1:nl) { 
    theta_k = theta0*Nobs2[i];
    target += lgamma(Nl[i] + theta_k) - lgamma(theta_k) - lgamma(Nl[i]+1);
  }
}
