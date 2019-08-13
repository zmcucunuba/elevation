data {
     int<lower=0> Nobs;
     int Npos [Nobs];
     int Ntrials [Nobs];
     real max_elevation;
     real z[Nobs];
     
}

parameters {
           real<lower=0,upper=1> a;
           real<lower=0,upper=max_elevation> b;
           real<lower=0,upper=1> c;
           real<lower=0,upper=1> p_zero;
}

transformed parameters {
  real P [Nobs];
  
 for (i in 1:Nobs){
    P[i] =  (( 1+exp( -a* b ) ) / (1+exp(a *( z[i] - b ) ))  * (1-c) + c) * p_zero;
  
 }


}

model {
  a ~ uniform (0, 1);
  b ~ uniform (0, max_elevation);
  c ~ uniform (0, 1);
  p_zero ~ uniform (0, 1);
  
  for (i in 1:Nobs){
  Npos[i] ~ binomial(Ntrials[i], P[i]) ;
  }
  
}
