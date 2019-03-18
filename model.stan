
data {
     int<lower=0> Nobs;
     int Npos [Nobs];
     int Ntrials [Nobs];
     real max_elevation;
     real z[Nobs];
     
}

parameters {
           real<lower=0,upper=2> a;
           real<lower=0,upper=2> z_half;
           real<lower=0,upper=2> c;
}

transformed parameters {
  real P [Nobs];
  
 for (i in 1:Nobs){
    P[i] = ((1+exp(a))/(1+exp(a*(1-z[i]/z_half))))*(1-c) + c;
 }


}

model {
  a ~ normal (-10, 2);
  z_half ~ uniform (0, max_elevation);
  c ~ uniform (0, 1);
  
  for (i in 1:Nobs){
  Npos[i] ~ binomial(Ntrials[i], P[i]) ;
  }
  
}


