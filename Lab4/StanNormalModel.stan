data {
  int<lower=0> T; // Number of observations 
  real x[T];
}

parameters {
  real mu; 
  real<lower=0> sigma2; 
  real<lower=-1, upper=1> phi;
}
 
model {
  //priors
  // mu ~ normal(10,2);
  // phi ~ uniform(-1,1);
  // sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 
  // 
  //likelihood
  for (t in 2:T)
    x[t] ~ normal(mu + phi * x[t-1], sqrt(sigma2));
}
