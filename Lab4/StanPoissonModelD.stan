data {
  int<lower=0> N; // Number of observations
  int c[N];
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  real<lower=0> x[N]; // Data points
}

transformed parameters {
  real lambda[N]; 
  lambda = exp(x);
}

model {
  //likelihood
  for(n in 1:N) 
    c[n] ~ poisson(lambda[n]); // Poisson
      
  //prior
  mu ~ normal(10,2);
  phi ~ normal(0, 1);
  sigma ~ gamma(1,1000);
  
  for(n in 2:N)
    x[n] ~ normal(mu + phi * x[n-1], sigma);
    
}
