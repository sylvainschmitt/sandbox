data {
  int<lower=1> N;  // # observations
  int<lower=1> S;  // # species
  vector[N] phi; // photosynthetic efficiency
  vector[N] age;  // leaf age
  vector[N] vpd;  // vapour pressure deficit
  vector[N] par;  // photosynthetically active radiation
  array[N] int<lower=1, upper=S> sp ; // species
}
parameters {
  real<lower=0, upper=1> ap; // intercept par
  real<lower=-10, upper=10> bp_log; // slope par
  real<lower=0> sigma; // residual error
}
transformed parameters {
  real bp=exp(bp_log);
  vector[N] mu_par = ap*exp(-bp*par);
}
model {
  phi ~ normal(mu_par, sigma);
}
generated quantities {
  real rmsep = sqrt(mean(square(mu_par-phi)));
  real brsq = variance(mu_par) / (variance(mu_par) + sigma^2);
}
