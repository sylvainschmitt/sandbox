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
  real av;
  real bv;
  real<lower=0> sigma; // residual error
}
transformed parameters {
  vector[N] mu_vpd = av+bv*vpd;
}
model {
  phi ~ normal(mu_vpd, sigma);
}
generated quantities {
  real rmsep = sqrt(mean(square(mu_vpd-phi)));
  real brsq = variance(mu_vpd) / (variance(mu_vpd) + sigma^2);
}
