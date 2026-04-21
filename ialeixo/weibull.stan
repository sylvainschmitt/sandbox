data {
  int<lower=1> N;  // # observations
  vector[N] t;  // leaf age
  vector[N] p; // photosynthetic efficiency
}
parameters {
  real<lower=0> a; // scale parameter
  real<lower=0, upper=365> b; // time parameter
  real<lower=1, upper=2> k; // shape parameter
  real<lower=0> sigma; // residual error
}
transformed parameters {
  vector[N] mu = a*pow(t/b, k-1).*exp(-pow(t/b, k));
}
model {
  p ~ normal(mu, sigma);
}
generated quantities {
  real tmax = b*((k-1)/k)^(1/k);
  real pmax = a*pow(tmax/b, k-1)*exp(-pow(tmax/b, k));
}
