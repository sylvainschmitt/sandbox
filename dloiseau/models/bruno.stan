data {
  int<lower=0> N; // # obs
  vector[N] dbh; // tree diameters
  vector[N] agr; // tree growths
}
parameters {
  real<lower=1, upper=10> gmax;
  real<lower=1, upper=1000> dopt;
  real<lower=0, upper=3> k;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu = gmax*exp(-0.5*(log(dbh/dopt)/k)^2);
}
model {
  log(agr) ~ normal(log(mu), sigma);
}
