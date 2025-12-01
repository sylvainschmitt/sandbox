data {
  int<lower=0> n;
  vector[n] light;
  array[n] int presence;
}
parameters {
  real<upper=0> a;
  real<lower=0, upper=100> O;
  real gamma;
}
model {
  presence ~ bernoulli(inv_logit(a*(light/100 - O/100)^2 + gamma));
}
