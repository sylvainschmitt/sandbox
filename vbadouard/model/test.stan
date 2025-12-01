data {
  int<lower=0> n;
  vector[n] x;
  array[n] int y;
}
parameters {
  real a;
  real b;
}
transformed parameters {
  vector[n] theta = inv_logit(a*x+b);
}
model {
  y ~ bernoulli(theta);
}
generated quantities {
  array[n] int y_p = bernoulli_rng(theta);
}
