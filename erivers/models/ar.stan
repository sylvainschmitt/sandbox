data {
  int<lower=0> n;
  vector[n] y;
  vector[n] clim;
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real<lower=0> sigma;
}
transformed parameters {
  vector[n - 1] mu = alpha + beta * y[1:(n - 1)] + gamma * clim[2:n];
}
model {
  y[2:n] ~ normal(mu, sigma);
}
generated quantities {
  vector[n] y_p;
  y_p[1] = y[1];
  y_p[2:n] = to_vector(normal_rng(mu, sigma));
}
