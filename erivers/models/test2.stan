data {
  int<lower=0> N;
  int<lower=0> Y;
  real y0;
  vector[N] y;
  array[N] int<lower=1, upper=Y> ysd;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[Y] mu;
  mu[1] = y0;
  for(i in 2:Y)
    mu[i] = alpha + beta * mu[i-1];
}
model {
  y ~ normal(mu[ysd], sigma);
  alpha ~ normal(-5, 5);
  beta ~ normal(1,1);
}
generated quantities {
  vector[N+1] y_p;
  y_p[1] = y0;
  y_p[2:N+1] = to_vector(normal_rng(mu[ysd], sigma));
}
