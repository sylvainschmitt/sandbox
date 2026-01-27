data {
  int<lower=0> N; // #obs
  int<lower=0> Y; // #years
  int<lower=0> P; // #plots
  vector[P] y0;
  vector[N] y;
  array[N] int<lower=0, upper=Y> ysd;
  array[N] int<lower=0, upper=P> plot;
}
parameters {
  vector[P] alpha;
  vector[P] beta;
  vector<lower=0> [P]sigma;
}
transformed parameters {
  matrix[P,Y] mu;
  mu[,1] = y0;
  for(t in 2:Y)
    for(p in 1:P)
      mu[p,t] = alpha[p] + beta[p]*mu[p,(t-1)];
}
model {
  for(i in 1:N)
    y[i] ~ normal(mu[plot[i],ysd[i]], sigma[plot[i]]);
}
generated quantities {
  vector[N] y_p;
  for(i in 1:N)
    y_p[i] = normal_rng(mu[plot[i],ysd[i]], sigma[plot[i]]);
}
