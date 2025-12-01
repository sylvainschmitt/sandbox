data {
  int<lower=0> n; // # obs
  int<lower=0> c; // # cells
  vector[n] y; // forest attribute
  vector[n] t; // time
  array[n] int<lower=0, upper=n> cell; // cell index
}
transformed data {
  vector[30] t_pred;
  for(i in 1:30)
    t_pred[i] = i;
}
parameters {
  vector<lower=10, upper=100>[c] alpha_c; // post disturbance value
  real<lower=10, upper=100> alpha;
  vector<lower=0, upper=0.5>[c] lambda_c; // recovery rate
  real<lower=0, upper=0.5> lambda;
  vector<lower=0>[3] sigma; // residual variation
}
model {
  y ~ normal(alpha_c[cell] + (100 - alpha_c[cell]) .* (1 - exp(-lambda_c[cell] .* t)), sigma[1]);
  alpha_c ~ normal(alpha, sigma[2]);
  lambda_c ~ normal(lambda, sigma[3]);
}
generated quantities {
  vector[c] tau_c=log(10)/lambda_c;
  real tau=log(10)/lambda;
  vector[30] y_pred = alpha + (100 - alpha) .* (1 - exp(-lambda .* t_pred));
}
