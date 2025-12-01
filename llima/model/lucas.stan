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
  real<lower=10, upper=100> mu_alpha;
  real<lower=0> sigma_alpha;
  vector<lower=0, upper=0.5>[c] lambda_c; // recovery rate
  real<lower=0, upper=0.5> mu_lambda;
  real<lower=0> sigma_lambda;
  real<lower=0> sigma; // residual variation
}
transformed parameters {
  vector[n] rec = 1 - exp(-lambda_c[cell] .* t);
  vector[n] mu = alpha_c[cell] + (100 - alpha_c[cell]) .* rec;
}
model {
  y ~ normal(mu, sigma);
  alpha_c ~ normal(mu_alpha, sigma_alpha);
  lambda_c ~ normal(mu_lambda, sigma_lambda);
}
generated quantities {
  vector[30] y_pred = mu_alpha + (100 - mu_alpha) .* (1 - exp(-mu_lambda .* t_pred));
}
