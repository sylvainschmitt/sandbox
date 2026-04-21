data {
  int<lower=1> N_rec; // # obs recovery
  int<lower=1> N_pre; // # obs predisturbance
  int<lower=1> P; // # plots
  vector[N_rec] agb; // biomass
  vector[N_rec] time; // time
  vector[N_pre] agb_pre; // biomass in predisturbance
  array[N_rec] int<lower=1, upper=P> plot; // plot index
  array[N_pre] int<lower=1, upper=P> plot_pre; // plot index in predisturbance
}
parameters {
  // equlibirum biomass
  vector<lower=100, upper=1000>[P] agb_inf_p; // plot level
  // disturbance intensity
  vector<lower=0, upper=2>[P] alpha_p; // plot level
  real<lower=0, upper=2> alpha; // mean
  // recovery rate
  vector<lower=0, upper=0.5>[P] lambda_p; // plot level
  real<lower=0, upper=0.5> lambda; // mean
  // variations
  real<lower=0> sigma; // residual variation
  real<lower=0> sigma_pre; // prelogging variation
  real<lower=0> sigma_a; // disturbance intensity variation
  real<lower=0> sigma_l; // recovery rate variation
}
transformed parameters {
  vector[P] agb_0_p = agb_inf_p - agb_inf_p .* alpha_p; // postdisturbance biomass
  vector[N_rec] mu = agb_0_p[plot] + (agb_inf_p[plot] - agb_0_p[plot]) .* (1-exp(-lambda_p[plot] .* time)); // predicted agb
}
model {
  agb ~ normal(mu, sigma);
  agb_pre ~ normal(agb_inf_p[plot_pre], sigma_pre);
  alpha_p ~ normal(alpha, sigma_a);
  lambda_p ~ normal(lambda, sigma_l);
}
