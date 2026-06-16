data {
  int<lower=1> N;  // # observations
  int<lower=1> S;  // # species
  vector[N] phi; // photosynthetic efficiency
  vector[N] age;  // leaf age
  vector[N] vpd;  // vapour pressure deficit
  vector[N] par;  // photosynthetically active radiation
  array[N] int<lower=1, upper=S> sp ; // species
  int<lower=1> N_p;  // # predictions
}
parameters {
  real<lower=0> a; // scale parameter
  real<lower=0, upper=10000> b; // time parameter
  real<lower=1, upper=2> k; // shape parameter
  vector<lower=0>[S] a_s; // scale parameter
  vector<lower=10, upper=10000> [S] b_s; // time parameter
  vector<lower=1, upper=2>[S] k_s; // shape parameter
  real<lower=0> sigma; // residual error
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real<lower=0> sigma_k;
}
transformed parameters {
  vector[N] mu_age = a_s[sp] .* pow(age ./ b_s[sp], k_s[sp]-1) .* exp(-pow(age ./ b_s[sp], k_s[sp]));
}
model {
  phi ~ normal(mu_age, sigma);
  a_s ~ normal(a, sigma_a);
  b_s ~ lognormal(log(b), sigma_b);
  k_s ~ normal(k, sigma_k);
}
generated quantities {
  real rmsep = sqrt(mean(square(mu_age-phi)));
  real brsq = variance(mu_age) / (variance(mu_age) + sigma^2);
}
