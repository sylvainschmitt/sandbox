data {
  int<lower=0> N; // # obs
  int<lower=0> G; // # genera
  vector[N] h; // tree heights
  vector[N] d; // tree diameters
  array[N] int<lower=1, upper=G> genus; // genus index
}
parameters {
  // max tree height
  vector<lower=1, upper=100>[G] alpha_g; // max tree hieght per genus
  real<lower=1, upper=100> alpha; // mean max tree height of the community
  real<lower=0> sigma_a; // variation of genus max tree height
  // growth in height speed
  vector<lower=1, upper=300>[G] beta_g; // growth in height speed per genus
  real<lower=1, upper=300> beta; // mean growth in height speed of the community
  real<lower=0> sigma_b; // variation of genus growth in height speed
  // residual variation
  real<lower=0> sigma; // residual variation
}
transformed parameters {
  vector[N] mu = (alpha_g[genus] .* d) ./ (beta_g[genus] + d); // modelled tree height
}
model {
  log(h) ~ normal(log(mu), sigma);
  alpha_g ~ normal(alpha, sigma_a); 
  beta_g ~ normal(beta, sigma_b); 
  [sigma_a, sigma_b, sigma] ~ exponential(1);
}
