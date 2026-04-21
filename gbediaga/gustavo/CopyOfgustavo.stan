data {
  int<lower=0> N; // # obs
  int<lower=0> S; // # species
  int<lower=0> U; // # ups
  vector[N] h; // tree heights
  vector[N] d; // tree diameters
  array[N] int<lower=1, upper=S> species; // species index
  array[N] int<lower=1, upper=U> up; // ups index
}
parameters {
  // max tree height
  vector<lower=1, upper=100>[S] alpha_s; // max tree hieght per species
  real<lower=1, upper=100> alpha; // mean max tree height of the community
  real<lower=0> sigma_a; // variation of species max tree height
  // growth in height speed
  vector<lower=1, upper=300>[S] beta_s; // growth in height speed per species
  real<lower=1, upper=300> beta; // mean growth in height speed of the community
  real<lower=0> sigma_b; // variation of species growth in height speed
  // alpha beta correlation
  corr_matrix[2] R;
  // ups effect
  vector[U] u;  // random effect of ups
  real<lower=0> sigma_u; // random variations of ups
  // residual variation
  real<lower=0> sigma; // residual variation
}
transformed parameters {
  vector[N] mu; // modelled tree height
  array[S] vector[2] v; // vector of species alpha and beta
  cov_matrix[2] Sigma; // species alpha beta covariance matrix
  mu = (alpha_s[species] .* d) ./ (beta_s[species] + d);
  for(s in 1:S) v[s] = [alpha_s[s], beta_s[s]]';
  Sigma = quad_form_diag(R, [sigma_a, sigma_b]);
}
model {
  log(h) ~ normal(log(mu) + u[up], sigma);
  u ~ normal(0, sigma_u); 
  v ~ multi_normal([alpha, beta], Sigma);  
  R ~ lkj_corr(2);
  [sigma_a, sigma_b, sigma_u, sigma] ~ exp(1);
}
