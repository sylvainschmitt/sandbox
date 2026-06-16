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
transformed data {
  vector[N_p] age_p;
  for(n in 1:N_p/10)
    age_p[n] = n; 
  age_p = age_p*10;
}
parameters {
  real<lower=0, upper=1> ap; // intercept par
  real<lower=-10, upper=10> bp_log; // slope par
  // real av; // intercept vpd
  real bv;  // slope vpd
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
  real bp=exp(bp_log);
  vector[N] mu_par = ap*exp(-bp*par);
  vector[N] mu_vpd = bv*vpd;
  vector[N] mu_age = a_s[sp] .* pow(age ./ b_s[sp], k_s[sp]-1) .* exp(-pow(age ./ b_s[sp], k_s[sp]));
  vector[N] mu = mu_age + mu_vpd + mu_par;
}
model {
  phi ~ normal(mu, sigma);
  a_s ~ normal(a, sigma_a);
  b_s ~ lognormal(log(b), sigma_b);
  k_s ~ normal(k, sigma_k);
}
generated quantities {
  real rmsep = sqrt(mean(square(mu-phi)));
  real brsq = variance(mu) / (variance(mu) + sigma^2);
  vector[S] tmax_s = b_s .* ((k_s-1) ./ k_s)^(1/k_s);
  vector[S] pmax_s = a_s .* pow(tmax_s ./ b, k_s-1) .* exp(-pow(tmax_s ./ b_s, k_s));
  array[S] vector[N_p] y_sp;
  for(s in 1:S)
    y_sp[s] = a_s[s] .* pow(age_p ./ b_s[s], k_s[s]-1) .* exp(-pow(age_p ./ b_s[s], k_s[s]));
}
