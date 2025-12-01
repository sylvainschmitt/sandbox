data {
  int<lower=0> n_obs; // # obs
  int<lower=0> n_sp; // # species
  vector[n_obs] dbh; // diameters at breast height
  vector[n_obs] cr; // crown radii
  array[n_obs] int<lower=1, upper=n_sp> species; // species index
  real cr_omax; // maximum crown radius at dbh=1cm
}
parameters {
  vector<lower=0>[n_sp] b;
  vector<upper=cr_omax>[n_sp] ap;
  real<lower=0> sigma;
  real mu_ap;
  real mu_b;
  real<lower=0> sigma_ap;
  real<lower=0> sigma_b;
}
transformed parameters {
  vector[n_sp] a = ap + b*log(100);
}
model {
  log(cr) ~ normal(ap[species] + b[species] .* log(dbh*100), sigma);
  ap ~ normal(mu_ap, sigma_ap);
  b ~ normal(mu_b, sigma_b);
  mu_b ~ normal(0.75, 1);
  mu_ap ~ normal(2.15-0.75*log(100), 1);
  sigma ~ std_normal();
  sigma_ap ~ std_normal();
  sigma_b ~ std_normal();
}

