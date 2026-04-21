data {
  int<lower=0> N; // # obs
  int<lower=0> Y; // # years
  int<lower=0> I; // # individuals
  int<lower=0> S; // # species
  vector[I] dbh0; // diameter at recruitment
  vector[S] dmax; // species maximum diameter
  vector[N] dbh; // tree diameters
  array[N] int<lower=1> year ; // time since recruitment
  array[N] int<lower=1> ind ; // individual
  array[I] int<lower=1, upper=S> indsp ; // species corresponding to individual
}
parameters {
  vector<lower=0.001, upper=10>[I] gmax_i ;
  vector<lower=0.001, upper=10>[S] gmax ;
  vector<lower=0.1, upper=1>[I] d_i;
  vector<lower=0.1, upper=1>[S] d;
  vector<lower=0.001, upper=3>[I] k_i;
  vector<lower=0.001, upper=3>[S] k;
  real<lower=0> sigma;
  real<lower=0> sigma_g;
  real<lower=0> sigma_d;
  real<lower=0> sigma_k;
}
transformed parameters {
  matrix[I,Y] DBH ;
  vector[N] mu;
  vector[I] dopt_i = d_i .* dmax[indsp];
  DBH[,1] = dbh0 ;
  for(y in 2:Y)
    DBH[,y] = DBH[,y-1] + gmax_i .* exp(-0.5 * square(log(DBH[,y-1] ./ dopt_i) ./ k_i)) ;
  for(n in 1:N)
    mu[n] = DBH[ind[n],year[n]] ;
}
model {
  dbh ~ normal(mu, sigma) ;
  gmax_i ~ normal(gmax[indsp], sigma_g) ;
  d_i ~ normal(d[indsp], sigma_d) ;
  k_i ~ normal(k[indsp], sigma_k) ;
  gmax ~ normal(0, 1) ;
  d ~ normal(0.4, 1) ;
  k ~ normal(0, 1) ;
  sigma ~ normal(0, 1) ;
  sigma_g ~ normal(0, 1) ;
  sigma_d ~ normal(0, 1) ;
  sigma_k ~ normal(0, 1) ;
}
