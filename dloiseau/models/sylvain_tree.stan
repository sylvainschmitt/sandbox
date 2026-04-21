data {
  int<lower=0> N; // # obs
  int<lower=0> Y; // # years
  real dbh0; // diameter at recruitment
  real dmax; // species maximum diameter
  vector[N] dbh; // tree diameters
  array[N] int<lower=1> year ; // time since recruitment
}
parameters {
  real<lower=0.001, upper=10> gmax ;
  real<lower=0.1, upper=1> d;
  real<lower=0.001, upper=3> k;
  real<lower=0> sigma;
}
transformed parameters {
  vector[Y] DBH;
  vector[N] mu;
  real dopt=d*dmax;
  DBH[1] = dbh0;
  for(y in 2:Y)
    DBH[y] = DBH[y-1]+gmax*exp(-0.5 * square(log(DBH[y-1]/dopt)/k));
  for(n in 1:N)
    mu[n] = DBH[year[n]];
}
model {
  dbh ~ normal(mu, sigma) ;
  gmax ~ normal(0, 1) ;
  d ~ normal(0.4, 1) ;
  k ~ normal(0, 1) ;
}
