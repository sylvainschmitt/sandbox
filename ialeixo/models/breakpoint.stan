data {
  int<lower=1> N;  // # observations
  vector[N] t;  // leaf age
  vector[N] p; // photosynthetic efficiency
}
parameters {
  real a; // intercept
  real<lower=0> b1; // first slope
  real<upper=0> b2; // second slope
  real<lower=min(t), upper=max(t)> tau; // breakpoint
  real<lower=0> sigma; // residual error
}
transformed parameters {
  vector[N] mu;
  for (n in 1:N) {
    if (t[n] < tau) {
      mu[n] = a+b1*t[n];
    } else {
      mu[n] = a+b1*tau+b2*(t[n]-tau);
    }
  }
}

model {
  p ~ normal(mu, sigma);
}
generated quantities {
  real tmax = tau;
  real pmax = a+b1*tau;
}
