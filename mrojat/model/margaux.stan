data{
  int<lower=1> N;
  vector[N] y;
}
parameters{
  real lambda;
}
model{
  (y+1) ~ exponential(exp(lambda));
}
