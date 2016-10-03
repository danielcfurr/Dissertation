data {
  int<lower=1> I;               // # items
  int<lower=1> J;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // item for n
  int<lower=1, upper=J> jj[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
  int<lower=1> K;               // # item predictors
  matrix[J, K] W;               // item covariate matrix
}
parameters {
  vector[K] gamma;
  vector[I] delta_unit;
  vector[J] zeta_unit;
  real<lower=0> tau;
  real<lower=0> sigma;
}
transformed parameters {
  vector[J] theta;
  vector[J] theta_fix;
  vector[J] zeta;
  vector[I] delta;
  zeta = zeta_unit * sigma;
  theta_fix = W*gamma;
  theta = theta_fix + zeta;
  delta = delta_unit * tau;
}
model {
  gamma ~ normal(0, 2);
  tau ~ exponential(.1);
  sigma ~ exponential(.1);
  zeta_unit ~ normal(0, 1);
  delta_unit ~ normal(0, 1);
  y ~ bernoulli_logit(theta[jj] - delta[ii]);
}
