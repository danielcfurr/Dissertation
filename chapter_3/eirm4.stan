data {
  int<lower=1> I;               // # questions
  int<lower=1> P;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // question for n
  int<lower=1, upper=P> pp[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
  int<lower=1> L;               // # item predictors
  matrix[I, L] X;               // item covariate matrix
  int<lower=1> K;               // # person predictors
  matrix[P, K] W;               // person covariate matrix
}
parameters {
  vector[L] beta;
  real<lower=0> tau;
  vector[I] epsilon;
  vector[K] gamma;
  real<lower=0> sigma;
  vector[P] zeta;
}
transformed parameters {
  vector[I] delta;
  vector[P] theta;
  delta <- X*beta + epsilon;
  theta <- W*gamma + zeta;
}
model {
  vector[N] eta;
  beta ~ normal(0, 2);
  gamma ~ normal(0, 2);
  tau ~ normal(0, 2);
  sigma ~ normal(0, 2);
  zeta ~ normal(0, sigma);
  epsilon ~ normal(0, tau);
  for (n in 1:N)
    eta[n] <- theta[pp[n]] - delta[ii[n]];
  y ~ bernoulli_logit(eta);
}
