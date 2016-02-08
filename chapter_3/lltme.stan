data {
  int<lower=1> I;               // # questions
  int<lower=1> P;               // # persons
  int<lower=1> N;               // # observations
  int<lower=1, upper=I> ii[N];  // question for n
  int<lower=1, upper=P> pp[N];  // person for n
  int<lower=0, upper=1> y[N];   // correctness for n
  int<lower=1> L;               // # item predictors
  matrix[I, L] X;               // item covariate matrix
}
parameters {
  vector[L] beta;
  real<lower=0> tau;
  vector[I] delta;
  real<lower=0> sigma;
  vector[P] theta;
}
model {
  vector[N] eta;
  beta ~ normal(0, 2);
  tau ~ normal(0, 2);
  sigma ~ normal(0, 2);
  theta ~ normal(0, sigma);
  delta ~ normal(X*beta, tau);
  for (n in 1:N)
    eta[n] <- theta[pp[n]] - delta[ii[n]];
  y ~ bernoulli_logit(eta);
}
generated quantities {
  vector[I] delta_star;
  vector[I] epsilon_star;
  vector[N] loglik_1;
  vector[N] loglik_2;
  for (i in 1:I) epsilon_star[i] <- normal_rng(0, tau);
  delta_star <- X*beta + epsilon_star;
  for (n in 1:N) {
    loglik_1[n] <- bernoulli_logit_log(y[n], theta[pp[n]] - delta[ii[n]]);
    loglik_2[n] <- bernoulli_logit_log(y[n], theta[pp[n]] - delta_star[ii[n]]);
  }
}
