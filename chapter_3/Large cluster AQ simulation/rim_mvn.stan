data {
  int<lower=1> I;                 // # obs per cluster
  int<lower=1> J;                 // # clusters
  int<lower=1, upper=I> ii[I*J];  // obs for n (1:I for each cluster)
  int<lower=1, upper=J> jj[I*J];  // cluster for n
  vector[I*J] y;                  // measurement for n
  int<lower=1> L;                 // # covariates
  matrix[J, L] X;                 // covariate matrix
  vector[I*J] y_rep;                  // measurement for n
  matrix[J, L] X_rep;                 // covariate matrix
}
transformed data{
  vector[I] y_vecs[J];
  vector[I] y_vecs_rep[J];
  for(n in 1:(I*J)) y_vecs[jj[n]][ii[n]] <- y[n];
  for(n in 1:(I*J)) y_vecs_rep[jj[n]][ii[n]] <- y_rep[n];
}
parameters {
  vector[L] beta;
  real<lower=0> sigma;
  real<lower=0> psi;
  vector[J] zeta;
}
model {
  vector[J] eta;
  eta <- X*beta;
  sigma ~ exponential(.1);
  psi ~ exponential(.1);
  zeta ~ normal(0, psi);
  y ~ normal(eta[jj] + zeta[jj], sigma);
}
generated quantities {

  vector[J] mll_j;
  real mll;
  vector[J] eta;

  vector[J] mll_j_rep;
  real mll_rep;

  {

    matrix[I, I] Omega;
    vector[J] eta_rep;

    Omega <- rep_matrix(psi^2, I, I);
    for(i in 1:I)
      Omega[i,i] <- psi^2 + sigma^2;

    eta <- X*beta;
    for(j in 1:J)
      mll_j[j] <- multi_normal_log(y_vecs[j], rep_vector(eta[j], I), Omega);
    mll <- sum(mll_j);

    eta_rep <- X_rep*beta;
    for(j in 1:J)
      mll_j_rep[j] <- multi_normal_log(y_vecs_rep[j], rep_vector(eta_rep[j], I),
                                       Omega);
    mll_rep <- sum(mll_j_rep);

  }

}
