data {
  int<lower=1> I;                 // # obs per cluster
  int<lower=1> J;                 // # clusters
  int<lower=1, upper=I> ii[I*J];  // obs for n (1:I for each cluster)
  int<lower=1, upper=J> jj[I*J];  // cluster for n
  vector[I*J] y;                  // measurement for n
  int<lower=1> L;                 // # covariates
  matrix[J, L] X;                 // covariate matrix
}
transformed data{
  vector[I] y_vecs[J];
  for(n in 1:(I*J)) y_vecs[jj[n]][ii[n]] = y[n];
}
parameters {
  vector[L] beta;
  real<lower=0> sigma;
  real<lower=0> psi;
  vector[J] zeta;
}
model {
  vector[J] eta;
  eta = X*beta;
  beta ~ normal(0, 2);
  sigma ~ exponential(.1);
  psi ~ exponential(.1);
  zeta ~ normal(0, psi);
  y ~ normal(eta[jj] + zeta[jj], sigma);
}
generated quantities {
  vector[J] eta;
  vector[I*J] cll_ij;
  vector[J] mll_j;
  eta = X*beta;
  for(n in 1:I*J)
    cll_ij[n] = normal_lpdf(y[n] | eta[jj[n]] + zeta[jj[n]], sigma);
  {
    matrix[I, I] Omega;
    Omega = rep_matrix(psi^2, I, I);
    for(i in 1:I)
      Omega[i,i] = psi^2 + sigma^2;
    for(j in 1:J)
      mll_j[j] = multi_normal_lpdf(y_vecs[j] | rep_vector(eta[j], I), Omega);
  }
}
