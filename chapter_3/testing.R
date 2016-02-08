library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)

sim_data <- function() {

  data_list <- list()
  data_list$I <- 32
  data_list$P <- 300
  data_list$N <- data_list$I*data_list$P
  data_list$ii <- rep(1:data_list$I, each = data_list$P)
  data_list$pp <- rep(1:data_list$P, times = data_list$I)

  item_vars <- expand.grid(0:1, 0:1, 0:1)
  item_vars <- item_vars[rep(1:nrow(item_vars), length.out = data_list$I), ]
  X <- model.matrix(~ 1 + Var1 + Var2 + Var3 + Var1:Var2, data = item_vars)
  beta <- c(-.5, .5, .5, .5, -.5)
  epsilon <- rnorm(nrow(X), 0, .3)
  delta <- X %*% beta + epsilon

  person_vars <- data.frame(Var1 = rep(0:1, length.out = data_list$P),
                            Var2 = rnorm(data_list$P, mean = 0, sd = 1))

  person_vars <- expand.grid(0:1, 0:1)
  person_vars <- person_vars[rep(1:nrow(person_vars), length.out = data_list$P), ]
  W <- model.matrix(~ 1 + Var1 + Var2, data = person_vars)
  W <- W[, -1]
  gamma <- rep(c(.5, .5), length.out = ncol(W))
  zeta <- rnorm(nrow(W), 0, 1)
  theta <- W %*% gamma + zeta

  eta <- theta[data_list$pp] - delta[data_list$ii]
  data_list$y <- as.numeric(runif(length(eta)) < plogis(eta))

  dl_2 <- data_list
  dl_2$L <- ncol(X)
  dl_2$X <- X
  dl_2$K <- ncol(W)
  dl_2$W <- W

  dl_1 <- dl_2
  X_1 <- model.matrix(~ 1 + Var1 + Var2 + Var3,
                      data = item_vars)
  dl_1$L <- ncol(X_1)
  dl_1$X <- X_1

  dl_3 <- dl_2
  X_3 <- model.matrix(~ 1 + Var1 + Var2 + Var3 + Var1:Var2 + Var2:Var3,
                      data = item_vars)
  dl_3$L <- ncol(X_3)
  dl_3$X <- X_3

  return(list(dl_1, dl_2, dl_3))

}

#stan_obj <- stan_model(file = "eirm.stan")
stan_obj <- stan_model(file = "lltme.stan")
ptm <- proc.time()
J <- 100
select <- matrix(NA, nrow = J, ncol = 2)
compare_1 <- list()
compare_2 <- list()
for(j in 1:J) {
  message("Starting ", j, " of ", J, ": ", Sys.time())
  d <- sim_data()
  loo1 <- list()
  loo2 <- list()
  for(m in 1:3) {
    fit <- sampling(stan_obj, data = d[[m]], chains = 4, iter = 500)
    ll <- extract(fit, pars = c("loglik_1", "loglik_2"))
    loo1[[m]] <- loo(ll$loglik_1, cores = 4)
    loo2[[m]] <- loo(ll$loglik_2, cores = 4)
    mon <- monitor(fit, print = FALSE)
    trouble <- mon[, "Rhat"] > 1.1
    message("  Model ", m, " complete: ", Sys.time(), ". ",
            sum(trouble), " Rhat > 1.1.")
    if(sum(trouble) > 0) message("!!! ", paste(rownames(mon)[trouble], collapse = " "))
  }
  compare_1[[j]] <- compare(loo1[[1]], loo1[[2]], loo1[[3]])
  select[j, 1] <- as.numeric(substr(row.names(compare_1[[j]])[1], 7, 7))
  compare_2[[j]] <- compare(loo2[[1]], loo2[[2]], loo2[[3]])
  select[j, 2] <- as.numeric(substr(row.names(compare_2[[j]])[1], 7, 7))
}
s <- proc.time() - ptm
s[3] / 60 / 60 * 100

table(select[,1])
table(select[,2])

save(compare_1, compare_2, select, file = "results.Rdata")




# -------------------

data(VerbAgg, package = "lme4")

X <- model.matrix(~ 1 + btype + situ + mode,
                  data = unique(VerbAgg[, c("item", "btype", "situ", "mode")]))
W_full <- model.matrix(~ 1 + Anger + Gender,
                       data = unique(VerbAgg[, c("id", "Anger", "Gender")]))
W <- W_full[, -1]
W[, "Anger"] <- W[, "Anger"] - mean(W[, "Anger"])

data_list <- list(I = nrow(X),
                  P = nrow(W),
                  N = nrow(VerbAgg),
                  ii = as.numeric(VerbAgg$item),
                  pp = as.numeric(VerbAgg$id),
                  y = as.numeric(VerbAgg$r2 == "Y"),
                  L = ncol(X),
                  X = X,
                  K = ncol(W),
                  W = W)

fit <- sampling(stan_obj, data = data_list, chains = 6, iter = 500,
                cores = 6)
print(fit, pars = c("beta", "tau", "delta", "gamma", "sigma"))
stan_rhat(fit, pars = "theta")

# -----

m <- monitor(fit, print = FALSE)
beta <- m[grep("^beta\\[.*]$", rownames(m)), "mean"]
tau <- m[grep("^tau$", rownames(m)), "mean"]
sim_delta <- data_list$X %*% beta + rnorm(data_list$I, 0, tau)
gamma <- m[grep("^gamma\\[.*]$", rownames(m)), "mean"]
sigma <- m[grep("^sigma$", rownames(m)), "mean"]
sim_theta <- data_list$W %*% gamma + rnorm(data_list$P, 0, sigma)
sim_eta <- sim_theta[data_list$pp] - sim_delta[data_list$ii]
sim_list <- data_list
sim_list$y <- as.numeric(runif(length(sim_eta)) < plogis(sim_eta))

fit_sim <- sampling(stan_obj, data = sim_list, chains = 6, iter = 500,
                    cores = 6)
print(fit_sim, pars = c("beta", "tau", "delta", "gamma", "sigma"))

# -----

sim_list2 <- data_list
sim_list2$X <- data_list$X[, c("(Intercept)", "situself", "modedo")]
sim_list2$X <- cbind(sim_list2$X,
                     var3 = rep(0:1, length.out = sim_list2$I),
                     interact = sim_list2$X[,"situself"] * sim_list2$X[,"modedo"])
beta2 <- c(beta[c(1, 4, 5)], 1.5, -.5)
sim_delta2 <- sim_list2$X %*% beta2 + rnorm(sim_list2$I, 0, tau)
sim_eta2 <- sim_theta[data_list$pp] - sim_delta2[data_list$ii]
sim_list2$y <- as.numeric(runif(length(sim_eta2)) < plogis(sim_eta2))

fit_sim2 <- sampling(stan_obj, data = sim_list2, chains = 6, iter = 500,
                     cores = 6)
print(fit_sim2, pars = c("beta", "tau", "delta", "gamma", "sigma"))

# -------------------

sim_list3 <- sim_list2
sim_list3$W[,"Anger"] <- as.numeric(sim_list3$W[,"Anger"] > 0)

fit_sim3 <- sampling(stan_obj, data = sim_list3, chains = 24, iter = 500,
                     cores = 6)
print(fit_sim3, pars = c("beta", "tau", "delta", "gamma", "sigma"))
stan_rhat(fit_sim3, pars = "theta")
monitor_sim_3 <- monitor(fit_sim3, print = FALSE)

stan_obj_alt <- stan_model(file = "eirm2.stan")
fit_sim_alt <- sampling(stan_obj_alt, data = sim_list3, chains = 24, iter = 500,
                        cores = 6)
print(fit_sim_alt, pars = c("beta", "tau", "delta", "gamma", "sigma"))
stan_rhat(fit_sim_alt, pars = "theta")
monitor_alt <- monitor(fit_sim_alt, print = FALSE)




# -------------------

d <- sim_data()[[1]]
wanted_stubs <- c("beta", "delta", "tau", "gamma", "theta", "sigma")
df_list <- list()

for(i in 1:4) {

  # Get .stan file name
  f <- paste0("eirm", i, ".stan")

  # Fit model and get time in seconds
  message("--- Started: ", f, " ---")
  fit <- stan(file = f, data = d, iter = 500, chains = 24, cores = 6)
  seconds <- sum(get_elapsed_time(fit))
  message("--- Run time in seconds: ", seconds, " ---")

  # Get results, extract parameter "stubs", prep filter
  m <- monitor(fit, print = FALSE)
  parameters <- rownames(m)
  existing_stubs <- gsub("\\[.*]$", "", parameters)
  keep = existing_stubs %in% wanted_stubs

  # Assemble result
  df_list[[i]] <- data.frame(stub = existing_stubs[keep],
                   parameter = parameters[keep],
                   n_eff = m[keep, "n_eff"] / seconds,
                   row.names = NULL)
  names(df_list[[i]])[3] <- paste0("Model_", i)

  print(fit, pars = c("beta", "tau", "gamma", "sigma"))
  message("--- Finished: ", f, " ---")

}

df <- merge(df_list[[1]], df_list[[2]])
df <- merge(df, df_list[[3]])
df <- merge(df, df_list[[4]])
stub_string <- as.character(df$stub)
stub_factor <- factor(stub_string, levels = wanted_stubs)
df$stub <- stub_factor


library(ggplot2)
ggplot(data = df) +
  aes(x = Model_1, y = Model_2) +
  facet_wrap(~ stub) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(data = df) +
  aes(x = Model_1, y = Model_3) +
  facet_wrap(~ stub) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(data = df) +
  aes(x = Model_1, y = Model_4) +
  facet_wrap(~ stub) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)



# -------------------

generated quantities {

  vector[I] delta_struct;
  vector[I] delta_resid;
  vector[P] theta_struct;
  vector[P] theta_resid;
  vector[N] loglik;
  vector[I] i_loglik1;
  vector[I] i_loglik2;
  vector[I] e_loglik;

  delta_struct <- X*beta;
  delta_resid <- delta - delta_struct;
  theta_struct <- W*gamma;
  theta_resid <- theta - theta_struct;

  i_loglik1 <- rep_vector(0, I);
  for (n in 1:N) {
    loglik[n] <- bernoulli_logit_log(y[n], theta[pp[n]] - delta[ii[n]]);
    i_loglik1[ii[n]] <- i_loglik1[ii[n]] + loglik[n];
  }
  for (i in 1:I) {
    e_loglik[i] <- normal_log(delta_resid[i], 0, tau);
    i_loglik2[i] <- i_loglik1[i] + e_loglik[i];
  }

}
