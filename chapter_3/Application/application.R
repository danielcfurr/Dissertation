n_posterior <- 10000
n_warmup <- 500
n_chains <- 5
n_iter <- n_posterior / n_chains + n_warmup
n_cores <- parallel::detectCores()
n_repeats <- 10

# Increase each node by 50% over last, rounding to odd number
n_agq_try <- c(7, rep(NA, times = 2))
for(i in 2:length(n_agq_try)) {
  x <- floor(1.5*n_agq_try[i-1])
  is_even <- x %% 2 == 0
  n_agq_try[i] <- x + is_even
}

monitor_pars <- c("delta_free", "theta", "sigma", "lambda_adj", "lp__")

start <- Sys.time()

variables_to_save <- ls()

# set.seed(72793333)


# Set up -----------------------------------------------------------------------

library(rstan)
library(edstan)
library(loo)
library(boot)
library(doParallel)
library(foreach)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = n_cores)
options(loo.cores = n_cores)

source("../functions.R")


# Set up data and models -------------------------------------------------------

# Person regression models to try
model_list <- list(~ 1,
                   ~ 1 + anger,
                   ~ 1 + male,
                   ~ 1 + anger + male,
                   ~ 1 + anger + male + I(anger*male))


# Select number of adaptive quadrature points ----------------------------------

dl <- irt_data(y = aggression$dich, jj = aggression$person,
               ii = aggression$item, covariates = aggression,
               formula = model_list[[4]])
# fit <- irt_stan(dl, "rasch_edstan_modified.stan", iter = n_iter,
#                 chains = n_chains, warmup = n_warmup)
fit <- irt_stan(dl, "rasch_edstan_modified.stan", iter = n_iter,
                chains = n_chains, warmup = n_warmup)

# Get Marginal IC results for each choice of number AGQ nodes
agq_try <- matrix(NA, ncol = 4, nrow = length(n_agq_try))
colnames(agq_try) <- c("nodes", "dic", "waic", "looic")
agq_try[, "nodes"] <- n_agq_try
agq_change <- agq_try[-1, ]
for(i in 1:length(n_agq_try)) {
  message(Sys.time(), ": i=", i)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg_j <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", agq_try[i])
  stopCluster(cl)
  agq_try[i, "dic"] <- dic(ll_marg_j)["dic"]
  agq_try[i, "waic"] <- waic(ll_marg_j$ll)[["waic"]]
  agq_try[i, "looic"] <- loo(ll_marg_j$ll, cores = n_cores)[["looic"]]
  if(i > 1) agq_change[i-1, 2:4] <- agq_try[i, 2:4] - agq_try[i-1, 2:4]
}

# Select first choice of N nodes with < .01 change
acceptable_change <- apply(agq_change[, -1], 1, function(x) all(abs(x) < .01))
first_acceptable <- min(which(acceptable_change == TRUE))
(n_agq <- agq_try[first_acceptable + 1])

variables_to_save <- c(variables_to_save, "agq_try", "n_agq")


# Fit the models and get conditional/marginal IC with bootstrapping ------------

cond_dic <- cond_waic <- cond_loo <- list()
marg_dic <- marg_waic <- marg_loo <- list()
bad_rhats <- list()

model_expand <- rep(1:length(model_list), times = n_repeats)

for(i in 1:length(model_expand)) {

  message(Sys.time(), ", Step ", i, " of ", length(model_expand), ": Starting fit")

  dl <- irt_data(y = aggression$dich, jj = aggression$person,
                 ii = aggression$item, covariates = aggression,
                 formula = model_list[[model_expand[i]]])

  retry <- TRUE
  while(retry) {
    fit <- irt_stan(dl, "rasch_edstan_modified.stan", iter = n_iter,
                    chains = n_chains, warmup = n_warmup)
    sum_mat <- summary(fit, pars = monitor_pars, probs = NA)[["summary"]]
    bad_rhats <- rownames(sum_mat)[sum_mat[, "Rhat"] > 1.1]
    if(length(bad_rhats) > 0) {
      nonconverge[[length(nonconverge) + 1]] <- row.names(sum_mat)[bad_rhats]
    } else {
      retry <- FALSE
    }
  }

  # Conditional IC
  message(Sys.time(), ", Step ", i, ": Starting conditional IC")
  ll_cond_ij <- cll(fit, dl, f_conditional_ij)
  cond_dic[[i]] <- c(i = i, model = model_expand[i], unlist(dic(ll_cond_ij)))
  cond_waic[[i]] <- c(i = i, model = model_expand[i], waic_wrapper(ll_cond_ij$ll))
  cond_loo[[i]] <- c(i = i, model = model_expand[i],
                     loo_wrapper(ll_cond_ij$ll, cores = n_cores))

  # Marginal IC
  message(Sys.time(), ", Step ", i, ": Starting marginal IC")
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", n_agq)
  stopCluster(cl)
  marg_dic[[i]] <- c(i = i, model = model_expand[i], unlist(dic(ll_marg)))
  marg_waic[[i]] <- c(i = i, model = model_expand[i], waic_wrapper(ll_marg$ll))
  marg_loo[[i]] <- c(i = i, model = model_expand[i],
                     loo_wrapper(ll_marg$ll, cores = n_cores))

}

# Format IC results

combine_focus <- function(clist, mlist) {
  df_c <- data.frame(focus = "Conditional", do.call(rbind, clist))
  df_m <- data.frame(focus = "Marginal", do.call(rbind, mlist))
  df_long <- melt(rbind(df_c, df_m), id.vars = c("focus", "model", "i"))
  return(df_long)
}

df_models <- rbind(combine_focus(cond_dic, marg_dic),
                   combine_focus(cond_waic, marg_waic),
                   combine_focus(cond_loo, marg_loo))

# Format summary() results

end <- Sys.time()
end - start

variables_to_save <- c(variables_to_save, "df_models", "bad_rhats",
                       "start", "end")

save(list = variables_to_save, file = "application.Rdata")

