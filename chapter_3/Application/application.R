n_posterior <- 10000
n_warmup <- 500
n_chains <- 5
n_iter <- n_posterior / n_chains + n_warmup
n_cores <- parallel::detectCores()
n_repeats <- 10

# Increase each node by 50% over last, rounding to odd number
n_agq_try <- c(7, rep(NA, times = 3))
for(i in 2:length(n_agq_try)) {
  x <- floor(1.5*n_agq_try[i-1])
  is_even <- x %% 2 == 0
  n_agq_try[i] <- x + is_even
}

monitor_pars <- c("delta_free", "theta", "sigma", "lambda_adj", "lp__")

start <- Sys.time()

variables_to_save <- ls()

set.seed(727933533)


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

# Fit model 4
dl <- irt_data(y = aggression$dich, jj = aggression$person,
               ii = aggression$item, covariates = aggression,
               formula = model_list[[4]])
fit <- irt_stan(dl, "rasch_edstan_modified.stan", iter = n_iter,
                chains = n_chains, warmup = n_warmup)

# Try the node counts until difference < .01 acheived.
i <- 1
success <- FALSE
agq_results <- matrix(NA, ncol = 3, nrow = length(n_agq_try))
colnames(agq_results) <- c("dic", "waic", "looic")
while(i <= length(n_agq_try) & !success) {
  message(Sys.time(), ": Trying n_agq = ", n_agq_try[i])
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg_j <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", n_agq_try[i])
  stopCluster(cl)
  agq_results[i, "dic"] <- dic(ll_marg_j)["dic"]
  agq_results[i, "waic"] <- waic(ll_marg_j$ll)[["waic"]]
  agq_results[i, "looic"] <- loo(ll_marg_j$ll, cores = n_cores)[["looic"]]
  if(i > 1) success <- all(abs(agq_results[i] - agq_results[i]) < .01)
  i <- i + 1
}

if(success) {
  n_agq <- n_agq_try[i]
} else {
  warning("Search for number of adaptive quadrature nodes failed.")
}

variables_to_save <- c(variables_to_save, "n_agq", "agq_results")


# Fit the models and get conditional/marginal IC with bootstrapping ------------

# Function to get approximate MC error for WAIC penalty term
# This was added by request for paper with Sophia and Ed. It does not
# necessarily work and should not really be used.
waic_p_mcerror <- function(ll, n_chains) {
  a <- array(ll, c(nrow(ll) / n_chains, n_chains, ncol(ll)))
  m <- monitor(a, warmup = 0, print = FALSE)
  n <- m[,"n_eff"]
  v <- m[,"sd"]^2
  # mcerror_variances <- v / (2*(n-1))
  mcerror_variances <- 2/(n-1) * v^2
  total_penalty_mcerror <- sqrt(sum(mcerror_variances))
  return(total_penalty_mcerror)
}

# Function to get approximate MC error for PART of the DIC penalty term
# The MC error for the posterior mean log-likelihood
# This was added by request for paper with Sophia and Ed. It does not
# necessarily work and should not really be used.
dic_p_mcerror <- function(ll, n_chains) {
  full_data_loglik_by_draw <- apply(ll, 1, sum)
  a <- array(full_data_loglik_by_draw,
             c(length(full_data_loglik_by_draw) / n_chains, n_chains, 1))
  m <- monitor(a, warmup = 0, print = FALSE)
  return(m[1, "se_mean"])
}

cond_dic <- cond_waic <- cond_loo <- list()
marg_dic <- marg_waic <- marg_loo <- list()
cond_dic_p_mce <- marg_dic_p_mce <- cond_waic_p_mce <- marg_waic_p_mce <- list()
nonconverge <- list()

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
      message("Non-converged model")
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

  # MC error approximations for DIC and WAIC penalty terms
  cond_dic_p_mce[[i]] <- c(i = i, model = model_expand[i],
                           mce = dic_p_mcerror(ll_cond_ij$ll, n_chains = n_chains))
  marg_dic_p_mce[[i]] <- c(i = i, model = model_expand[i],
                           mce = dic_p_mcerror(ll_marg$ll, n_chains = n_chains))
  cond_waic_p_mce[[i]] <- c(i = i, model = model_expand[i],
                            mce = waic_p_mcerror(ll_cond_ij$ll, n_chains = n_chains))
  marg_waic_p_mce[[i]] <- c(i = i, model = model_expand[i],
                            mce = waic_p_mcerror(ll_marg$ll, n_chains = n_chains))

}

# Format IC results

combine_focus <- function(clist, mlist, id_vars = c("focus", "model", "i")) {
  df_c <- data.frame(focus = "Conditional", do.call(rbind, clist))
  df_m <- data.frame(focus = "Marginal", do.call(rbind, mlist))
  df_long <- melt(rbind(df_c, df_m), id.vars = id_vars)
  return(df_long)
}

df_models <- rbind(combine_focus(cond_dic, marg_dic),
                   combine_focus(cond_waic, marg_waic),
                   combine_focus(cond_loo, marg_loo))

id_vars <-  c("focus", "model", "i")
df_penalty_mce_dic <- combine_focus(cond_dic_p_mce, marg_dic_p_mce, id_vars)
df_penalty_mce_waic <- combine_focus(cond_waic_p_mce, marg_waic_p_mce, id_vars)

# Format summary() results

end <- Sys.time()
end - start

variables_to_save <- c(variables_to_save, "df_models", "bad_rhats",
                       "df_penalty_mce_dic", "df_penalty_mce_waic",
                       "start", "end")

save(list = variables_to_save, file = "application.Rdata")

