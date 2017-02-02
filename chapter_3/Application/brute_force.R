n_brute_force_reps <- 100
variables_to_save <- c("n_brute_force_reps")

# set.seed(72793333)


# Set up -----------------------------------------------------------------------

library(rstan)
library(edstan)
library(loo)
library(boot)
library(parallel)
library(doParallel)
library(foreach)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Borrow the simulation parameters from application.R
load("application.Rdata")
source("../functions.R")


# Set up data and models -------------------------------------------------------

# Person regression models to try
model_list <- list(~ 1,
                   ~ 1 + anger,
                   ~ 1 + male,
                   ~ 1 + anger + male,
                   ~ 1 + anger + male + I(anger*male))

dl <- irt_data(y = aggression$dich, jj = aggression$person,
               ii = aggression$item, covariates = aggression,
               formula = model_list[[4]])

sm <- stan_model("rasch_edstan_modified.stan")


# Fit the models and get conditional/marginal IC with bootstrapping ------------

cond_dic <- cond_waic <- list()
marg_dic <- marg_waic <- list()
max_rhat_brute_force <- list()

start <- Sys.time()
for(i in 1:n_brute_force_reps) {

  message(Sys.time(), ", replicate ", i)

  fit <- sampling(sm, dl, iter = n_iter, chains = n_chains, warmup = n_warmup)

  # Track convergence
  s <- summary(fit)[[1]]
  d <- data.frame(parameter = rownames(s), s)
  max_rhat_brute_force[[i]] <- d[d[, "Rhat"] == max(d[, "Rhat"]),]

  # Conditional IC
  ll_cond_ij <- cll(fit, dl, f_conditional_ij)
  cond_dic[[i]] <- unlist(dic(ll_cond_ij))
  cond_waic[[i]] <- waic_wrapper(ll_cond_ij$ll)

  # Marginal IC
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", n_agq)
  stopCluster(cl)
  marg_dic[[i]] <- unlist(dic(ll_marg))
  marg_waic[[i]] <- waic_wrapper(ll_marg$ll)

}
end <- Sys.time()
end - start

# Format results

combine_focus <- function(clist, mlist) {
  df_c <- data.frame(focus = "Conditional", model = 4,
                     replicate = 1:length(clist), do.call(rbind, clist))
  df_m <- data.frame(focus = "Marginal", model = 4,
                     replicate = 1:length(mlist), do.call(rbind, mlist))
  df_long <- melt(rbind(df_c, df_m), id.vars = c("focus", "model", "replicate"))
  return(df_long)
}

df_brute_force <- rbind(combine_focus(cond_dic, marg_dic),
                        combine_focus(cond_waic, marg_waic))

df_brute_force_max_rhat <- do.call(rbind, max_rhat_brute_force)


# Finish

variables_to_save <- c(variables_to_save, "df_brute_force",
                       "df_brute_force_max_rhat")

save(list = variables_to_save, file = "brute_force.Rdata")


max(df_brute_force_max_rhat$Rhat)
ggplot(subset(df_brute_force, variable %in% c("dic", "waic"))) +
  aes(x = value) +
  geom_histogram(bins = 10) +
  facet_wrap(~focus + variable, scales = "free_x")


