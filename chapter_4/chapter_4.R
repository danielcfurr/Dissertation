n_posterior <- 10000
n_warmup <- 500
n_chains <- 20
n_iter <- n_posterior / n_chains + n_warmup
n_agq <- 20
n_cores <- parallel::detectCores()

set.seed(72793333)


# Set up -----------------------------------------------------------------------

library(rstan)
library(edstan)
library(loo)
library(doParallel)
library(foreach)
options(mc.cores = n_cores)
options(loo.cores = n_cores)

source("../functions.R")

# Wrappers for IC functions

waic_wrapper <- function(ll) {
  w <- waic(ll)
  p_04 <- sum(w$pointwise[, "p_waic"] > .4)
  return(c(unlist(w)[1:6], p_04 = p_04))
}

loo_wrapper <- function(ll) {
  l <- loo(ll)
  pk_05 <- sum(l$pareto_k > 0.5)
  pk_10 <- sum(l$pareto_k > 1.0)
  return(c(unlist(l)[1:6], pk_05 = pk_05, pk_10 = pk_10))
}

# Likelihood functions

f_conditional_ij <- function(iter, data_list, draws) {
  theta_vec <- draws$theta[iter, data_list$jj]
  delta_vec <- draws$delta[iter, data_list$ii]
  p <- boot::inv.logit(theta_vec - delta_vec)
  dbinom(data_list$y, 1, p, log = TRUE)
}

f_conditional_j <- function(iter, data_list, draws) {
  theta_vec <- draws$theta[iter, data_list$jj]
  delta_vec <- draws$delta[iter, data_list$ii]
  p <- boot::inv.logit(theta_vec - delta_vec)
  ll_ij <- dbinom(data_list$y, 1, p, log = TRUE)
  ll_j <- tapply(ll_ij, data_list$jj, sum)
  return(ll_j)
}

f_marginal_j <- function(node, r, iter, data_list, draws) {
  y <- data_list$y[data_list$jj == r]
  theta_fix <- draws$theta_fix[iter, r]
  delta <- draws$delta[iter, data_list$ii[data_list$jj == r]]
  p <- boot::inv.logit(theta_fix + node - delta)
  sum(dbinom(y, 1, p, log = TRUE))
}


# Set up data and models -------------------------------------------------------

compiled_model <- stan_model("rasch_reg.stan")

start <- Sys.time()

# Make data frame for item covariates
covars <- unique(aggression[, c("person", "anger", "male")])
covars <- covars[sort(covars$person), -1]
covars$anger <- with(covars, (anger - mean(anger)) / sd(anger))

# Person regression models to try
model_list <- list(~ 1,
                   ~ 1 + anger,
                   ~ 1 + male,
                   ~ 1 + anger + male,
                   ~ 1 + anger + male + I(anger*male))

# Make incomplete data list for Stan. Okay to ignore warnings.
data_list <- irt_data(y = aggression$dich,
                      jj = aggression$person,
                      ii = aggression$item)


# Fit models and store results -------------------------------------------------

df <- data.frame(focus = rep("", times = length(model_list)*(n_chains+1)*3),
                 model = NA, chain = NA,
                 elpd_dic = NA, p_dic = NA, dic = NA,
                 best_lpd = NA, mean_lpd = NA,
                 elpd_waic = NA, p_waic = NA, waic = NA,
                 se_elpd_waic = NA, se_p_waic = NA, se_waic = NA, p_04 = NA,
                 elpd_loo = NA, p_loo = NA, loo = NA,
                 se_elpd_loo = NA, se_p_loo = NA, se_loo = NA,
                 pk_05 = NA, pk_10 = NA,
                 stringsAsFactors = FALSE)

summary_list <- list()

i <- 0
for(m in 1:length(model_list)) {

  message(Sys.time(), ", Model ", m, ": Starting fit")

  fit <- rasch_reg(model_list[[m]], data_list, covars, compiled_model,
                   iter = n_iter, chains = n_chains, warmup = n_warmup)
  summary_list[[m]] <- summary(fit, probs = c())[[1]]

  message(Sys.time(), ", Model ", m, ": Starting conditional")

  # Conditional IC not aggregated

  ll_cond_ij <- cll(fit, data_list, f_conditional_ij)
  dic_cond_ij <- unlist(dic(ll_cond_ij))
  waic_cond_ij <- waic_wrapper(ll_cond_ij$ll)
  loo_cond_ij <- loo_wrapper(ll_cond_ij$ll)
  i <- i + 1
  df[i, 1] <- "Conditional"
  df[i, -1] <- c(m, 0, dic_cond_ij, waic_cond_ij, loo_cond_ij)

  message(Sys.time(), ", Model ", m, ": Starting conditional per chain")

  for(cc in 1:n_chains) {
    start <- (cc - 1)*(n_posterior/n_chains) + 1
    end <- start + (n_posterior/n_chains) - 1
    ll <- ll_cond_ij$ll[start:end, ]
    best_ll <- cll(fit, data_list, f_conditional_ij, iters = start:end,
                   best_only = TRUE)
    ll_cond_c <- list(ll = ll, best_ll = best_ll)
    dic_cond_c <- unlist(dic(ll_cond_c))
    waic_cond_c <- waic_wrapper(ll_cond_c$ll)
    loo_cond_c <- loo_wrapper(ll_cond_c$ll)
    i <- i + 1
    df[i, 1] <- "Conditional"
    df[i, -1] <- c(m, cc, dic_cond_c, waic_cond_c, loo_cond_c)
  }

  # Conditional IC aggregated to each j

  message(Sys.time(), ", Model ", m, ": Starting aggregated conditional")

  ll_cond_j <- cll(fit, data_list, f_conditional_j)
  dic_cond_j <- unlist(dic(ll_cond_j))
  waic_cond_j <- waic_wrapper(ll_cond_j$ll)
  loo_cond_j <- loo_wrapper(ll_cond_j$ll)
  i <- i + 1
  df[i, 1] <- "Aggregated conditional"
  df[i, -1] <- c(m, 0, dic_cond_j, waic_cond_j, loo_cond_j)

  message(Sys.time(), ", Model ", m, ": Starting aggregated conditional per chain")

  for(cc in 1:n_chains) {
    start <- (cc - 1)*(n_posterior/n_chains) + 1
    end <- start + (n_posterior/n_chains) - 1
    ll <- ll_cond_j$ll[start:end, ]
    best_ll <- cll(fit, data_list, f_conditional_j, iters = start:end,
                   best_only = TRUE)
    ll_cond_c <- list(ll = ll, best_ll = best_ll)
    dic_cond_c <- unlist(dic(ll_cond_c))
    waic_cond_c <- waic_wrapper(ll_cond_c$ll)
    loo_cond_c <- loo_wrapper(ll_cond_c$ll)
    i <- i + 1
    df[i, 1] <- "Aggregated conditional"
    df[i, -1] <- c(m, cc, dic_cond_c, waic_cond_c, loo_cond_c)
  }

  # Marginal IC aggregated to each j

  message(Sys.time(), ", Model ", m, ": Starting marginal IC")

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg <- mll_parallel(fit, data_list, f_marginal_j, "zeta", "sigma", n_agq)
  stopCluster(cl)
  dic_marg <- unlist(dic(ll_marg))
  waic_marg <- waic_wrapper(ll_marg$ll)
  loo_marg <- loo_wrapper(ll_marg$ll)
  i <- i + 1
  df[i, 1] <- "Marginal"
  df[i, -1] <- c(m, 0, dic_marg, waic_marg, loo_marg)

  message(Sys.time(), ", Model ", m, ": Starting marginal IC per chain")

  for(cc in 1:n_chains) {
    start <- (cc - 1)*(n_posterior/n_chains) + 1
    end <- start + (n_posterior/n_chains) - 1
    ll <- ll_marg$ll[start:end, ]
    cl <- makeCluster(1)
    registerDoParallel(cl)
    best_ll <- mll_parallel(fit, data_list, f_m, "zeta", "sigma", n_agq,
                            iters = start:end, best_only = TRUE)
    stopCluster(cl)
    ll_marg_c <- list(ll = ll, best_ll = best_ll)
    dic_marg_c <- unlist(dic(ll_marg_c))
    waic_marg_c <- waic_wrapper(ll_marg_c$ll)
    loo_marg_c <- loo_wrapper(ll_marg_c$ll)
    i <- i + 1
    df[i, 1] <- "Marginal"
    df[i, -1] <- c(m, cc, dic_marg_c, waic_marg_c, loo_marg_c)
  }

}

end <- Sys.time()
end - start

save(n_posterior, n_warmup, n_chains, n_iter, n_agq, n_cores, df, summary_list,
     covars, data_list, model_list, start, end, file = "chapter_4.Rdata")


# ----------

library(reshape2)
library(ggplot2)

sum_mat <- do.call(rbind, summary_list)
max_rhat <- max(sum_mat[, "Rhat"])
max_rhat

df_long <- melt(df, id.vars = c("model", "chain", "focus"))
df_mcerror <- aggregate(value ~ model + focus + variable,
                        data = subset(df_long, chain > 0),
                        FUN = function(x) sd(x)/sqrt(length(x)))
names(df_mcerror)[names(df_mcerror) == "value"] <- "mcerror"
df_r <- merge(subset(df_long, chain == 0), df_mcerror)
df_r$chain <- NULL

df_plot <- subset(df_r, variable %in% c("dic", "waic", "loo"))
df_plot$IC <- ic <- sub("^.*_", "", toupper(df_plot$variable))
df_plot$IC <- factor(df_plot$IC, c("LPD", "DIC", "WAIC", "LOO"))
df_plot$Model <- as.factor(df_plot$model)
ggplot(df_plot) +
  aes(x = Model, y = value, ymin = value - mcerror, ymax = value + mcerror,
      color = IC, pch = IC) +
  geom_pointrange(position = position_dodge(width = .5)) +
  facet_wrap(~focus, scales = "free_y") +
  labs(y = NULL, colour = NULL, pch = NULL)

df_tab <- subset(df_r, variable %in% c("p_dic", "p_waic", "p_loo"))
tab <- dcast(df_tab, focus + variable ~ model, value.var = "value")
knitr::kable(tab, format = "markdown", digits = 1)

df_ratio <- subset(df_r, variable %in% c("se_waic", "se_loo"))
names(df_ratio)[names(df_ratio) == "value"] <- "se"
df_ratio$variable <- sub("^se_", "", df_ratio$variable)
df_ratio$mcerror <- NULL
df_ratio <- merge(df_plot, df_ratio)
df_ratio$ratio <- with(df_ratio, mcerror / se)
max_ratio <- aggregate(ratio ~ focus + variable, data = df_ratio, FUN = max)
max_ratio

