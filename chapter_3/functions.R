# Functions --------------------------------------------------------------------

# Replacement for rstan::get_posterior_means() that returns object with same
# structure as rstan::extract()

better_posterior_means <- function(draws) {
  f <- function(x) {
    if(is.matrix(x)) {
      matrix(apply(x, 2, mean), nrow = 1)
    } else {
      mean(x)
    }
  }
  lapply(draws, f)
}


# Replacement for rstan::extract() that takes a subset of iterations

extract_wrapper <- function(fit, iters = NULL, ...) {
  draws <- extract(fit, ...)
  if(is.null(iters)) {
    return(draws)
  } else {
    f <- function(x, iters) {
      if(is.matrix(x)) {
        x[iters,]
      } else {
        x[iters]
      }
    }
    lapply(draws, f, iters = iters)
  }
}


# Function to get conditional log-likelihoods at each iteration and for
# posterior means of parameters.

cll <- function(stan_fit = NULL, data_list, CFUN, draws = NULL, ll_par = NULL,
                iters = NULL, best_only = FALSE) {

  if(is.null(draws)) {
    draws <- extract_wrapper(stan_fit, iters)
  }

  n_iter <- length(draws$lp__)

  post_means <- better_posterior_means(draws)
  best_ll <- f_c(1, data_list = data_list, draws = post_means)

  if(best_only) {
    return(best_ll)
  } else {
    if(is.null(ll_par)) {
      ll <- t(sapply(1:n_iter, CFUN, data_list = data_list, draws = draws))
    } else {
      ll <- draws[[ll_par]]
    }
    return(list(ll = ll, best_ll = best_ll))
  }

}

# Function to calculate DIC

dic <- function(ll_obj) {
  full_ll <- apply(ll_obj$ll, 1, sum)
  full_best <- sum(ll_obj$best_ll)
  mean_lpd <-  mean(full_ll)
  pdic <- 2 * (full_best - mean_lpd)
  elpd_dic <- full_best - pdic
  c(elpd_dic = elpd_dic, p_dic = pdic, dic = -2*elpd_dic,
       best_lpd = full_best, mean_lpd = mean_lpd)
}


# Function to get marginal likelihoods with parallel processing.

mll_parallel <- function(stan_fit, data_list, MFUN, resid_name, sd_name, n_nodes,
                         iters = NULL, best_only = FALSE) {

  library(statmod)       # For gauss.quad.prob()
  library(matrixStats)   # For logSumExp()

  draws <- extract_wrapper(stan_fit, iters)
  n_iter <- ifelse(best_only, 0, length(draws$lp__))
  post_means <- better_posterior_means(draws)

  # Seperate out draws for residuals and their SD
  resid <- apply(draws[[resid_name]], 2, mean)
  stddev <- apply(draws[[resid_name]], 2, sd)

  # Get standard quadrature points
  std_quad <- gauss.quad.prob(n_nodes, "normal", mu = 0, sigma = 1)
  std_log_weights <- log(std_quad$weights)

  # Extra iteration is to evaluate marginal log-likelihood at parameter means.
  ll <- foreach(i = 1:(n_iter + 1), .combine = rbind,
                .packages = "matrixStats") %dopar% {

    ll_j <- matrix(NA, nrow = 1, ncol = ncol(draws[[resid_name]]))

    for(j in 1:ncol(ll_j)) {

      # Set up adaptive quadrature using SD for residuals either from draws or
      # posterior mean (for best_ll).
      sd_i <- ifelse(i <= n_iter, draws[[sd_name]][i], post_means[[sd_name]])
      adapt_nodes <- resid[j] + stddev[j] * std_quad$nodes
      log_weights <- log(sqrt(2*pi)) + log(stddev[j]) + std_quad$nodes^2/2 +
        dnorm(adapt_nodes, sd = sd_i, log = TRUE) + std_log_weights

      # Evaluate mll with adaptive quadrature. If at n_iter + 1, evaluate
      # marginal likelihood at posterior means.
      if(i <= n_iter) {
        loglik_by_node <- sapply(adapt_nodes, FUN = MFUN, r = j, iter = i,
                                 data_list = data_list, draws = draws)
        weighted_loglik_by_node <- loglik_by_node + log_weights
        ll_j[1,j] <- logSumExp(weighted_loglik_by_node)
      } else {
        loglik_by_node <- sapply(adapt_nodes, FUN = MFUN, r = j, iter = 1,
                                 data_list = data_list, draws = post_means)
        weighted_loglik_by_node <- loglik_by_node + log_weights
        ll_j[1,j] <- logSumExp(weighted_loglik_by_node)
      }

    }

    ll_j

  }

  if(best_only) {
    return(ll[nrow(ll), ])
  } else {
    return(list(ll = ll[-nrow(ll), ], best_ll = ll[nrow(ll), ]))
  }

}


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

