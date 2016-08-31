# ------------------------------------------------------------------------------
# Simulate dataset using the random intercept model

rim_simulate <- function(I = 20, J = 100, sigma = 1, psi = 1,
                         beta = c(-.5, .5, .5, .5, -.5, 0), seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  return_list <- base_list <- list()
  base_list$I <- I
  base_list$J <- J
  base_list$ii <- rep(1:I, each = J)
  base_list$jj <- rep(1:J, times = I)
  for(i in 1:3) return_list[[i]] <- base_list
  names(return_list) <- paste0("stan_list_", 1:3)

  # Get covariates
  X <- matrix(rnorm(J*3, mean = 0, sd = 1), ncol = 3, nrow = J)
  X_rep <- matrix(rnorm(J*3, mean = 0, sd = 1), ncol = 3, nrow = J)
  colnames(X) <- colnames(X_rep) <- paste0("v", 1:3)
  X <- as.matrix(model.matrix(~ 1 + v1*v2 + v2*v3, data = as.data.frame(X)))
  X_rep <- as.matrix(model.matrix(~ 1 + v1*v2 + v2*v3, data = as.data.frame(X_rep)))

  # Get cluster means
  zeta <- rnorm(J, mean = X %*% beta, sd = psi)
  zeta_rep <- rnorm(J, mean = X_rep %*% beta, sd = psi)

  # Get level 1 residuals
  epsilon <- rnorm(I*J, mean = 0, sd = sigma)
  epsilon_rep <- rnorm(I*J, mean = 0, sd = sigma)

  # Get response variables
  y <- zeta[base_list$jj] + epsilon
  y_rep <- zeta_rep[base_list$jj] + epsilon_rep

  L <- 4:6
  for(i in 1:3) {
    return_list[[i]]$L <- L[i]
    return_list[[i]]$X <- X[, 1:L[i]]
    return_list[[i]]$X_rep <- X_rep[, 1:L[i]]
    return_list[[i]]$y <- y
    return_list[[i]]$y_rep <- y_rep
  }

  return_list[["df"]] <- as.data.frame(cbind(y = y,
                                             cluster = base_list$jj,
                                             X[base_list$jj, 2:4]))
  return_list[["df_rep"]] <- as.data.frame(cbind(y = y_rep,
                                                 cluster = base_list$jj,
                                                 X_rep[base_list$jj, 2:4]))

  return(return_list)

}


# ------------------------------------------------------------------------------
# Get marginal likelihoods using adaptive quadrature

mll_adapt <- function(data_list, stan_fit, n_nodes) {

  get_loglik_by_node <- function(node, y, eta, sigma) {
    sum(dnorm(y, mean = eta + node, sd = sigma, log = TRUE))
  }

  post_sum <- summary(stan_fit, pars = "zeta",  use_cache = FALSE)[[1]]
  post_draws <- extract(stan_fit)
  n_post <- length(post_draws$lp__)
  mll <- matrix(NA, nrow = n_post, ncol = data_list$J)
  
  std_quad <- gauss.quad.prob(n_nodes, "normal", mu = 0, sigma = 1) 
  
  for(i in 1:nrow(mll)) {
    for(j in 1:ncol(mll)) {
      
      adapt_nodes <- post_sum[j, "mean"] + post_sum[j, "sd"] * std_quad$nodes
      adapt_weights <- sqrt(2*pi) * post_sum[j, "sd"] * exp(std_quad$nodes^2/2) *
        dnorm(adapt_nodes, sd = post_draws$psi[i]) * std_quad$weights
      log_weights <- log(adapt_weights)
      
      loglik_by_node <- sapply(adapt_nodes,
                               FUN = get_loglik_by_node,
                               y = data_list$y[data_list$jj == j],
                               eta = post_draws$eta[i,j],
                               sigma = post_draws$sigma[i])
      weighted_loglik_by_node <- loglik_by_node + log_weights
      mll[i,j] <- logSumExp(weighted_loglik_by_node)
      
    }
  }

  return(mll)

}


# ------------------------------------------------------------------------------
# Get marginal likelihoods using adaptive quadrature: integrate()

mll_integrate <- function(data_list, stan_fit) {

  get_loglik_by_zeta <- function(zeta, y, eta, sigma, psi) {
    lik <- dnorm(y, mean = eta + zeta, sd = sigma, log = FALSE)
    prior <- dnorm(zeta, mean = 0, sd = psi, log = FALSE)
    return(prod(lik, prior))
  }

  post_draws <- extract(stan_fit)
  n_post <- length(post_draws$lp__)
  mll <- matrix(NA, nrow = n_post, ncol = data_list$J)

  for(i in 1:nrow(mll)) {

    for(j in 1:ncol(mll)) {
      integ <- integrate(Vectorize(get_loglik_by_zeta, vectorize.args = "zeta"),
                         -Inf, Inf,
                         y = data_list$y[data_list$jj == j],
                         eta = post_draws$eta[i,j],
                         sigma = post_draws$sigma[i],
                         psi = post_draws$psi[i])
      mll[i,j] <- log(integ$value)
    }

  }

  return(mll)

}

