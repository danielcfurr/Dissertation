library(statmod)       # For gauss.quad.prob()
library(matrixStats)   # For logSumExp()
library(rstan)
library(reshape2)
options(mc.cores = parallel::detectCores())
source("functions.R")
compile <- stan_model(file = "rim_mvn.stan")

# ------------------------------------------------------------------------------
# Simulate dataset using the random intercept model

rim_simulate <- function(I = 20, J = 100, sigma = 1, psi = 1,
                         beta = c(-.5, .5, .5, .5, -.5, 0),
                         link = "normal", seed = NULL) {

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
  colnames(X) <- paste0("v", 1:3)
  X <- as.matrix(model.matrix(~ 1 + v1*v2 + v2*v3, data = as.data.frame(X)))

  # Get cluster means
  zeta <- rnorm(J, mean = X %*% beta, sd = psi)

  # Get level 1 residuals
  epsilon <- rnorm(I*J, mean = 0, sd = sigma)

  # Get response variables
  if(link == "logit") {
    y <- rbinom(I*J, 1, boot::inv.logit(zeta[base_list$jj]))
  } else {
    y <- zeta[base_list$jj] + epsilon
  }

  L <- 4:6
  for(i in 1:3) {
    return_list[[i]]$L <- L[i]
    return_list[[i]]$X <- X[, 1:L[i]]
    return_list[[i]]$y <- y
  }

  return_list[["df"]] <- as.data.frame(cbind(y = y,
                                             cluster = base_list$jj,
                                             X[base_list$jj, 2:4]))

  return(return_list)

}



start <- Sys.time()

sim_nodes <- c(10, 20, 40, 80)
sim_I <- c(50, 100, 200)
sim_J <- 100
sim_sigma = 1
sim_psi = 1
sim_beta = c(-.5, .5, .5, .5, -.5, 0)
result_list <- list()

for(i in 1:length(sim_I)) {

  message(Sys.time(), " starting i = ", i)

  sim_list <- rim_simulate(sim_I[i], sim_J, sim_sigma, sim_psi, sim_beta)
  data_list <- sim_list[[2]]
  fit <- sampling(compile, data_list, chains = 5, iter = 200)

  # Get results from Stan fit itself
  post_draws <- extract(fit)
  mvn_df <- as.data.frame(post_draws$mll_j)
  mvn_df$method <- "MVN"
  mvn_df$I <- sim_I[i]
  mvn_df$iter <- 1:nrow(mvn_df)
  mvn_df$secs <- 0
  result_list[[length(result_list) + 1]] <- mvn_df

  # Get results from integrate/quadpack quadrature
  message(Sys.time(), " starting quadpack quadrature")
  secs <- system.time(
    integ <- mll_integrate(data_list, fit)
  )
  qpk_df <- as.data.frame(integ)
  qpk_df$method <- "QPK/AQ"
  qpk_df$I <- sim_I[i]
  qpk_df$iter <- 1:nrow(qpk_df)
  qpk_df$secs <- secs[3]
  result_list[[length(result_list) + 1]] <- qpk_df

  for(n in 1:length(sim_nodes)) {
    message(Sys.time(), " starting adaptive quadrature, n_nodes = ", sim_nodes[n])
    secs <- system.time(
      quad <- mll_adapt(fit, data_list, sim_nodes[n], f, "zeta", "psi"))
    )
    quad_df <- as.data.frame(quad)
    quad_df$method <- paste("Adapt Quad:", sim_nodes[n])
    quad_df$I <- sim_I[i]
    quad_df$iter <- 1:nrow(quad_df)
    quad_df$secs <- secs[3]
    result_list[[length(result_list) + 1]] <- quad_df
  }

}

end <- Sys.time()
seconds <- end - start

df <- do.call(rbind, result_list)
df <- melt(df, id.vars = c("method", "I", "iter", "secs"), value.name = "mll")
df$j <- as.numeric(sub("V", "", df$variable))
df$variable <- NULL

save(sim_nodes, sim_I, sim_J, sim_sigma, sim_psi, sim_beta, seconds, df,
     file = "sim_data.Rdata")
