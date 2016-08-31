library(statmod)       # For gauss.quad.prob()
library(matrixStats)   # For logSumExp()
library(rstan)
library(reshape2)
options(mc.cores = parallel::detectCores())
source("functions.R")
compile <- stan_model(file = "rim_mvn.stan")

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
      quad <- mll_adapt(data_list, fit, sim_nodes[n])
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
