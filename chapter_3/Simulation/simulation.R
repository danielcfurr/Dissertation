sim_nodes <- round(7*1.5^(0:5))
sim_I <- c(25, 50, 100)
sim_J <- 500
sim_sigma = 1
sim_psi = 1
sim_beta = c(-.5, .5, .5, .5, -.5, 0)

n_posterior <- 1000
n_warmup <- 250
n_chains <- 4
n_iter <- n_posterior / n_chains + n_warmup
n_cores <- parallel::detectCores()

variables_to_save <- ls()

set.seed(451680)


# Setup ------------------------------------------------------------------------

library(statmod)       # For gauss.quad.prob()
library(matrixStats)   # For logSumExp()
library(rstan)
library(loo)
library(doParallel)
library(foreach)
library(mvtnorm)
library(reshape2)

options(loo.cores = n_cores)
options(mc.cores = n_cores)

compile <- stan_model(file = "rim_mvn.stan")

source("../functions.R")


# Functions --------------------------------------------------------------------

# Function to simulate datasets using the random intercept model

rim_simulate <- function(I = 20, J = 100, sigma = 1, psi = 1,
                         beta = c(-.5, .5, .5, .5, -.5, 0),
                         link = "normal", seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  # Start Stan data list
  data_list <- list()
  data_list$I <- I
  data_list$J <- J
  data_list$ii <- rep(1:I, each = J)
  data_list$jj <- rep(1:J, times = I)

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
    y <- rbinom(I*J, 1, boot::inv.logit(zeta[data_list$jj]))
  } else {
    y <- zeta[data_list$jj] + epsilon
  }

  # Complete Stan data list
  data_list$L <- ncol(X)
  data_list$X <- X
  data_list$y <- y

  return(data_list)

}


# Function to obtain marginal likelihood with mll_parallel()

f_marginal <- function(node, r, iter, data_list, draws) {
  y <- data_list$y[data_list$jj == r]
  eta <- draws$eta[iter, r]
  zeta <- draws$eta[iter, r]
  sigma <- draws$sigma[iter]
  sum(dnorm(y, mean = eta + node, sd = sigma, log = TRUE))
}

variables_to_save <- c(variables_to_save, "rim_simulate", "f_marginal")


# Simulation -------------------------------------------------------------------

# result_list <- list()
dic_list <- list()
waic_list <- list()
loo_list <- list()
l <- 0

start <- Sys.time()
for(i in 1:length(sim_I)) {

  message(Sys.time(), " starting I = ", sim_I[i])

  data_list <- rim_simulate(sim_I[i], sim_J, sim_sigma, sim_psi, sim_beta)
  fit <- sampling(compile, data_list, chains = n_chains, iter = n_iter,
                  warmup = n_warmup)

  # Get results from Stan fit itself
  post_draws <- extract(fit)
  ll_obj <- list(ll = post_draws$mll_j, best_ll = numeric(sim_J))

  # Complications required for DIC - get marginal likelihood at posterior means
  post_means <- get_posterior_mean(fit)[, "mean-all chains"]
  beta_mean <- post_means[grepl("^beta\\[", names(post_means))]
  eta <- data_list$X %*% beta_mean
  sigma_mean <- post_means[grepl("^sigma$", names(post_means))]
  psi_mean <- post_means[grepl("^psi$", names(post_means))]
  Omega <- matrix(psi_mean^2, nrow = sim_I[i], ncol = sim_I[i]) +
    diag(sigma_mean^2, nrow = sim_I[i], ncol = sim_I[i])
  for(j in 1:sim_J) {
    ll_obj$best_ll[j] <- dmvnorm(data_list$y[data_list$jj == j],
                                 mean = rep(eta[j], times = sim_I[i]),
                                 sigma = Omega, log = TRUE)
  }

  # Store IC results
  l <- l + 1
  dic_list[[l]] <- c(nagq = 0, I = sim_I[i], dic(ll_obj))
  waic_list[[l]] <- c(nagq = 0, I = sim_I[i], waic_wrapper(ll_obj$ll))
  loo_list[[l]] <- c(nagq = 0, I = sim_I[i], loo_wrapper(ll_obj$ll))

  for(n in 1:length(sim_nodes)) {

    message("  ", Sys.time(), " starting AGQ with ", sim_nodes[n], " nodes")

    secs <- system.time({
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      quad <- mll_parallel(fit, data_list = data_list, MFUN = f_marginal,
                           resid_name = "zeta", sd_name = "psi",
                           n_nodes = sim_nodes[n])
      stopCluster(cl)
    })

    # Store IC results
    l <- l + 1
    dic_list[[l]] <- c(nagq = sim_nodes[n], I = sim_I[i], dic(quad))
    waic_list[[l]] <- c(nagq = sim_nodes[n], I = sim_I[i], waic_wrapper(quad$ll))
    loo_list[[l]] <- c(nagq = sim_nodes[n], I = sim_I[i], loo_wrapper(quad$ll))

  }

}

end <- Sys.time()
end - start


# Assemble simulation data and save --------------------------------------------

# Modify the lists from the simulation into data frames
dic_df <- melt(as.data.frame(do.call(rbind, dic_list)),
               id.vars = c("nagq", "I"))
dic_df$IC <- "DIC"
waic_df <- melt(as.data.frame(do.call(rbind, waic_list)),
                id.vars = c("nagq", "I"))
waic_df$IC <- "WAIC"
loo_df <- melt(as.data.frame(do.call(rbind, loo_list)),
               id.vars = c("nagq", "I"))
loo_df$IC <- "PSIS-LOO"

# Combine those data frames and clean up variable names
df_combine <- rbind(dic_df, waic_df, loo_df)
df_combine$variable <- gsub("_*(dic|waic|loo|looic)", "", df_combine$variable)
df_combine$variable[df_combine$variable == ""] <- "dev"

# Modify data frame so AGQ and MVN results are on same line
df <- subset(df_combine,
             variable %in% c("elpd", "p", "dev", "mean_lpd", "best_lpd"))
df_mvn <- subset(df, nagq == 0)
df_mvn$nagq <- NULL
names(df_mvn)[names(df_mvn) == "value"] <- "mvn"
names(df)[names(df) == "value"] <- "agq"
df <- merge(subset(df, nagq > 0), df_mvn)

# Sort data frame, get differences between AGQ and MVN, and N nodes
df <- df[with(df, order(I, IC, variable, nagq)), ]
df$dif_mvn <- df$agq - df$mvn
df$dif_nagq <- NA
df$dif_nagq[2:nrow(df)] <- df$agq[2:nrow(df)] - df$agq[2:nrow(df) - 1]
df$dif_nagq[df$nagq == min(sim_nodes)] <- NA

variables_to_save <- c(variables_to_save, "df", "df_combine", "start", "end")

save(list = variables_to_save, file = "simulation.Rdata")
