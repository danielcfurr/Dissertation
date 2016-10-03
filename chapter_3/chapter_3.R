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

opts_list <- ls()

set.seed(451680)


# Setup ------------------------------------------------------------------------

library(statmod)       # For gauss.quad.prob()
library(matrixStats)   # For logSumExp()
library(rstan)
library(loo)
library(doParallel)
library(foreach)

options(loo.cores = n_cores)
options(mc.cores = n_cores)

compile <- stan_model(file = "rim_mvn.stan")

source("functions.R")


# Functions --------------------------------------------------------------------

# Function to simulate datasets using the random intercept model

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


# Function to obtain marginal likelihood with mll_parallel()

f_marginal <- function(node, r, iter, data_list, draws) {
  y <- data_list$y[data_list$jj == r]
  eta <- draws$eta[iter, r]
  zeta <- draws$eta[iter, r]
  sigma <- draws$sigma[iter]
  sum(dnorm(y, mean = eta + node, sd = sigma, log = TRUE))
}


# Simulation -------------------------------------------------------------------

result_list <- list()

start <- Sys.time()
for(i in 1:length(sim_I)) {

  message(Sys.time(), " starting I = ", sim_I[i])

  sim_list <- rim_simulate(sim_I[i], sim_J, sim_sigma, sim_psi, sim_beta)
  data_list <- sim_list[[2]]
  fit <- sampling(compile, data_list, chains = n_chains, iter = n_iter,
                  warmup = n_warmup)

  # Get results from Stan fit itself
  post_draws <- extract(fit)
  mvn_df <- as.data.frame(post_draws$mll_j)
  mvn_df$method <- "MVN"
  mvn_df$I <- sim_I[i]
  mvn_df$iter <- 1:nrow(mvn_df)
  mvn_df$secs <- 0
  result_list[[length(result_list) + 1]] <- mvn_df

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

    quad_df <- as.data.frame(quad$ll)
    quad_df$method <- paste("Adapt Quad:", sim_nodes[n])
    quad_df$I <- sim_I[i]
    quad_df$iter <- 1:nrow(quad_df)
    quad_df$secs <- secs[3]
    result_list[[length(result_list) + 1]] <- quad_df

  }

}

df_wide <- do.call(rbind, result_list)

waic_list <- list()
loo_list <- list()
idx <- unique(df_wide[, c("method", "I", "secs")])
for(i in 1:nrow(idx)) {
  keep <- df_wide$method == idx$method[i] & df_wide$I == idx$I[i]
  ll <- as.matrix(df_wide[keep, 1:sim_J])
  waic_list[[i]] <- waic_wrapper(ll)
  loo_list[[i]] <- loo_wrapper(ll)
}

df_waic <- cbind(idx, do.call(rbind, waic_list))
df_loo <- cbind(idx, do.call(rbind, loo_list))

end <- Sys.time()

end - start

save(df_waic, df_loo, start, end, list = opts_list, file = "chapter_3.Rdata")

# df_long <- melt(df_wide, id.vars = c("method", "I", "iter", "secs"), value.name = "mll")
# df_long$j <- as.numeric(sub("V", "", df_long$variable))
# df_long$variable <- NULL
# df_mvn <- subset(df_long, method == "MVN")
# df_mvn$method <- df_mvn$secs <- NULL
# names(df_mvn)[names(df_mvn) == "mll"] <- "mvn"
# df_agq <- subset(df_long, method != "MVN")
# df_agq$nodes <- as.numeric(gsub("^.*:", "", df_agq$method))
# df_agq$method <- NULL
# names(df_agq)[names(df_agq) == "mll"] <- "agq"
# df_long <- merge(df_mvn, df_agq)
# df_long$dif <- df_long$agq - df_long$mvn

# #
#
# ggplot(df) +
#   aes(x = as.factor(I), y = dif, fill = as.factor(nodes)) +
#   geom_boxplot() +
#   facet_wrap(~I, scales = "free_y")
#
# #
#
# pct_threshold <- function(x, t) sum(abs(x) > t) / length(x) * 100
# df_agg <- aggregate(dif ~ I + nodes, df, pct_threshold, t = .001)
#
# ggplot(df_agg) +
#   aes(x = as.factor(nodes), y = dif, fill = as.factor(nodes)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~I)
#
# #

# mean_and_var <- function(x) c(mean = mean(x), var = var(x))
# df_agg <- aggregate(cbind(mvn, agq) ~ j + I + nodes, df, mean_and_var)
# df_agg$dif_mean <- df_agg$agq[, "mean"] - df_agg$mvn[, "mean"]
# df_agg$dif_var <- df_agg$agq[, "var"] - df_agg$mvn[, "var"]
# pct_threshold <- function(x, t) sum(abs(x) > t) / length(x) * 100
# df_pct <- aggregate(cbind(dif_mean, dif_var) ~ I + nodes, df_agg,
#                     pct_threshold, t = .00001)
#
# # Dif in mean of posterior marginal likelihood
# ggplot(df_pct) +
#   aes(x = as.factor(nodes), y = dif_mean, fill = as.factor(nodes)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   facet_wrap(~I)
#
# # Dif in variance of posterior marginal likelihood
# ggplot(df_pct) +
#   aes(x = as.factor(nodes), y = dif_var, fill = as.factor(nodes)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   facet_wrap(~I)
#
# #




# Quick demo -------------------------------------------------------------------

# data_list <- rim_simulate(J = 10)$stan_list_2
# fit <- sampling(compile, data = data_list, chains = 4, iter = 100)
#
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
# m <- mll_parallel(fit, data_list, f_marginal, "zeta", "psi", 100)
# stopCluster(cl)
#
# ex <- extract(fit)
# ex$mll_j[1:3, 1:5]
# m$ll[1:3, 1:5]





