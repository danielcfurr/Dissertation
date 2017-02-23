sim_I <- c(25, 50, 100)
sim_J <- 200
sim_sigma = 1
sim_psi = 1
sim_beta = c(-.5, .5, .5, .5, -.5, 0)

# Increase each node by 50% over last, rounding to odd number
sim_nodes <- c(7, rep(NA, times = 7))
for(i in 2:length(sim_nodes)) {
  x <- floor(1.5*sim_nodes[i-1])
  is_even <- x %% 2 == 0
  sim_nodes[i] <- x + is_even
}

n_posterior <- 2000
n_warmup <- 500
n_chains <- 4
n_iter <- n_posterior / n_chains + n_warmup
n_cores <- parallel::detectCores()

monitor_pars = c("beta", "sigma", "psi", "zeta", "lp__", "mll_j", "cll_ij")

variables_to_save <- ls()

set.seed(451667)


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
dic_marg_list <- dic_cond_list <- list()
waic_marg_list <- waic_cond_list <- list()
loo_marg_list <- loo_cond_list <- list()
mll_mvn_list <- mll_agq_list <- list()
l <- 0
nonconverge <- list()

start <- Sys.time()
for(i in 1:length(sim_I)) {

  message(Sys.time(), " starting I = ", sim_I[i])

  data_list <- rim_simulate(sim_I[i], sim_J, sim_sigma, sim_psi, sim_beta)


  retry <- TRUE
  while(retry) {
    fit <- sampling(compile, data_list, chains = n_chains, iter = n_iter,
                    warmup = n_warmup)
    sum_mat <- summary(fit, pars = monitor_pars, probs = NA)[["summary"]]
    bad_rhats <- rownames(sum_mat)[sum_mat[, "Rhat"] > 1.1]
    if(length(bad_rhats) > 0) {
      nonconverge[[length(nonconverge) + 1]] <- row.names(sum_mat)[bad_rhats]
    } else {
      retry <- FALSE
    }
  }

  # Get results from Stan fit itself
  post_draws <- extract(fit)

  # Store marginal results for WAIC, LOO, and MLL matrix
  l <- l + 1
  waic_marg_list[[l]] <- c(nagq = 0, I = sim_I[i], waic_wrapper(post_draws$mll_j))
  loo_marg_list[[l]] <- c(nagq = 0, I = sim_I[i], loo_wrapper(post_draws$mll_j))
  mll_mvn_list[[i]] <- list(I = sim_I[i], mll = post_draws$mll_j)

  # Store marginal results for DIC
  post_means <- better_posterior_means(post_draws)
  Omega <- matrix(post_means$psi^2, nrow = sim_I[i], ncol = sim_I[i]) +
    diag(post_means$sigma^2, nrow = sim_I[i], ncol = sim_I[i])
  mll_obj <- list(ll = post_draws$mll_j, best_ll = numeric(sim_J))
  for(j in 1:sim_J) {
    mll_obj$best_ll[j] <- dmvnorm(data_list$y[data_list$jj == j],
                                  mean = rep(post_means$eta[j], times = sim_I[i]),
                                  sigma = Omega, log = TRUE)
  }
  dic_marg_list[[l]] <- c(nagq = 0, I = sim_I[i], dic(mll_obj))

  # Store conditional results for WAIC, LOO
  waic_cond_list[[i]] <- c(nagq = 0, I = sim_I[i], waic_wrapper(post_draws$cll_ij))
  loo_cond_list[[i]] <- c(nagq = 0, I = sim_I[i], loo_wrapper(post_draws$cll_ij))

  # Store conditional results for DIC
  cll_obj <- list(ll = post_draws$cll_ij, best_ll = numeric(sim_I[i]*sim_J))
  cll_obj$best_ll <- dnorm(data_list$y,
                  post_means$eta[data_list$jj] + post_means$zeta[data_list$jj],
                  post_means$sigma, log = TRUE)
  dic_cond_list[[i]] <- c(nagq = 0, I = sim_I[i], dic(cll_obj))


  for(n in 1:length(sim_nodes)) {

    message("  ", Sys.time(), " starting AGQ with ", sim_nodes[n], " nodes")

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    quad <- mll_parallel(fit, data_list = data_list, MFUN = f_marginal,
                         resid_name = "zeta", sd_name = "psi",
                         n_nodes = sim_nodes[n])
    stopCluster(cl)

    # Store IC results
    l <- l + 1
    dic_marg_list[[l]] <- c(nagq = sim_nodes[n], I = sim_I[i], dic(quad))
    waic_marg_list[[l]] <- c(nagq = sim_nodes[n], I = sim_I[i], waic_wrapper(quad$ll))
    loo_marg_list[[l]] <- c(nagq = sim_nodes[n], I = sim_I[i], loo_wrapper(quad$ll))
    mll_agq_list[[length(mll_agq_list) + 1]] <-
      list(nagq = sim_nodes[n], I = sim_I[i], mll = quad$ll)

  }

}

end <- Sys.time()
end - start

variables_to_save <- c(variables_to_save, "mll_mvn_list", "mll_agq_list",
                       "start", "end")


# Assemble IC quadrature comparisons -------------------------------------------

reformat_df <- function(list_obj, ic, focus) {
  df <- as.data.frame(do.call(rbind, list_obj))
  df_melt <- melt(df, id.vars = c("nagq", "I"))
  df_melt$IC <- ic
  df_melt$focus <- focus
  df_melt$variable <- gsub("_*(dic|waic|loo|looic)", "", df_melt$variable)
  df_melt$variable[df_melt$variable == ""] <- "dev"
  return(df_melt)
}

df_combine <- rbind(reformat_df(dic_marg_list, "DIC", "Marginal"),
                    reformat_df(waic_marg_list, "WAIC", "Marginal"),
                    reformat_df(loo_marg_list, "PSIS-LOO", "Marginal"),
                    reformat_df(dic_cond_list, "DIC", "Conditional"),
                    reformat_df(waic_cond_list, "WAIC", "Conditional"),
                    reformat_df(loo_cond_list, "PSIS-LOO", "Conditional"))

# Modify data frame so AGQ and MVN results are on same line
df_quad_compare <- subset(df_combine, focus == "Marginal" &
                            variable %in% c("elpd", "p", "dev", "mean_lpd", "best_lpd"))
df_mvn <- subset(df_quad_compare, nagq == 0)
df_mvn$nagq <- NULL
names(df_mvn)[names(df_mvn) == "value"] <- "mvn"
names(df_quad_compare)[names(df_quad_compare) == "value"] <- "agq"
df_quad_compare <- merge(subset(df_quad_compare, nagq > 0), df_mvn)

# Sort data frame, get differences between AGQ and MVN, and N nodes
df_quad_compare <- df_quad_compare[with(df_quad_compare,
                                        order(I, IC, variable, nagq)), ]
df_quad_compare$dif_mvn <- df_quad_compare$agq - df_quad_compare$mvn
df_quad_compare$dif_nagq <- NA
df_quad_compare$dif_nagq[2:nrow(df_quad_compare)] <-
  df_quad_compare$agq[2:nrow(df_quad_compare)] -
  df_quad_compare$agq[2:nrow(df_quad_compare) - 1]
df_quad_compare$dif_nagq[df_quad_compare$nagq == min(sim_nodes)] <- NA

variables_to_save <- c(variables_to_save, "df_combine", "df_quad_compare")

save(list = variables_to_save, file = "simulation.Rdata")



#

log_mean_prob_of_logprobs <- function(logs) {
  logSumExp(logs) + log(1/length(logs))
}

mvn_means <- mvn_vars <- list()
for(i in 1:length(mll_mvn_list)) {
  mvn_means[[i]] <- apply(mll_mvn_list[[i]][["mll"]], 2, log_mean_prob_of_logprobs)
  mvn_vars[[i]] <- apply(mll_mvn_list[[i]][["mll"]], 2, var)
}

r <- list()
for(i in 1:length(mll_agq_list)) {
  which_I <- which(mll_agq_list[[i]][["I"]] == sim_I)
  aqq_means <- apply(mll_agq_list[[i]][["mll"]], 2, log_mean_prob_of_logprobs)
  aqq_vars <- apply(mll_agq_list[[i]][["mll"]], 2, var)
  r[[i]] <- data.frame(j = 1:sim_J,
                       mvn_mean = mvn_means[[which_I]],
                       agq_mean = aqq_means,
                       mvn_var = mvn_vars[[which_I]],
                       agq_var = aqq_vars,
                       I = mll_agq_list[[i]][["I"]],
                       nagq = mll_agq_list[[i]][["nagq"]])
}
r <- do.call(rbind, r)
r$dif_mean <- r$agq_mean - r$mvn_mean
r$dif_var <- r$agq_var - r$mvn_var
ggplot(r) + aes(x = as.factor(nagq), y = dif_mean) + geom_boxplot() +
  facet_wrap(~as.factor(I), scales = "free_y")
ggplot(r) + aes(x = as.factor(nagq), y = dif_var) + geom_boxplot() +
  facet_wrap(~as.factor(I), scales = "free_y")

ggplot(r) + aes(x = as.factor(nagq), y = dif_mean) +
  geom_jitter(height = 0, alpha = .3) +
  facet_wrap(~as.factor(I), scales = "free_y")
ggplot(r) + aes(x = as.factor(nagq), y = dif_var) +
  geom_jitter(height = 0, alpha = .3) +
  facet_wrap(~as.factor(I), scales = "free_y")

r_long <- melt(r[, c("I", "nagq", "dif_mean", "dif_var")],
               id.vars = c("I", "nagq"))
ggplot(r_long) +
  aes(x = as.factor(nagq), y = value) +
  geom_jitter(height = 0, alpha = .3) +
  facet_grid(as.factor(I) ~ variable, scales = "free_y")
ggplot(subset(r_long, nagq > 20 & I == 100)) +
  aes(x = as.factor(nagq), y = value) +
  geom_jitter(height = 0, alpha = .3) +
  facet_grid(as.factor(I) ~ variable, scales = "free_y")


# variables_to_save <- c(variables_to_save, "df_quad_compare", "df_combine",
#                        "start", "end")
#
# save(list = variables_to_save, file = "simulation.Rdata")
