sim_bs_r <- 500
sim_bs_l <- rep(1:100, times = 2)
sim_bs_I <- 25

variables_to_save <- ls()

# set.seed(451680)
source("../functions.R")
load("simulation.Rdata")


# Setup ------------------------------------------------------------------------

library(statmod)       # For gauss.quad.prob()
library(matrixStats)   # For logSumExp()
library(rstan)
library(loo)
library(doParallel)
library(foreach)
library(mvtnorm)
library(reshape2)
library(boot)

options(loo.cores = 1)
options(mc.cores = n_cores)

compile <- stan_model(file = "rim_mvn.stan")

data_list <- rim_simulate(sim_bs_I, sim_J, sim_sigma, sim_psi, sim_beta)

waic_wrapper <- function(ll) {
  w <- waic(ll)
  c(lpd_waic = w$elpd_waic + w$p_waic, p_waic = w$p_waic,
    elpd_waic = w$elpd_waic, waic = w$waic)
}


# Simulation on bootstrap ------------------------------------------------------

waic_marg_list <- waic_cond_list <- vector("list", length(sim_bs_l))
nonconverge <- list()

sim_bs_start <- Sys.time()
for(i in 1:length(sim_bs_l)) {

  message(Sys.time(), " block sizes; i = ", i, ", l = ", sim_bs_l[i])

  retry <- TRUE
  while(retry) {
    fit <- sampling(compile, data_list, chains = n_chains, iter = n_iter,
                    warmup = n_warmup)
    sum_mat <- summary(fit, pars = monitor_pars, probs = NA)[["summary"]]
    bad_rhats <- sum_mat[, "Rhat"] > 1.1
    if(sum(bad_rhats) > 0) {
      nonconverge[[length(nonconverge) + 1]] <- row.names(sum_mat)[bad_rhats]
    } else {
      retry <- FALSE
    }
  }

  post_draws <- extract(fit)
  m <- summary(fit, pars = "mll_j", probs = NA)
  hist(m$summary[,"Rhat"])
  m <- summary(fit, pars = "cll_ij", probs = NA)
  hist(m$summary[,"Rhat"])

  # Marginal WAIC. bs$t0 is result for full data. mcerror is boostrap MC error.
  cl <- makeCluster(n_cores, type='PSOCK')
  clusterExport(cl, "waic")
  bs <- tsboot(post_draws$mll_j, statistic = waic_wrapper, l = sim_bs_l[i],
               R = sim_bs_r, sim = "fixed", endcorr = TRUE,
               cl = cl, parallel = "snow", ncpus = n_cores)
  stopCluster(cl)
  mcerror <- apply(bs$t, 2, sd)
  names(mcerror) <- paste0("mce_", names(bs$t0))
  waic_marg_list[[i]] = c(l = sim_bs_l[i], bs$t0, mcerror)

  # Conditional WAIC
  cl <- makeCluster(n_cores, type='PSOCK')
  clusterExport(cl, "waic")
  bs <- tsboot(post_draws$cll_ij, statistic = waic_wrapper, l = sim_bs_l[i],
               R = sim_bs_r, sim = "fixed", endcorr = TRUE,
               cl = cl, parallel = "snow", ncpus = n_cores)
  stopCluster(cl)
  mcerror <- apply(bs$t, 2, sd)
  names(mcerror) <- paste0("mce_", names(bs$t0))
  waic_cond_list[[i]] = c(l = sim_bs_l[i], bs$t0, mcerror)

  # Marginal DIC
#   post_means <- get_posterior_mean(fit)[, "mean-all chains"]
#   beta_mean <- post_means[grepl("^beta\\[", names(post_means))]
#   eta <- data_list$X %*% beta_mean
#   sigma_mean <- post_means[grepl("^sigma$", names(post_means))]
#   psi_mean <- post_means[grepl("^psi$", names(post_means))]
#   Omega <- matrix(psi_mean^2, nrow = sim_bs_I, ncol = sim_bs_I) +
#     diag(sigma_mean^2, nrow = sim_bs_I, ncol = sim_bs_I)
#   mll_obj <- list(ll = post_draws$mll_j, best_ll = numeric(sim_J))
#   for(j in 1:sim_J) {
#     mll_obj$best_ll[j] <- dmvnorm(data_list$y[data_list$jj == j],
#                                  mean = rep(eta[j], times = sim_bs_I),
#                                  sigma = Omega, log = TRUE)
#   }
#
#   colnames(mll_obj$ll) <- paste0("loglik[", 1:ncol(mll_obj$ll), "]")
#   mat <- cbind(as.matrix(fit), mll_obj$ll)
#
#   dic_cond_ij_wrapper(mat, data_list)
#
#   cl <- makeCluster(n_cores, type='PSOCK')
#   clusterExport(cl, c("dic", "cll", "f_conditional_ij",
#                       "better_posterior_means"))
#   bs <- tsboot(mat, statistic = dic_cond_ij_wrapper, l = block_sizes[b],
#                R = n_bs_r, sim = tsboot_sim, data_list = data_list,
#                cl = cl, parallel = "snow", ncpus = n_cores)
#   stopCluster(cl)
#
# }
#
# dic_marg_ij_wrapper <- function(mat, data_list) {
#   ll <- mat[, grepl("^loglik\\[.*]$", colnames(mat))]
#   eta_draws <- mat[, grepl("^eta\\[.*]$", colnames(mat))]
#   eta_means <- apply(eta_draws, 2, mean)
#
}

sim_bs_end <- Sys.time()
sim_bs_end - sim_bs_start

df_waic_marg <- data.frame(do.call(rbind, waic_marg_list),
                           focus = "Marginal", IC = "WAIC")
df_waic_cond <- data.frame(do.call(rbind, waic_cond_list),
                           focus = "Conditional", IC = "WAIC")
df_bootstrap <- rbind(df_waic_marg, df_waic_cond)

variables_to_save <- c(variables_to_save, "sim_bs_start", "sim_bs_end",
                       "nonconverge", "df_bootstrap")

save(list = variables_to_save, file = "bootstrap.Rdata")


#

est_names <- c("lpd_waic", "p_waic", "waic")
mce_names <- paste0("mce_", est_names)

# Brute force results histograms
df_bootstrap$i <- 1:nrow(df_bootstrap)
df_bootstrap_long <- melt(df_bootstrap, id.vars = c("focus", "IC", "l", "i"))
ggplot(subset(df_bootstrap_long, variable %in% est_names)) +
  aes(x = value) +
  geom_histogram(bins = 10) +
  facet_wrap(~focus + variable, nrow = 2, scales = "free_x") #+
  #my_theme

# Boostrap results with brute force sd for comparison
df_brute_force <- aggregate(value ~ focus + IC + variable,
                            subset(df_bootstrap_long, variable %in% est_names),
                            FUN = sd)
names(df_brute_force)[names(df_brute_force) == "value"] <- "brute_force"
df_bs_mcerror <- subset(df_bootstrap_long, variable %in% mce_names)
df_bs_mcerror$variable <- sub("mce_", "", df_bs_mcerror$variable)
names(df_bs_mcerror)[names(df_bs_mcerror) == "value"] <- "bootstrap"
df_combine <- merge(df_bs_mcerror, df_brute_force)
ggplot(df_combine) +
  aes(x = l, y = bootstrap) +
  geom_hline(aes(yintercept = brute_force), linetype = "dashed") +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~focus + variable, scales = "free_y")
