sim_bs_r <- 500
sim_bs_l <- 1:100
sim_bs_I <- 25

brute_force_reps <- 100

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

options(loo.cores = n_cores)
options(mc.cores = n_cores)

compile <- stan_model(file = "rim_mvn.stan")

data_list <- rim_simulate(sim_bs_I, sim_J, sim_sigma, sim_psi, sim_beta)

# Simulation on bootstrap ------------------------------------------------------

waic_marg_list <- waic_cond_list <- list()

sim_bs_start <- Sys.time()
for(i in 1:length(sim_bs_l)) {

  message(Sys.time(), " block sizes; i = ", i, ", l = ", sim_bs_l[i])

  fit <- sampling(compile, data_list, chains = n_chains, iter = n_iter,
                  warmup = n_warmup)
  post_draws <- extract(fit)

  # Marginal WAIC
  cl <- makeCluster(n_cores, type='PSOCK')
  bs <- tsboot(post_draws$mll_j, statistic = waic_wrapper, l = sim_bs_l[i],
               R = sim_bs_r, sim = "fixed", endcorr = TRUE,
               cl = cl, parallel = "snow", ncpus = n_cores)
  stopCluster(cl)
  mcerror <- apply(bs$t[,1:3], 2, sd)
  names(mcerror) <- paste0("mce.", names(bs$t0[1:3]))
  waic_marg_list[[i]] = c(l = sim_bs_l[i], bs$t0[1:3], mcerror)

  # Conditional WAIC
  cl <- makeCluster(n_cores, type='PSOCK')
  bs <- tsboot(post_draws$cll_ij, statistic = waic_wrapper, l = sim_bs_l[i],
               R = sim_bs_r, sim = "fixed", endcorr = TRUE,
               cl = cl, parallel = "snow", ncpus = n_cores)
  stopCluster(cl)
  mcerror <- apply(bs$t[,1:3], 2, sd)
  names(mcerror) <- paste0("mce.", names(bs$t0[1:3]))
  waic_cond_list[[i]] = c(l = sim_bs_l[i], bs$t0[1:3], mcerror)

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
                        "df_bootstrap")


# Brute force approximation for mc error ---------------------------------------

brute_force_waic_marg_list <- brute_force_waic_cond_list <- list()

brute_force_start <- Sys.time()
for(i in 1:brute_force_reps) {

  message(Sys.time(), " brute force; i = ", i)

  fit <- sampling(compile, data_list, chains = n_chains, iter = n_iter,
                  warmup = n_warmup)
  post_draws <- extract(fit, pars = c("mll_j", "cll_ij"))

  brute_force_waic_marg_list[[i]] <- waic_wrapper(post_draws$mll_j)[1:3]
  brute_force_waic_cond_list[[i]] <- waic_wrapper(post_draws$cll_ij)[1:3]

}

brute_force_end <- Sys.time()
brute_force_end - brute_force_start

df_waic_marg <- data.frame(do.call(rbind, brute_force_waic_marg_list),
                           focus = "Marginal", IC = "WAIC")
df_waic_cond <- data.frame(do.call(rbind, brute_force_waic_cond_list),
                           focus = "Conditional", IC = "WAIC")
df_brute_force <- rbind(df_waic_marg, df_waic_cond)

variables_to_save <- c(variables_to_save, "brute_force_start",
                        "brute_force_end", "df_brute_force")

save(list = variables_to_save, file = "bootstrap.Rdata")


#

df_brute_force$i <- 1:nrow(df_brute_force)
df_brute_force_long <- melt(df_brute_force, id.vars = c("focus", "IC", "i"))
p_brute_marg <- ggplot(subset(df_brute_force_long, focus == "Marginal")) +
  aes(x = value) +
  geom_histogram(bins = 10) +
  facet_wrap(~focus + variable, scales = "free_x")
p_brute_cond <- p_brute_marg %+% subset(df_brute_force_long,
                                        focus == "Conditional")
p_brute_marg
p_brute_cond

df_bootstrap_long <- melt(df_bootstrap, id.vars = c("focus", "IC", "l"))
p_boot_marg <- ggplot(subset(df_bootstrap_long, focus == "Marginal")) +
  aes(x = l, y = value) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~focus + variable, scales = "free_y")
p_boot_cond <- p_boot_marg %+% subset(df_bootstrap_long, focus == "Conditional")
p_boot_marg
p_boot_cond

df_brute_force_mce <- dcast(df_brute_force_long, focus + IC + variable ~ .,
                            fun.aggregate = sd, value.var = "value")
names(df_brute_force_mce)[names(df_brute_force_mce) == "."] <- "brute_force"
df_brute_force_mce$variable <- paste0("mce.", df_brute_force_mce$variable)
df_combine <- merge(df_bootstrap_long, df_brute_force_mce)
p_comb_marg <- ggplot(subset(df_combine, focus == "Marginal")) +
  aes(x = l, y = value) +
  geom_hline(aes(yintercept = brute_force), linetype = "dashed") +
  geom_point() +
  geom_smooth() +
  facet_wrap(~focus + variable, scales = "free_y")
p_boot_cond <- p_comb_marg %+% subset(df_combine, focus == "Conditional")
p_comb_marg
p_boot_cond


