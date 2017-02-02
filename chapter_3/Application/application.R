n_posterior <- 5000
n_warmup <- 500
n_chains <- 5
n_iter <- n_posterior / n_chains + n_warmup
n_cores <- 5

n_agq_try <- round(7*1.5^(0:2))

n_bs_r <- 500
n_block_size_reps <- 5
tsboot_sim <- "fixed"

start <- Sys.time()

variables_to_save <- ls()

# set.seed(72793333)


# Set up -----------------------------------------------------------------------

library(rstan)
library(edstan)
library(loo)
library(boot)
library(doParallel)
library(foreach)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = n_cores)
options(loo.cores = 1)

source("../functions.R")


# Set up data and models -------------------------------------------------------

# Person regression models to try
model_list <- list(~ 1,
                   ~ 1 + anger,
                   ~ 1 + male,
                   ~ 1 + anger + male,
                   ~ 1 + anger + male + I(anger*male))


# Select number of adaptive quadrature points ----------------------------------

dl <- irt_data(y = aggression$dich, jj = aggression$person,
               ii = aggression$item, covariates = aggression,
               formula = model_list[[4]])
fit <- irt_stan(dl, "rasch_edstan_modified.stan", iter = n_iter,
                chains = n_chains, warmup = n_warmup)

# Get Marginal IC results for each choice of number AGQ nodes
agq_try <- matrix(NA, ncol = 4, nrow = length(n_agq_try))
colnames(agq_try) <- c("nodes", "dic", "waic", "looic")
agq_try[, "nodes"] <- n_agq_try
agq_change <- agq_try[-1, ]
for(i in 1:length(n_agq_try)) {
  message(Sys.time(), ": i=", i)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg_j <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", agq_try[i])
  stopCluster(cl)
  agq_try[i, "dic"] <- dic(ll_marg_j)["dic"]
  agq_try[i, "waic"] <- waic(ll_marg_j$ll)[["waic"]]
  agq_try[i, "looic"] <- loo(ll_marg_j$ll, cores = n_cores)[["looic"]]
  if(i > 1) agq_change[i-1, 2:4] <- agq_try[i, 2:4] - agq_try[i-1, 2:4]
}

# Select first choice of N nodes with < .01 change
acceptable_change <- apply(agq_change[, -1], 1, function(x) all(abs(x) < .01))
first_acceptable <- min(which(acceptable_change == TRUE))
(n_agq <- agq_try[first_acceptable])

variables_to_save <- c(variables_to_save, "agq_try", "n_agq")


# Repeat bootstrap with different block sizes for one model --------------------

# Choose block sizes to test
block_size <- round((n_posterior)^(1/3))
# block_sizes <- c(block_size, block_size - 5, block_size + 5)
# block_sizes <- round(block_size*c(1/2, 1, 2))
block_sizes <- c(1, block_size, 50, 100, 200)


# Get conditional liklihood
ll_cond_ij <- cll(fit, dl, f_conditional_ij)

# Get marginal liklihood from N nodes chosen previously
cl <- makeCluster(n_cores)
registerDoParallel(cl)
ll_marg_j <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", n_agq)
stopCluster(cl)


# Conditional and marginal WAIC

waic_cond_result <- waic_wrapper(ll_cond_ij$ll)

waic_cond_mcerror <- matrix(NA, nrow = n_block_size_reps*length(block_sizes),
                            ncol = length(waic_cond_result) + 1)
colnames(waic_cond_mcerror) <- c("block_size", names(waic_cond_result))
i <- 0
for(r in 1:n_block_size_reps) {
  for(b in 1:length(block_sizes)) {
    i <- i + 1
    message(Sys.time(), ": r=", r, " b=", b, " i=", i)
    bs <- tsboot(ll_cond_ij$ll, statistic = waic_wrapper, l = block_sizes[b],
                 R = n_bs_r, sim = tsboot_sim, parallel = "snow", ncpus = n_cores)
    mcerrors <- apply(bs$t, 2, sd)
    waic_cond_mcerror[i, ] <- c(block_sizes[b], mcerrors)
  }
}

waic_marg_mcerror <- waic_cond_mcerror
waic_marg_mcerror[,] <- NA
i <- 0
for(r in 1:n_block_size_reps) {
  for(b in 1:length(block_sizes)) {
    i <- i + 1
    message(Sys.time(), ": r=", r, " b=", b, " i=", i)
    bs <- tsboot(ll_marg_j$ll, statistic = waic_wrapper, l = block_sizes[b],
                 R = n_bs_r, sim = tsboot_sim, parallel = "snow", ncpus = n_cores)
    mcerrors <- apply(bs$t, 2, sd)
    waic_marg_mcerror[i, ] <- c(block_sizes[b], mcerrors)
  }
}


# Conditional DIC

# Format draws for compatibility with tsboot()
colnames(ll_cond_ij$ll) <- paste0("loglik[", 1:ncol(ll_cond_ij$ll), "]")
mat <- cbind(as.matrix(fit), ll_cond_ij$ll)
dic_cond_ij_wrapper(mat, dl)

# Compare to dic()
(dic_cond_result <- dic(ll_cond_ij))

# Try each block size for one model several times
dic_cond_mcerror <- matrix(NA, nrow = n_block_size_reps*length(block_sizes),
                            ncol = length(dic_cond_result) + 1)
colnames(dic_cond_mcerror) <- c("block_size", names(dic_cond_result))
i <- 0
for(r in 1:n_block_size_reps) {
  for(b in 1:length(block_sizes)) {
    i <- i + 1
    message(Sys.time(), ": r=", r, " b=", b, " i=", i)
    cl <- makeCluster(n_cores, type='PSOCK')
    clusterExport(cl, c("dic", "cll", "f_conditional_ij",
                        "better_posterior_means"))
    bs <- tsboot(mat, statistic = dic_cond_ij_wrapper, l = block_sizes[b],
                 R = n_bs_r, sim = tsboot_sim, data_list = dl,
                 cl = cl, parallel = "snow", ncpus = n_cores)
    stopCluster(cl)
    mcerrors <- apply(bs$t, 2, sd)
    dic_cond_mcerror[i, ] <- c(block_sizes[b], mcerrors)
  }
}


# Marginal DIC

dic_marg_result <- dic(ll_marg_j)

# Format draws for compatibility with tsboot()
colnames(ll_marg_j$ll) <- paste0("loglik[", 1:ncol(ll_marg_j$ll), "]")
mat <- cbind(as.matrix(fit), ll_marg_j$ll)
dic_marg_wrapper(mat, dl)

# Compare to dic()
(dic_marg_result <- dic(ll_marg_j))

dic_marg_mcerror <- dic_cond_mcerror
dic_marg_mcerror[,] <- NA
i <- 0
for(r in 1:n_block_size_reps) {
  for(b in 1:length(block_sizes)) {
    i <- i + 1
    message(Sys.time(), ": r=", r, " b=", b, " i=", i)
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("dic", "mll_serial", "f_marginal_j",
                        "better_posterior_means", "n_agq", "dl"))
    bs <- tsboot(mat, statistic = dic_marg_wrapper, l = block_sizes[b],
                 R = n_bs_r, sim = tsboot_sim, data_list = dl,
                 cl = cl, parallel = "snow", ncpus = n_cores)
    stopCluster(cl)
    mcerrors <- apply(bs$t, 2, sd)
    dic_marg_mcerror[i, ] <- c(block_sizes[b], mcerrors)
  }
}

df_mcerror <- as.data.frame(
  cbind(rbind(dic_cond_mcerror, dic_marg_mcerror),
        rbind(waic_cond_mcerror[,-1], waic_marg_mcerror[,-1]))
)
df_mcerror$focus <- rep(c("Conditional", "Marginal"),
                        each = nrow(dic_cond_mcerror))

variables_to_save <- c(variables_to_save, "df_mcerror", "block_size",
                       "block_sizes")


# Fit the models and get conditional/marginal IC with bootstrapping ------------

cond_dic <- cond_waic <- cond_loo <- list()
marg_dic <- marg_waic <- marg_loo <- list()
summary_list <- list()

i <- 0
for(m in 1:length(model_list)) {

  message(Sys.time(), ", Model ", m, ": Starting fit")

  dl <- irt_data(y = aggression$dich, jj = aggression$person,
                 ii = aggression$item, covariates = aggression,
                 formula = model_list[[m]])
  fit <- irt_stan(dl, "rasch_edstan_modified.stan", iter = n_iter,
                  chains = n_chains, warmup = n_warmup)

  summary_list[[m]] <- summary(fit, probs = c())[[1]]

  message(Sys.time(), ", Model ", m, ": Starting conditional IC")

  # Conditional IC (not aggregated)

  ll_cond_ij <- cll(fit, dl, f_conditional_ij)

  # Conditional DIC
  dic_hold <- unlist(dic(ll_cond_ij))
  # Format draws for compatibility with tsboot()
  colnames(ll_cond_ij$ll) <- paste0("loglik[", 1:ncol(ll_cond_ij$ll), "]")
  mat <- cbind(as.matrix(fit), ll_cond_ij$ll)
  dic_cond_ij_wrapper(mat, dl)
  cl <- makeCluster(n_cores, type="PSOCK")
  clusterExport(cl, c("dic", "cll", "f_conditional_ij",
                      "better_posterior_means"))
  bs <- tsboot(mat, statistic = dic_cond_ij_wrapper, l = block_size,
               R = n_bs_r, sim = tsboot_sim, data_list = dl,
               cl = cl, parallel = "snow", ncpus = n_cores)
  stopCluster(cl)
  mcerrors <- apply(bs$t, 2, sd)
  names(mcerrors) <- paste0(names(dic_hold), "_mcerror")
  cond_dic[[m]] <- c(dic_hold, mcerrors)

  # Conditional WAIC
  waic_hold <- waic_wrapper(ll_cond_ij$ll)
  bs <- tsboot(ll_cond_ij$ll, statistic = waic_wrapper, l = block_size,
               R = n_bs_r, sim = tsboot_sim, parallel = "snow", ncpus = n_cores)
  mcerrors <- apply(bs$t, 2, sd)
  names(mcerrors) <- paste0(names(waic_hold), "_mcerror")
  cond_waic[[m]] <- c(waic_hold, mcerrors)

  # Conditional LOO-IC (No bootstrap performed)
  cond_loo[[m]] <- loo_wrapper(ll_cond_ij$ll, cores = n_cores)

  # Marginal IC aggregated to each j

  message(Sys.time(), ", Model ", m, ": Starting marginal IC")

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  ll_marg <- mll_parallel(fit, dl, f_marginal_j, "zeta", "sigma", n_agq)
  stopCluster(cl)

  # Marginal DIC
  dic_hold <- unlist(dic(ll_marg))
  # Format draws for compatibility with tsboot()
  colnames(ll_marg_j$ll) <- paste0("loglik[", 1:ncol(ll_marg_j$ll), "]")
  mat <- cbind(as.matrix(fit), ll_marg_j$ll)
  dic_marg_wrapper(mat, dl)
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("dic", "mll_serial", "f_marginal_j",
                      "better_posterior_means", "n_agq", "dl"))
  bs <- tsboot(mat, statistic = dic_marg_wrapper, l = block_size,
               R = n_bs_r, sim = tsboot_sim, data_list = dl,
               cl = cl, parallel = "snow", ncpus = n_cores)
  stopCluster(cl)
  mcerrors <- apply(bs$t, 2, sd)
  names(mcerrors) <- paste0(names(dic_hold), "_mcerror")
  marg_dic[[m]] <- c(dic_hold, mcerrors)

  # Marginal WAIC
  waic_hold <- waic_wrapper(ll_marg$ll)
  bs <- tsboot(ll_marg$ll, statistic = waic_wrapper, l = block_size,
               R = n_bs_r, sim = tsboot_sim, parallel = "snow", ncpus = n_cores)
  mcerrors <- apply(bs$t, 2, sd)
  names(mcerrors) <- paste0(names(waic_hold), "_mcerror")
  marg_waic[[m]] <- c(waic_hold, mcerrors)

  # Marginal LOO-IC (No bootstrap performed)
  marg_loo[[m]] <- loo_wrapper(ll_marg$ll, cores = n_cores)

}

# Format IC results

combine_focus <- function(clist, mlist) {
  df_c <- data.frame(focus = "Conditional", model = 1:length(clist),
                     do.call(rbind, clist))
  df_m <- data.frame(focus = "Marginal", model = 1:length(mlist),
                     do.call(rbind, mlist))
  df_long <- melt(rbind(df_c, df_m), id.vars = c("focus", "model"))
  return(df_long)
}

df_models <- rbind(combine_focus(cond_dic, marg_dic),
                   combine_focus(cond_waic, marg_waic),
                   combine_focus(cond_loo, marg_loo))

# Format summary() results

df_list <- list()
for(i in 1:length(summary_list)) {
  df_list[[i]] <- as.data.frame(summary_list[[i]])
  df_list[[i]]$par <- gsub("\\[.*]", "", rownames(df_list[[i]]))
  rownames(df_list[[i]]) <- NULL
  df_list[[i]]$model = i
}
df_summary <- do.call(rbind, df_list)

max_rhat <- max(df_summary$Rhat)

end <- Sys.time()
end - start

variables_to_save <- c(variables_to_save, "df_models", "df_summary", "max_rhat",
                       "start", "end")


save(list = variables_to_save, file = "application.Rdata")


ggplot(df_summary) +
  aes(x = par, color = par, y = Rhat) +
  geom_jitter(height = 0, width = .5, show.legend = FALSE) +
  facet_wrap(~model)

ggplot(df_summary) +
  aes(x = par, color = par, y = n_eff) +
  geom_jitter(height = 0, width = .5, show.legend = FALSE) +
  facet_wrap(~model)


