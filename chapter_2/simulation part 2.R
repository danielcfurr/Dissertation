library(readstata13)
library(reshape2)


# Import data from Stata simulation --------------------------------------------

# Function to assemble a data frame from individual .dta files
append_stata_files <- function(pattern, 
                               id.vars = c("seed", "condition", "model")) {
  dta_files <- dir(pattern = pattern)
  data_list <- lapply(dta_files, read.dta13)
  data_stacked <- do.call(rbind, data_list)
  if(length(id.vars) > 0) {
    data_melted <- melt(data_stacked, id.vars = id.vars)
    return(data_melted)
  } else {
    return(data_stacked)
  }
}

# Assemble df describing conditions
conditions_df <- read.dta13("conditions.dta")

# Assemble df of initial Rasch fits
sim_df <- append_stata_files("^results_sim_[0-9]*.dta$", c("seed", "condition"))

# Assemble df of results
lltm_df <- append_stata_files("^results_lltm_[0-9]*.dta$")
lltm_df <- cbind(method = "lltm", lltm_df)
twostage_df <- append_stata_files("^results_twostage_[0-9]*.dta$")
twostage_df <- cbind(method = "twostage", twostage_df)
df <- merge(conditions_df, rbind(lltm_df, twostage_df))
df$over_Rsq <- df$condition %in% c(1, 2, 3)
df$over_nitems <- df$condition %in% c(2, 4, 5)
df[, c("Vsq", "upsilon", "tau", "b", "sigma", "npersons")] <- list(NULL)

# Get a vector of seeds for non-converging iterations
v_1 <- subset(df, variable == "converge_fails", "value")
s_1 <- df$seed[v_1 > 0 | is.na(v_1)]
v_2 <- subset(sim_df, variable == "converge_fails", "value")
s_2 <- sim_df$seed[v_2 > 0 | is.na(v_2)]
bad_seeds <- unique(c(s_1, s_2))
length(bad_seeds)

# Remove non-converging iterations and reduce each condition to 500 replications.
# (Stata simulation runs a few extra replications to get 500 good replications.)
seeds_df <- subset(sim_df, ! seed %in% bad_seeds)
kept_seeds <- tapply(seeds_df$seed, list(seeds_df$condition), 
                     function(x) x[1:500])
kept_seeds <- unlist(kept_seeds)
sim_df <- subset(sim_df, seed %in% kept_seeds)
df <- subset(df, seed %in% kept_seeds)

# Get number of betas by model
max_beta <- function(x) max(as.numeric(gsub("(beta|_est)", "", x)))
n_betas <- with(subset(df, !is.na(value) & grepl("^beta[0-9]_est$", variable)), 
                tapply(variable, model, max_beta))
n_pars <- n_betas + 1

variables_to_save <- c("conditions_df", "bad_seeds", "n_pars", "n_betas")


# Adjust AIC for twostage method -----------------------------------------------

aic_c <- function(dev, n_obs, n_pars) {
  dev + 2*(n_pars+ 1) * n_obs / (n_obs - n_pars - 2)
}

aic_c_pick <- df$method == "twostage" & df$variable == "aic"
# Remove standard penalty and add corrected one
df_aic <- df[aic_c_pick,]
aic_c_values <- aic_c(dev = df_aic$value - 2*n_pars[df_aic$model], 
                      n_obs = df_aic$nitems, n_pars = n_pars[df_aic$model])
df[aic_c_pick, "value"] <- aic_c_values

variables_to_save <- c(variables_to_save, "aic_c")


# LR tests: descriptive vs explanatory models ----------------------------------

# Make dataframe of LLTM and LLTM-E deviances in training data
df_lr1 <- subset(df, variable == "dev_training")
df_lr1$variable <- NULL
names(df_lr1)[names(df_lr1) == "value"] <- "dev_explanatory"

# Make dataframe of Rasch deviances in training data
df_lr2 <- subset(sim_df, variable %in% c("dev_meta", "dev_rasch"))
df_lr2$method <- ifelse(df_lr2$variable == "dev_meta", "twostage", "lltm")
df_lr2$variable <- NULL
names(df_lr2)[names(df_lr2) == "value"] <- "dev_rasch"

# Combine and do likelihood ratio tests
df_lr <- merge(df_lr1, df_lr2)
df_lr$dif <- df_lr$dev_explanatory - df_lr$dev_rasch
df_lr$df <- df_lr$nitems - (n_pars[df_lr$model])
df_lr$pval <- 1 - pchisq(df_lr$dif, df = df_lr$df)
df_lr$reject <- df_lr$pval < .05

# What proportion of times explanatory model is rejected in favor of Rasch?
aggregate(reject ~ Rsq + nitems + model + method, df_lr, 
          function(x) mean(x)*100)


# Plot labels ------------------------------------------------------------------

refactor <- function(x, label_list) {
  factor(label_list[as.character(x)], levels = label_list)
}

method_labels <- c("lltm" = "LLTM", "twostage" = "LLTM-E2S")

# nitems_labels <- paste("I =", sort(unique(df$nitems)))
# names(nitems_labels) <- sort(unique(df$nitems))

# tau_labels <- paste("tau =", sort(unique(df$tau)))
# names(tau_labels) <- sort(unique(df$tau))


# Plots for parameter recovery -------------------------------------------------

# Set up data frame for parameter estimates and SE for model 2
df_recov <- subset(df, grepl("_(est|se)$", variable) & !is.na(value))
df_recov$par <- sub("_(est|se)$", "", df_recov$variable)
df_recov$variable <- sub("^.*_", "", df_recov$variable)
df_recov <- dcast(df_recov, ... ~ variable, value.var = "value")

# Change variance estimates to SDs, remove SE for these
is_variance <- grepl("sq$", df_recov$par)
df_recov$est[is_variance] <- sqrt(df_recov$est[is_variance])
df_recov$se[is_variance] <- NA
df_recov$par <- sub("sq$", "", df_recov$par)

# Change beta estimates to difficulty rather than easiness scale
df_recov$est[!is_variance] <- -1 * df_recov$est[!is_variance]

# Append generating values and get difference with estimated values and z
df_gen <- conditions_df
df_gen$beta1 <- 0
for(i in 2:n_betas[2]) df_gen[, paste0("beta", i)] <- df_gen$b
for(i in (n_betas[2]+1):n_betas[3]) df_gen[, paste0("beta", i)] <- 0
df_gen <- melt(df_gen, id.vars = "condition",
               measure.vars = c(paste0("beta", 1:n_betas[3]), "sigma", "tau"))
names(df_gen)[names(df_gen) == "value"] <- "generating"
names(df_gen)[names(df_gen) == "variable"] <- "par"
df_recov <- merge(df_recov, df_gen)           
df_recov$dif <- df_recov$est - df_recov$generating
df_recov$z <- df_recov$est / df_recov$se
df_recov$Method <- refactor(df_recov$method, method_labels)


# Data frame for assessing bias in model 2
mean_and_ci <- function(x) {
  n <- length(x)
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  return(c(mean = m, lower = m - 1.96*se, upper = m + 1.96*se))
}
df_bias <- aggregate(dif ~ .,
                     data = subset(df_recov, model == 2,
                                   select = -c(seed, est, se, generating, z)),
                     FUN = mean_and_ci)
df_bias <- cbind(df_bias[, -ncol(df_bias)], df_bias$dif)

# Data frame for assessing standard errors for last beta (beta = 0)
last_beta <- paste0("beta", n_betas[3])
df_qq <- subset(df_recov, par == last_beta)


variables_to_save <- c(variables_to_save, "df_bias", "df_qq")


# Plots for model selection ----------------------------------------------------

# Make data frame for selection by HV, CV, and IC
df_hvcv <- subset(df, grepl("^(aic|bic|dev_loo|dev_(new|same)items)", variable))
df_min <- aggregate(value ~ method + seed + variable, data = df_hvcv, FUN = min)
names(df_min)[names(df_min) == "value"] <- "minimum"
df_hvcv <- merge(df_hvcv, df_min)
df_hvcv$selected <- df_hvcv$value == df_hvcv$minimum
df_hvcv[, c("value", "minimum")] <-list(NULL)
names(df_hvcv)[names(df_hvcv) == "variable"] <- "selector"
df_hvcv$selector <- gsub("^dev_", "", df_hvcv$selector)

# Create and append results for likelihood ratio test
chi2_select <- function(devs, df = 1) {
  n <- length(devs)
  differences <- devs[1:(n-1)] - devs[2:n]
  p <- 1 - pchisq(differences, df = df)
  non_sig <- which(p > .05) # May have length of zero
  chosen_model <- ifelse(length(non_sig) == 0, n, min(non_sig))
  chosen_vector <- 1:n == chosen_model
  return(chosen_vector)
}
df_chi2 <- subset(df, variable == "dev_training")
# Ensure is ordered by model
df_chi2 <- df_chi2[with(df_chi2, order(seed, method, model)), ]
df_chi2$model <- NULL
df_chi2 <- aggregate(value ~ ., data = df_chi2, FUN = chi2_select)
df_chi2 <- cbind(df_chi2[, -ncol(df_chi2)], df_chi2$value)

# Reformat chi2 results to match df_hvcv and append
df_chi2 <- melt(df_chi2, measure.vars = c("1", "2", "3"),
                variable.name = "model", value.name = "selected")
df_chi2$model <- as.numeric(df_chi2$model)
df_chi2$selector <- "chi2"
df_chi2 <- df_chi2[, names(df_hvcv)]

# Combine hv/cv and LR test data frames, make dataframe for plot
df_sel <- rbind(df_hvcv, df_chi2)
percents_df <- aggregate(selected ~ ., subset(df_sel, select = -seed),
                         function(x) mean(x)*100)

# Format data frame for plotting
selector_labels <- c("sameitems" = "HV same items",
                     "newitems" = "HV new items",
                     "aic" = "AIC",
                     "loo" = "LOCO-CV",
                     "chi2" = "LR test",
                     "bic" = "BIC")
percents_df$Selector <- refactor(percents_df$selector, selector_labels)
percents_df$Method <- refactor(percents_df$method, method_labels)

variables_to_save <- c(variables_to_save, "percents_df")


# Plots for IC/CV penalties ----------------------------------------------------

# Make data frame for loco-cv penalties
df_lococv <- subset(df, grepl("^dev_(training|loo)$", variable) & 
                      method == "twostage")
df_lococv$variable <- ifelse(df_lococv$variable == "dev_training", "in_dev", 
                             "out_dev")
df_lococv$type <- "loco"

# Make data frame for holdout validation penalties
df_hv <- subset(df, grepl("^dev_(ineval|evaluation)$", variable))
df_hv$variable <- ifelse(df_hv$variable == "dev_evaluation", "out_dev", "in_dev")
df_hv$type <- "hv"

# Combine the two, aggregate to mean penalty
df_pen <- rbind(df_lococv, df_hv)
df_pen <- dcast(df_pen, ... ~ variable, value.var = "value")
df_pen$p <- df_pen$out_dev - df_pen$in_dev
df_pen <- aggregate(p ~ ., 
                    data = subset(df_pen, select = -c(seed, in_dev, out_dev)), 
                    mean)

# Attach AIC penalties to penalty df
df_aic <- subset(df_pen, type == "hv" & method == "lltm")
df_aic$type = "aic"
df_aic$p <- (n_pars[df_aic$model]) * 2
df_aic_c <- subset(df_pen, type == "hv" & method == "twostage")
df_aic_c$type = "aic"
df_aic_c$p <- with(df_aic_c, 2*(n_pars[model] + 1) * nitems / 
                     (nitems - n_pars[model] - 2))
df_pen <- rbind(df_pen, df_aic, df_aic_c)

# Format df
df_pen$group <- paste(df_pen$method, df_pen$type)
df_pen$Model <- factor(df_pen$model, 1:3)
df_pen$type[df_pen$type == "hv"] <- "HV new items"
df_pen$type[df_pen$type == "aic"] <- "AIC"
df_pen$type[df_pen$type == "loco"] <- "LOCO-CV"
df_pen$Type <- factor(df_pen$type, c("HV new items", "AIC", "LOCO-CV"))
df_pen$method <- ifelse(df_pen$method == "lltm", "LLTM", "LLTM-E2S")

variables_to_save <- c(variables_to_save, "df_pen")


# Data frames for prediction ----------------------------------------------------

# Get mean RMSE for selected model by selection strategy
df_rmse <- merge(subset(df_sel, selected), subset(df, variable == "rmse_full"))
df_rmse <- aggregate(value ~ ., subset(df_rmse, select = -c(seed, model)), 
                     mean)
df_rmse$Selector <- refactor(df_rmse$selector, selector_labels)
df_rmse$Method <- refactor(df_rmse$method, method_labels)
df_rmse <- merge(df_rmse, conditions_df)

# Get mean RMSE for predictions by model by selection strategy
df_rmse_m <- subset(df, variable == "rmse_full")
df_rmse_m$variable <- NULL
df_rmse_m <- aggregate(value ~ ., subset(df_rmse_m, select = -seed), mean)
df_rmse_m$Model <- as.factor(df_rmse_m$model)
df_rmse_m$Method <- refactor(df_rmse_m$method, method_labels)
df_rmse_m <- merge(df_rmse_m, conditions_df)

variables_to_save <- c(variables_to_save, "df_rmse_m", "df_rmse")
save(list = variables_to_save, file = "simulation part 2.Rdata")


#  ------

# # LR test "penalty"
# x <- seq(from = 2, to = 4, by = .05)
# plot(x, pchisq(x, 1, lower.tail = FALSE))
# abline(h = .05)
# 
# # BIC penalties
# log(500*32)
# log(32)

