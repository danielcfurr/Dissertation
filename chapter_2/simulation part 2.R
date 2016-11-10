library(readstata13)
library(reshape2)
library(ggplot2)
library(xtable)
library(RColorBrewer)


# Import data from Stata simulation --------------------------------------------

# Function to assemble a data frame from individual .dta files
append_stata_files <- function(pattern, id.vars = c("seed", "condition",
                                    "model", "tau", "nitems", "npersons")) {
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

# Assemble df of initial Rasch fits
rasch_df <- append_stata_files("^results_rasch_[0-9]*.dta$",
                  id.vars = c("seed", "condition", "tau", "nitems", "npersons"))
rasch_df[, c("condition", "npersons")] <- list(NULL)
sum(rasch_df[rasch_df$variable == "converge_fails", "value"])

# Assemble df of results
lltm_df <- append_stata_files("^results_lltm_[0-9]*.dta$")
lltm_df <- cbind(method = "lltm", lltm_df)
twostage_df <- append_stata_files("^results_twostage_[0-9]*.dta$")
twostage_df <- cbind(method = "twostage", twostage_df)
# df <- subset(rbind(lltm_df, twostage_df), ! seed %in% errors$seed)
df <- rbind(lltm_df, twostage_df)
df[, c("condition", "npersons")] <- list(NULL)

# Get a vector of seeds for non-converging iterations
v_1 <- df[df$variable == "converge_fails", "value"]
s_1 <- df$seed[v_1 > 0 | is.na(v_1)]
v_2 <- rasch_df[rasch_df$variable == "converge_fails", "value"]
s_2 <- rasch_df$seed[v_2 > 0 | is.na(v_2)]
bad_seeds <- unique(c(s_1, s_2))

# Remove non-converging iterations and reduce each condition to 500 replications.
# (Stata simulation runs a few extra replications to get 500 good replications.)
seeds_df <- subset(rasch_df, ! seed %in% bad_seeds)
kept_seeds <- tapply(seeds_df$seed,
                     list(seeds_df$tau, seeds_df$nitems),
                     function(x) x[1:500])
kept_seeds <- unlist(kept_seeds)
rasch_df <- subset(rasch_df, seed %in% kept_seeds)
df <- subset(df, seed %in% kept_seeds)

variables_to_save <- c("bad_seeds")


# LR tests: descriptive vs explanatory models ----------------------------------

df_lr1 <- subset(df, variable == "dev_train")
df_lr1$variable <- NULL
names(df_lr1)[names(df_lr1) == "value"] <- "dev_exp"

df_lr2 <- subset(rasch_df, variable %in% c("dev_meta", "dev_rasch"))
df_lr2$method <- ifelse(df_lr2$variable == "dev_meta", "twostage", "lltm")
df_lr2$variable <- NULL
names(df_lr2)[names(df_lr2) == "value"] <- "dev_des"

df_lr <- merge(df_lr1, df_lr2)
df_lr$dif <- df_lr$dev_exp - df_lr$dev_des
df_lr$df <- df_lr$nitems - (4 + df_lr$model)
df_lr$pval <- 1 - pchisq(df_lr$dif, df = df_lr$df)
df_lr$reject <- df_lr$pval < .05

aggregate(reject ~ tau + nitems + model + method, df_lr, function(x) mean(x)*100)


# Plot labels ------------------------------------------------------------------

refactor <- function(x, label_list) {
  factor(label_list[as.character(x)], levels = label_list)
}

method_labels <- c("lltm" = "LLTM", "twostage" = "LLTM-E2S")

nitems_labels <- paste("I =", sort(unique(df$nitems)))
names(nitems_labels) <- sort(unique(df$nitems))

tau_labels <- paste("tau =", sort(unique(df$tau)))
names(tau_labels) <- sort(unique(df$tau))


# Plots for parameter recovery -------------------------------------------------

# Set up data frame for parameter estimates
df_recov <-subset(df, grepl("_(est|se)$", variable) & !is.na(value))
df_recov$par <- sub("_(est|se)$", "", df_recov$variable)
df_recov$variable <- sub("^.*_", "", df_recov$variable)
df_recov <- dcast(df_recov, method + seed + model + tau + nitems + par ~ variable)

# Change variance estimates to SDs, remove SE for these
is_variance <- grepl("sq$", df_recov$par)
df_recov$est[is_variance] <- sqrt(df_recov$est[is_variance])
df_recov$se[is_variance] <- NA
df_recov$par <- sub("sq$", "", df_recov$par)

# Change beta estimates to difficulty rather than easiness scale
df_recov$est[!is_variance] <- -1 * df_recov$est[!is_variance]

# Append generating values and get difference with estimated values and z
genvalues <- c(beta1 = -.5, beta2 = 1, beta3 = .5, beta4 = .5, beta5 = -.5,
               beta6 = 0, sigma = 1, tau = NA)
df_recov$gen <- genvalues[df_recov$par]
df_recov$gen[df_recov$par == "tau"] <- df_recov$tau[df_recov$par == "tau"]
df_recov$dif <- df_recov$est - df_recov$gen
df_recov$z <- df_recov$est / df_recov$se

# Plot formatting
df_recov$Method <- refactor(df_recov$method, method_labels)
par_colors <- brewer.pal(7, "Set1")
par_labels <- list()
for(i in 1:5) par_labels[[i]] <- bquote(beta[.(i)])
par_labels[[6]] <- bquote(sigma)
par_labels[[7]] <- bquote(tau)

# Aggregate data frame for bias plots
mean_and_ci <- function(x) {
  n <- length(x)
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  return(c(mean = m, lower = m - 1.96*se, upper = m + 1.96*se))
}
df_bias <- aggregate(dif ~ Method + model + tau + nitems + par,
                     data = subset(df_recov, model == 2),
                     FUN = mean_and_ci)
df_bias <- cbind(df_bias[, -ncol(df_bias)], df_bias$dif)

variables_to_save <- c(variables_to_save, "df_recov", "df_bias", "par_colors",
                       "par_labels")


# Plots for model selection ----------------------------------------------------

# Make data frame for selection by HV, CV, and IC
df_hvcv <- subset(df, grepl("^(aic|bic|dev_loo|dev_(new|same)items)", variable))
df_min <- aggregate(value ~ method + seed + variable, data = df_hvcv, FUN = min)
names(df_min)[names(df_min) == "value"] <- "minimum"
df_hvcv <- merge(df_hvcv, df_min)
df_hvcv$selected <- df_hvcv$value == df_hvcv$minimum
df_hvcv[, c("value", "minimum")] <-list(NULL)
names(df_hvcv)[names(df_hvcv) == "variable"] <- "selector"
df_hvcv$selector <- gsub("^dev_validation_", "", df_hvcv$selector)

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
df_chi2 <- subset(df, variable == "dev_train")
# Ensure is ordered by model
df_chi2 <- df_chi2[with(df_chi2, order(seed, method, model)), ]
df_chi2 <- aggregate(value ~ method + seed + tau + nitems,
                     data = df_chi2, FUN = chi2_select)

# Reformat chi2 results to match df_hvcv and append
df_chi2 <- cbind(df_chi2[, -ncol(df_chi2)], df_chi2$value)
df_chi2 <- melt(df_chi2, id.vars = c("method", "seed", "tau", "nitems"),
                variable.name = "model", value.name = "selected")
df_chi2$model <- as.numeric(df_chi2$model)
df_chi2$selector <- "chi2"
df_chi2 <- df_chi2[, names(df_hvcv)]

# Combine hv/cv and LR test data frames, make dataframe for plot
df_sel <- rbind(df_hvcv, df_chi2)
percents_df <- aggregate(selected ~ method + selector + model + tau + nitems,
                         df_sel, function(x) mean(x)*100)

# Format data frame for plotting
selector_labels <- c("dev_sameitems" = "HV same items",
                     "dev_newitems" = "HV new items",
                     "aic" = "AIC",
                     "dev_loo" = "LOCO-CV",
                     "chi2" = "LR test",
                     "bic" = "BIC")
percents_df$Selector <- refactor(percents_df$selector, selector_labels)
percents_df$Method <- refactor(percents_df$method, method_labels)

variables_to_save <- c(variables_to_save, "df_sel", "percents_df")


# Plots for IC/CV penalties ----------------------------------------------------

# Make data frame for loco-cv penalties
df_lococv <- subset(df, grepl("^dev_(train|loo)$", variable) & method == "twostage")
df_lococv$variable <- ifelse(df_lococv$variable == "dev_train", "in_dev", "out_dev")
df_lococv$type <- "loco"

# Make data frame for holdout validation penalties
df_hv <- subset(df, grepl("^dev_(intest|test)$", variable))
df_hv$variable <- ifelse(df_hv$variable == "dev_intest", "in_dev", "out_dev")
df_hv$type <- "hv"

# Combine the two, aggregate to mean penalty +- 1 SD
df_pen <- rbind(df_lococv, df_hv)
df_pen <- dcast(df_pen, ... ~ variable, value.var = "value")
df_pen$p <- df_pen$out_dev - df_pen$in_dev
mean_and_ci <- function(x) {
  m <- mean(x)
  s <- sd(x)
  return(c(mean = m, lower = m-s, upper = m+s))
}
df_pen <- aggregate(p ~ method + model + tau + nitems + type,
                    data = df_pen, mean_and_ci)
df_pen <- cbind(df_pen[, -ncol(df_pen)], df_pen[, ncol(df_pen)])

# Assemble data frame for plotting. Top is means without SD, bottom is means
# and SD for only Model 2.
side_labels <- c("Est" = "Estimated penalties",
                 "SD" = "\u00B1 1 SD for Model 2 penalty")
df_pen1 <- df_pen
df_pen1$lower <- df_pen1$upper <- NA
df_pen1$side <- "Est"
df_pen2 <- subset(df_pen, model == 2)
df_pen2$side <- "SD"
df_penstack <- rbind(df_pen1, df_pen2)

# Add labels to data frame to plot
# df_penstack_labels <- c("lltm hv" = "LLTM: HV (new items)",
#                         "twostage hv" = "LLTM-E2S: HV (new items)",
#                         "twostage loco" = "LLTM-E2S: LOCO-CV")
df_penstack_labels <- c("lltm hv" = "LLTM",
                        "twostage hv" = "LLTM-E2S")
df_penstack$Type <- refactor(paste(df_penstack$method, df_penstack$type),
                             df_penstack_labels)
df_penstack$text <- ifelse(df_penstack$nitems == 64 | df_penstack$tau == 1,
                           paste0("M", df_penstack$model), "")
df_penstack$side <- refactor(df_penstack$side, side_labels)

# Add AIC penalty to data frame to plot
aic_penalty <- 2*(5:7)
df_penstack$aic_penalty <- aic_penalty[df_penstack$model]

# Remove LOCO-CV from the plots
# df_penstack <- subset(df_penstack, type == "hv")

df_penstack <- df_penstack[,1:6]
df_aic <- subset(df_penstack, type == "hv")
df_aic$type = "aic"
df_aic$mean <- (df_aic$model + 4) * 2
df_penstack <- rbind(df_penstack, df_aic)
df_penstack$group <- paste(df_penstack$method, df_penstack$type)
df_penstack$Model <- factor(df_penstack$model, 1:3)
df_penstack$type[df_penstack$type == "hv"] <- "HV new items"
df_penstack$type[df_penstack$type == "aic"] <- "AIC"
df_penstack$type[df_penstack$type == "loco"] <- "LOCO-CV"
df_penstack$Type <- factor(df_penstack$type, c("HV new items", "AIC", "LOCO-CV"))
df_penstack$method <- ifelse(df_penstack$method == "lltm", "LLTM", "LLTM-E2S")

variables_to_save <- c(variables_to_save, "df_penstack")


# Data frame for prediction ----------------------------------------------------

df_frmse_m <- subset(df, variable == "rmse_fix")
df_frmse_m$variable <- NULL
df_frmse_m <- aggregate(value ~ method + tau + nitems + model, df_frmse_m, mean)
df_frmse_m$Model <- as.factor(df_frmse_m$model)
df_frmse_m$Method <- refactor(df_frmse_m$method, method_labels)

df_frmse <- merge(subset(df_sel, selected), subset(df, variable == "rmse_fix"))
df_frmse <- aggregate(value ~ selector + tau + nitems + method, df_frmse, mean)
df_frmse$Selector <- refactor(df_frmse$selector, selector_labels)
df_frmse$Method <- refactor(df_frmse$method, method_labels)

variables_to_save <- c(variables_to_save, "df_frmse_m", "df_frmse")
save(list = variables_to_save, file = "simulation part 2.Rdata")


#  ------

# LR test "penalty"
x <- seq(from = 2, to = 4, by = .05)
plot(x, pchisq(x, 1, lower.tail = FALSE))
abline(h = .05)

# BIC penalties
log(500*32)
log(32)

