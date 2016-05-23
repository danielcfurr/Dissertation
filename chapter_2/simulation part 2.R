library(readstata13)
library(reshape2)
library(ggplot2)

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

# Assemble df of failed iterations
errors <- append_stata_files("^results_errors_[0-9]*.dta$", id.vars = c())

# Assemble df of results
lltm_df <- append_stata_files("^results_lltm_[0-9]*.dta$")
lltm_df <- cbind(method = "lltm", lltm_df)
twostage_df <- append_stata_files("^results_twostage_[0-9]*.dta$")
twostage_df <- cbind(method = "twostage", twostage_df)
df <- subset(rbind(lltm_df, twostage_df), ! seed %in% errors$seed)

# Sanity check. Should have same number of results per seed.
table(tapply(df$seed, df$seed, length))

# Check on number of failed iterations
length(unique(errors$seed))

# Add column to df for which model was selected
minimums <- aggregate(df$value,
                      by = list(method = df$method,
                                seed = df$seed,
                                variable = df$variable),
                      FUN = min)
names(minimums)[names(minimums) == "x"] <- "min"
df <- merge(df, minimums)
df$relative <- df$value - df$min
df$selected <- df$relative == 0

# Sanity check. Should equal 1/3
mean(df$relative == 0, na.rm=T)

# Create and append results for likelihood ratio test
chi2_select <- function(devs, df = 1) {
  n <- length(devs)
  differences = devs[1:(n-1)] - devs[2:n]
  p = 1 - pchisq(differences, df = df)
  non_sig <- which(p > .05) # May have length of zero
  chosen_model <- ifelse(length(non_sig) == 0, n, min(non_sig))
  #chosen_vector <- 1:n == chosen_model
  return(chosen_model)
}
chi2_df <- subset(df, variable == "dev_insample")
agg_df <- aggregate(chi2_df$value,
                    by = list(method = chi2_df$method,
                              seed = chi2_df$seed,
                              variable = chi2_df$variable),
                    FUN = chi2_select)
chi2_df <- merge(chi2_df, agg_df)
chi2_df$variable <- "chi2"
chi2_df$selected <- chi2_df$model == chi2_df$x
chi2_df$value <- chi2_df$min <- chi2_df$relative <- NA
chi2_df$x <- NA
chi2_df <- chi2_df[, names(df)]
df <- rbind(df, chi2_df)


# Attatch labels to main data.frame --------------------------------------------

refactor <- function(x, label_list) {
  factor(label_list[as.character(x)], levels = label_list)
}

variable_labels <- c("dev_insample" = "In-sample deviance",
                     "chi2" = "LR test",
                     "aic" = "AIC",
                     "bic" = "BIC",
                     "dev_newpersons" = "HV over persons",
                     "dev_newitems" = "HV over items",
                     "dev_loo" = "LOO CV",
                     "rmse_insample" = "rmse_insample",
                     "rmse_newitems" = "rmse_newitems",
                     "rmse_newpersons" = "rmse_newpersons",
                     "rmse_loo" = "rmse_loo",
                     "rmsea" = "rmsea")
df$Variable <- refactor(df$variable, variable_labels)

method_labels <- c("lltm" = "LLTM", "twostage" = "Two stage")
df$Method <- refactor(df$method, method_labels)

nitems_labels <- paste("I =", sort(unique(df$nitems)))
names(nitems_labels) <- sort(unique(df$nitems))
df$Nitems <- refactor(df$nitems, nitems_labels)

tau_labels <- paste("tau =", sort(unique(df$tau)))
names(tau_labels) <- sort(unique(df$tau))
df$Tau <- refactor(df$tau, tau_labels)


# Plots of selection -----------------------------------------------------------

percents_df <- subset(df, grepl("^(dev|aic|bic|chi2)", variable) &
                      variable != "dev_insample")
percents_df$seed <- percents_df$value <- percents_df$min <-
  percents_df$relative <- NULL
percents_df <- aggregate(selected ~ ., percents_df, function(x) mean(x)*100)

ggplot(subset(percents_df, percents_df$nitems == 32)) +
  aes(as.factor(tau), selected, fill = factor(model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(list(x = expression(paste("Residual item SD (", tau, ")")),
            y = "Percentage of times selected", fill = "Model")) +
  facet_grid(Variable ~ Method)

ggplot(subset(percents_df, percents_df$tau == .5)) +
  aes(factor(nitems), selected, fill = factor(model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(list(x = "Number of items",
            y = "Percentage of times selected", fill = "Model")) +
  facet_grid(variable ~ method)


# Plots of prediction ----------------------------------------------------------

df_rmse <- subset(df, grepl("^rmse_new(items|persons)", variable) & selected)
df_rmse$min <- df_rmse$relative <- df_rmse$Variable <-  df_rmse$model <- NULL
df_rmse_points <- dcast(df_rmse, ... ~ variable, value.var = "value")
df_rmse$seed <- NULL
df_rmse_means <- dcast(df_rmse, ... ~ variable, fun.aggregate = mean,
                       value.var = "value")

limits <- range(df_rmse_points[df_rmse_points$nitems == 32,
                               c("rmse_newitems", "rmse_newpersons")])
ggplot(subset(df_rmse_points, nitems == 32)) +
  aes(x = rmse_newitems, y = rmse_newpersons) +
  geom_point(alpha = .1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = rmse_newitems, y = rmse_newpersons),
             subset(df_rmse_means, nitems == 32), col = "red", pch = 3) +
  coord_fixed(xlim = limits, ylim = limits) +
  labs(list(x = "RMSE with HV over items",
            y = "RMSE with HV over persons")) +
  facet_grid(Tau ~ Method)

limits <- range(subset(df_rmse_points, tau == .5,
                       select = c("rmse_newitems", "rmse_newpersons")))
ggplot(subset(df_rmse_points, tau == .5)) +
  aes(x = rmse_newitems, y = rmse_newpersons) +
  geom_point(alpha = .1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = rmse_newitems, y = rmse_newpersons),
             subset(df_rmse_means, tau == .5), col = "red", pch = 3) +
  coord_fixed(xlim = limits, ylim = limits) +
  labs(list(x = "RMSE with HV over items",
            y = "RMSE with HV over persons")) +
  facet_grid(Nitems ~ Method)


# Plot tau against R^2 ---------------------------------------------------------

# Get variance of fixed part of items
ex <- read.dta13("example.dta")
xB <- ex[ex$person == 1, "xB"]
# Using population variance:
upsilon_sq <- (length(xB)-1) / length(xB) * var(xB)

# Data frame and function to plot
r_sq <- function(tau, upsilon_sq) upsilon_sq / (tau^2 + upsilon_sq)
tau <- sort(unique(df$tau))
pts <- data.frame(tau = tau, r_sq = r_sq(tau, upsilon_sq))

xlim <- c(0, 1.5)
ggplot() +
  scale_x_continuous(limits = xlim) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(data = data.frame(pts),
             mapping = aes(x = tau, y = r_sq),
             size = 2) +
  stat_function(data = data.frame(x = xlim),
                mapping = aes(x),
                fun = r_sq,
                args = list(upsilon_sq = upsilon_sq)) +
  xlab(expression(paste("Residual item SD (", tau, ")"))) +
  ylab(expression(paste("Explained variance (", R^2, ")")))
# ggsave("figs/rsq_vs_tau.pdf", family = "Times", width = 3, height = 2,
#        units = "in", pointsize = 9)



