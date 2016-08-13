library(readstata13)
library(reshape2)
library(ggplot2)
library(xtable)
library(RColorBrewer)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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

# Assemble df of failed iterations
errors <- append_stata_files("^results_errors_[0-9]*.dta$", id.vars = c())

# Assemble df of results
lltm_df <- append_stata_files("^results_lltm_[0-9]*.dta$")
lltm_df <- cbind(method = "lltm", lltm_df)
twostage_df <- append_stata_files("^results_twostage_[0-9]*.dta$")
twostage_df <- cbind(method = "twostage", twostage_df)
df <- subset(rbind(lltm_df, twostage_df), ! seed %in% errors$seed)
df[, c("condition", "npersons")] <- list(NULL)
df$method[df$variable == "dev_rasch"] <- "lltm"


# Rasch lr test ----------------------------------------------------------------

df_rasch <- subset(df, method == "lltm" &
                       variable %in% c("dev_invalidation", "dev_rasch") &
                       model == 2)
df_rasch <- dcast(df_rasch, seed + tau + nitems ~ variable)
df_rasch$pars <- df_rasch$nitems + 1 - 6
df_rasch$dif <- with(df_rasch, dev_invalidation - dev_rasch)
df_rasch$pval <- 1 - pchisq(df_rasch$dif, df = df_rasch$pars)
df_rasch$reject <- df_rasch$pval < .05

# Percent rejected lltm Model 2
aggregate(reject ~ tau + nitems, df_rasch, function(x) mean(x)*100)
# Count of non-rejections
aggregate(reject ~ tau + nitems, df_rasch, function(x) sum(!x))


# Plot labels ------------------------------------------------------------------

refactor <- function(x, label_list) {
  factor(label_list[as.character(x)], levels = label_list)
}

method_labels <- c("lltm" = "LLTM", "twostage" = "LLTM-E2S")

nitems_labels <- paste("I =", sort(unique(df$nitems)))
names(nitems_labels) <- sort(unique(df$nitems))

tau_labels <- paste("tau =", sort(unique(df$tau)))
names(tau_labels) <- sort(unique(df$tau))


# Recovery ----------------------------------------------------------------

# List of generating values
genvalues <- c(beta1 = -.5, beta2 = 1, beta3 = .5, beta4 = .5, beta5 = -.5,
               sigma = 1, tau = NA)

# DF for recovery info
# df_recov <-subset(df, grepl("_(est|se)_", variable) & model == 2 & !is.na(value))
df_recov <-subset(df, grepl("_(est|se)_", variable) & !is.na(value))
df_recov <-subset(df_recov, ! grepl("(sigma|tau)sq_se_.*$", variable))
df_recov$par <- sub("_(est|se)_.*$", "", df_recov$variable)
is_variance <- grepl("sq$", df_recov$par)
df_recov$value[is_variance] <- sqrt(df_recov$value[is_variance])
df_recov$par <- sub("sq$", "", df_recov$par)
df_recov$batch <- sub("^.*_(est|se)_", "", df_recov$variable)
df_recov$variable <- ifelse(grepl("^.*_est_.*$", df_recov$variable), "est", "se")
# df_recov <- dcast(df_recov, method + seed + batch + tau + nitems + par ~ variable)
df_recov <- dcast(df_recov, method + seed + batch + model + tau + nitems + par ~ variable)
df_recov$est[!is.na(df_recov$se)] <- -1 * df_recov$est[!is.na(df_recov$se)]
df_recov$gen <- genvalues[df_recov$par]
is_tau <- grepl("^tau$", df_recov$par)
df_recov$gen[is_tau] <- df_recov$tau[is_tau]
df_recov$dif <- df_recov$est - df_recov$gen
df_recov$capture <- abs(df_recov$dif)/df_recov$se < 1

# Plot formatting
df_recov$Method <- refactor(df_recov$method, method_labels)
par_colors <- brewer.pal(7, "Set1")
par_labels <- list()
for(i in 1:5) par_labels[[i]] <- bquote(beta[.(i)])
par_labels[[6]] <- bquote(sigma)
par_labels[[7]] <- bquote(tau)

# Point estimates
# p1 <- ggplot(subset(df_recov, nitems == 32 & model == 2)) +
#   aes(x = factor(tau), y = dif, color = par) +
#   geom_boxplot() +
#   facet_grid(. ~ Method) +
#   labs(list(x = expression(paste("Residual item SD (", tau, ")")),
#             y = "Difference",
#             color = NULL)) +
#   scale_colour_manual(values = par_colors, labels = par_labels)
# p2 <- ggplot(subset(df_recov, tau == .5 & model == 2)) +
#   aes(x = factor(nitems), y = dif, color = par) +
#   geom_boxplot() +
#   facet_grid(. ~ Method) +
#   labs(list(x = expression(paste("Number of items (", italic(I), ")")),
#             y = "Difference",
#             color = NULL)) +
#   scale_colour_manual(values = par_colors, labels = par_labels)
# png("../figs/2-recovery.png", family = "serif", width = 6, height = 7,
#     units = "in", res = 200, pointsize = 9)
# multiplot(p1, p2, cols=1)
# dev.off()

# Bias
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
p1 <- ggplot(subset(df_bias, nitems == 32 & model == 2)) +
  aes(x = as.factor(tau), y = mean, ymin = lower, ymax = upper, 
      color = factor(par)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(. ~ Method) +
  labs(list(x = expression(paste("Residual item SD (", tau, ")")),
            y = "Bias", color = NULL)) +
  scale_colour_manual(values = par_colors, labels = par_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
p2 <- ggplot(subset(df_bias, tau == .5 & model == 2)) +
  aes(x = as.factor(nitems), y = mean, ymin = lower, ymax = upper, 
      color = factor(par)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(. ~ Method) +
  labs(list(x = expression(paste("Number of items (", italic(I), ")")),
            y = "Bias", color = NULL)) +
  scale_colour_manual(values = par_colors, labels = par_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
png("../figs/2-bias.png", family = "serif", width = 6, height = 7,
    units = "in", res = 200, pointsize = 9)
multiplot(p1, p2, cols=1)
dev.off()

# # Intervals/standard errors
# df_recov_agg <- aggregate(capture ~ Method + model + tau + nitems + par, df_recov,
#                           function(x) mean(x)*100)
# p1 <- ggplot(subset(df_recov_agg, nitems == 32 & model == 2)) +
#   aes(x = as.factor(tau), y = capture, fill = factor(par)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_hline(yintercept = 68, lty = "dashed") +
#   facet_grid(. ~ Method) +
#   labs(list(x = expression(paste("Residual item SD (", tau, ")")),
#             y = "Percentage of replications",
#             fill = NULL)) +
#   scale_fill_manual(values = par_colors[1:5], labels = par_labels[1:5]) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())
# p2 <- ggplot(subset(df_recov_agg, tau == .5 & model == 2)) +
#   aes(x = as.factor(nitems), y = capture, fill = factor(par)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_hline(yintercept = 68, lty = "dashed") +
#   facet_grid(. ~ Method) +
#   labs(list(x = expression(paste("Number of items (", italic(I), ")")),
#             y = "Percentage of replications",
#             fill = NULL)) +
#   scale_fill_manual(values = par_colors[1:5], labels = par_labels[1:5]) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())
# png("../figs/2-capture.png", family = "serif", width = 6, height = 7,
#     units = "in", res = 200, pointsize = 9)
# multiplot(p1, p2, cols=1)
# dev.off()

# Q-Q plots for standard errors
df_recov$z <- df_recov$est / df_recov$se
p1 <- ggplot(subset(df_recov, par == "beta6" & nitems == 32)) +
  aes(sample = z) +
  geom_point(stat = "qq", size = .5) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Normal theoretical quantiles") + ylab("Observed quantiles") +
  facet_grid(tau + nitems ~ Method, labeller = 
             label_bquote(rows = list(tau == .(tau), italic(I) == .(nitems)) ))
p2 <- ggplot(subset(df_recov, par == "beta6" & tau == .5)) +
  aes(sample = z) +
  geom_point(stat = "qq", size = .5) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Normal theoretical quantiles") + ylab("Observed quantiles") +
  facet_grid(tau + nitems ~ Method, labeller = 
             label_bquote(rows = list(tau == .(tau), italic(I) == .(nitems)) ))
png("../figs/2-qqplot.png", family = "serif", width = 6, height = 4,
    units = "in", res = 200, pointsize = 9)
multiplot(p1, p2, cols=2)
dev.off()


# Determine selected models by deviance, aic, bic, lr test ---------------------

df_sel <- subset(df, grepl("^(aic|bic|dev_loo|dev_validation_(new|same)items)",
                           variable))
df_min <- aggregate(value ~ method + seed + variable, data = df_sel, FUN = min)
names(df_min)[names(df_min) == "value"] <- "minimum"
df_sel <- merge(df_sel, df_min)
df_sel$selected <- df_sel$value == df_sel$minimum
df_sel[, c("value", "minimum")] <-list(NULL)
names(df_sel)[names(df_sel) == "variable"] <- "selector"
df_sel$selector <- gsub("^dev_validation_", "", df_sel$selector)
df_sel$basis <- gsub("^(aic|bic|chi2|dev_loo)", "validation", df_sel$selector)


# Append selected models by likelihood ratio test

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
chi2_df <- subset(df, variable == "dev_invalidation")
# Ensure is ordered by model
chi2_df <- chi2_df[with(chi2_df, order(seed, method, model)), ]
chi2_df <- aggregate(value ~ method + seed + tau + nitems,
                     data = chi2_df, FUN = chi2_select)

# Reformat chi2 results to match df_sel and append
chi2_df <- cbind(chi2_df[, -ncol(chi2_df)], chi2_df$value)
chi2_df <- melt(chi2_df, id.vars = c("method", "seed", "tau", "nitems"),
                variable.name = "model", value.name = "selected")
chi2_df$model <- as.numeric(chi2_df$model)
chi2_df$selector <- "chi2"
chi2_df$basis <- "twosample"
chi2_df <- chi2_df[, names(df_sel)]
df_sel <- rbind(df_sel, chi2_df)


# Plots of selection -----------------------------------------------------------

percents_df <- df_sel
percents_df$seed <- NULL
percents_df <- aggregate(selected ~ ., percents_df, function(x) mean(x)*100)

selector_labels <- c("chi2" = "LR test",
                     "aic" = "AIC",
                     "bic" = "BIC",
                     "sameitems" = "HV same items",
                     "newitems" = "HV new items",
                     "dev_loo" = "LOCO-CV")
percents_df$Selector <- refactor(percents_df$selector, selector_labels)
percents_df$Method <- refactor(percents_df$method, method_labels)

ggplot(subset(percents_df, percents_df$nitems == 32)) +
  aes(as.factor(tau), selected, fill = factor(model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(list(x = expression(paste("Residual item SD (", tau, ")")),
            y = "Percentage of times selected", fill = "Model")) +
  facet_grid(Selector ~ Method)
ggsave("../figs/2-selection-tau.png", family = "serif", width = 6, height = 7,
       units = "in", pointsize = 9)

ggplot(subset(percents_df, percents_df$tau == .5)) +
  aes(factor(nitems), selected, fill = factor(model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(list(x = expression(paste("Number of items (", italic(I), ")")),
            y = "Percentage of times selected", fill = "Model")) +
  facet_grid(Selector ~ Method)
ggsave("../figs/2-selection-nitems.png", family = "serif", width = 6, height = 7,
       units = "in", pointsize = 9)


# Penalty -------

df_pen <- subset(df, variable %in% c("dev_invalidation",
                                     "dev_intest",
                                     "dev_test_newitems",
                                     "dev_loo"))
df_pen <- dcast(df_pen, ... ~ variable, value.var = "value")
df_pen$hv <- with(df_pen, dev_test_newitems - dev_intest)
df_pen$loo <- with(df_pen, dev_loo - dev_invalidation)

df_pen <- melt(df_pen, id.vars = c("method", "model", "tau", "nitems"),
                     measure.vars = c("hv", "loo"), na.rm = TRUE)
df_pen$method <- as.character(df_pen$method)
df_pen$method[df_pen$variable == "loo"] <- "loo"

mean_and_ci <- function(x) {
 m <- mean(x)
 s <- sd(x)
 return(c(mean = m, lower = m-s, upper = m+s))
}
df_pen_means <- aggregate(value ~ ., data = df_pen, mean_and_ci)
df_pen_means <- cbind(df_pen_means[, -ncol(df_pen_means)], 
                      df_pen_means[, ncol(df_pen_means)])
label <- c("lltm" = "LLTM HV new items", "twostage" = "LLTM-E2S HV new items",
           "loo" = "LLTM-E2S LOCO-CV")
df_pen_means$Method <- refactor(df_pen_means$method, label)
df_pen_means$text <- ifelse(df_pen_means$nitems == 64 | df_pen_means$tau == 1,
                            paste0("M", df_pen_means$model), "")

aic_penalty <- 2*(5:7)
df_pen_means$aic_penalty <- aic_penalty[df_pen_means$model]

p1 <- ggplot(subset(df_pen_means, nitems == 32)) +
  aes(x = tau, y = mean, color = as.factor(model), label = text) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_line(aes(y = aic_penalty), size = .5, linetype = 2, show.legend = FALSE) +
  geom_text(hjust = "left", nudge_x = (1-.2)*.015, show.legend = FALSE) +
  expand_limits(x = 1 + (1-.2)*.1) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  xlab(expression(paste("Residual item SD (", tau, ")"))) + ylab("Penalty")
p2 <- ggplot(subset(df_pen_means, tau == .5)) +
  aes(x = nitems, y = mean, color = as.factor(model), label = text) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_line(aes(y = aic_penalty), size = .5, linetype = 2, show.legend = FALSE) +
  geom_text(hjust = "left", nudge_x = (64-16)*.015, show.legend = FALSE) +
  expand_limits(x = 64 + (64-16)*.1) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  xlab(expression(paste("Number of items (", italic(I), ")"))) + ylab("Penalty")
png("../figs/2-penalty.png", family = "serif", width = 6, height = 7,
    units = "in", res = 200, pointsize = 9)
multiplot(p1, p2, cols=2)
dev.off()

ggplot(subset(df_pen_means, nitems == 32 & model == 2)) +
  aes(x = tau, y = mean, ymin = lower, ymax = upper, 
      color = as.factor(model), fill = as.factor(model), label = text) +
  geom_point(show.legend = FALSE) +
  geom_ribbon(alpha = .25, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_line(aes(y = aic_penalty), size = .5, linetype = 2, show.legend = FALSE) +
  geom_text(hjust = "left", nudge_x = (1-.2)*.015, show.legend = FALSE) +
  expand_limits(x = 1 + (1-.2)*.1) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  xlab(expression(paste("Residual item SD (", tau, ")"))) + ylab("Penalty")


# ? --------------------
# 
# df_rel <- subset(df, variable %in% c("dev_test_newitems", "dev_loo", "aic") &
#                    method == "twostage")
# df_rel_adj <- aggregate(value ~ method + seed + tau + nitems + variable,
#                         df_rel, mean)
# names(df_rel_adj)[names(df_rel_adj) == "value"] <- "mean"
# df_rel <- merge(df_rel, df_rel_adj)
# df_rel$adj <- with(df_rel, value - mean)
# ggplot(subset(df_rel, nitems == 32)) +
#   aes(x = factor(tau), y = adj, fill = factor(model)) +
#   geom_violin() +
#   facet_wrap(~variable)
# ggplot(subset(df_rel, nitems == 32)) +
#   aes(x = factor(tau), y = adj, fill = factor(model)) +
#   geom_boxplot() +
#   facet_wrap(~variable)
# ggplot(subset(df_rel, tau == .5)) +
#   aes(x = factor(nitems), y = adj, fill = factor(model)) +
#   geom_violin() +
#   facet_wrap(~variable)
# ggplot(subset(df_rel, tau == .5)) +
#   aes(x = factor(nitems), y = adj, fill = factor(model)) +
#   geom_boxplot(coef = 100) +
#   facet_wrap(~variable)
# 
# 
# 
# 
# 
# df_pen$aic_penalty <- aic_penalty[df_pen$model]
# ggplot(subset(df_pen, method != "lltm" & nitems == 32)) +
#   aes(x= factor(tau), y = value, fill = factor(model)) +
#   geom_hline(aes(yintercept = aic_penalty,  color = factor(model)),
#             size = .5, linetype = 2, show.legend = FALSE) +
#   geom_violin(scale = "width") +
#   facet_wrap(~variable)
# ggplot(subset(df_pen, method != "lltm" & tau == .5)) +
#   aes(x= factor(nitems), y = value, fill = factor(model)) +
#   geom_hline(aes(yintercept = aic_penalty,  color = factor(model)),
#              size = .5, linetype = 2, show.legend = FALSE) +
#   geom_violin(scale = "width") +
#   facet_wrap(~variable)
# #
# 
# df_pen_p <- aggregate(value ~ ., df_pen, mean)
# names(df_pen_p)[names(df_pen_p) == "value"] <- "mean"
# df_pen_p <- merge(df_pen, df_pen_p)
# df_pen_p$adj <- with(df_pen_p, value - mean)
# df_pen_p$aic_penalty <- (aic_penalty-mean(aic_penalty))[df_pen_p$model]
# ggplot(subset(df_pen_p, method != "lltm" & nitems == 32)) +
#   aes(x = factor(tau), y = adj, fill = factor(model)) +
#   geom_hline(aes(yintercept = aic_penalty,  color = factor(model)),
#              size = .5, linetype = 2, show.legend = FALSE) +
#   geom_violin(scale = "area") +
#   facet_wrap(~variable)
# ggplot(subset(df_pen_p, method != "lltm" & tau == .5)) +
#   aes(x = factor(nitems), y = adj, fill = factor(model)) +
#   geom_hline(aes(yintercept = aic_penalty,  color = factor(model)),
#              size = .5, linetype = 2, show.legend = FALSE) +
#   geom_violin(scale = "area") +
#   facet_wrap(~variable)


# Plots of HV prediction -------------------------------------------------------

df_dev <- subset(df, grepl("^dev_test_(new|same)items", variable))
names(df_dev)[names(df_dev) == "value"] <- "deviance"
df_dev$basis <- gsub("dev_test_", "", df_dev$variable)
df_dev$variable <- NULL

df_pred <- merge(subset(df_sel, selected & basis %in% c("newitems", "sameitems")),
                 df_dev)
df_pred[, c("model", "basis", "selected")] <- list(NULL)
df_pred <- dcast(df_pred, ... ~ selector, value.var = "deviance")
df_pred$dif <- df_pred$newitems - df_pred$sameitems

# df_agree <- subset(df_sel, selected & basis %in% c("newitems", "sameitems"))
# df_agree[, c("basis", "selected")] <- list(NULL)
# df_agree <- dcast(df_agree, ... ~ selector, value.var = "model")
# df_agree$agree <- ifelse(df_agree$newitems == df_agree$sameitems,
#                          "Agree", "Disagree")
# df_agree[, c("newitems", "sameitems")] <- list(NULL)
# for(n in c("newitems", "sameitems")) {
#   names(df_agree)[names(df_agree) == n] <- paste0("m_", n)
# }

# df_pred <- merge(df_agree, df_pred)
df_pred$Method <- refactor(df_pred$method, method_labels)

p1 <- ggplot(subset(df_pred, nitems == 32 & method == "twostage")) +
  aes(dif) +
  geom_vline(xintercept = 0) +
  geom_density(alpha = .2) +
  facet_grid(.~tau + nitems, labeller =
               label_bquote(cols = list(tau == .(tau), italic(I) == .(nitems)) )) +
  xlab("Difference in deviance (new item HV minus same item HV)") +
  ylab("Density") +
  labs(color = NULL, fill = NULL)
p2 <- ggplot(subset(df_pred, tau == .5 & method == "twostage")) +
  aes(dif) +
  geom_vline(xintercept = 0) +
  geom_density(alpha = .2) +
  facet_grid(.~tau + nitems, labeller =
               label_bquote(cols = list(tau == .(tau), italic(I) == .(nitems)) )) +
  xlab("Difference in deviance (new item HV minus same item HV)") +
  ylab("Density") +
  labs(color = NULL, fill = NULL)
png("../figs/2-hvpred-twostage.png", family = "serif", width = 6, height = 5,
    units = "in", res = 200, pointsize = 9)
multiplot(p1, p2, cols=1)
dev.off()










# p1 <- ggplot(subset(df_pred, nitems == 32 & method == "twostage")) +
#   aes(x = newitems, y = sameitems) +
#   coord_fixed() +
#   geom_point(alpha = .05) +
#   geom_abline(slope = 1, intercept = 0) +
#   facet_grid(.~Tau) +
#   labs(list(x = "Deviance for HV with new items",
#             y = "Deviance for HV with same items"))
# p2 <- ggplot(subset(df_pred, tau == .5 & method == "twostage")) +
#   aes(x = newitems, y = sameitems) +
#   coord_fixed() +
#   geom_point(alpha = .05) +
#   geom_abline(slope = 1, intercept = 0) +
#   facet_grid(.~Nitems) +
#   labs(list(x = "Deviance for HV with new items",
#             y = "Deviance for HV with same items"))
# png("../figs/2-hvpred-twostage.png", family = "serif", width = 6, height = 5,
#     units = "in", res = 200, pointsize = 9)
# multiplot(p1, p2, cols=1)
# dev.off()
#
# ggplot(subset(df_pred, nitems == 64 & method == "twostage")) +
#   aes(x = newitems, y = sameitems) +
#   coord_fixed() +
#   geom_point(alpha = .2) +
#   geom_abline(slope = 1, intercept = 0) +
#   facet_grid(m_newitems~m_sameitems, switch = "both") +
#   labs(list(x = "Deviance for HV with new items",
#             y = "Deviance for HV with same items"))












# ggplot(subset(df_pred, nitems == 32 & method == "twostage")) +
#   aes(dif, fill = agree, color = agree) +
#   geom_vline(xintercept = 0) +
#   geom_density(alpha = .1, show.legend = FALSE) +
#   facet_grid(~Tau) +
#   xlab("Difference in deviance")
# ggsave("../figs/2-hvpred-twostage-left.png", family = "serif", width = 3-.4,
#        height = 4, units = "in", pointsize = 9)
#
# ggplot(subset(df_pred, tau == .5 & method == "twostage")) +
#   aes(dif, fill = agree, color = agree) +
#   geom_vline(xintercept = 0) +
#   geom_density(alpha = .1) +
#   facet_grid(~Nitems) +
#   xlab("Difference in deviance") +
#   labs(color = NULL, fill = NULL)
# ggsave("../figs/2-hvpred-twostage-left.png", family = "serif", width = 3+.4,
#   height = 4, units = "in", pointsize = 9)


# Plots of unvalidated prediction ----------------------------------------------

# df_hvm <- dcast(subset(df_dev, basis == "newitems"),
#                 ... ~ model, value.var = "deviance")
# df_hvm$dif_1 <- df_hvm[, "2"] - df_hvm[, "1"]
# df_hvm$dif_3 <- df_hvm[, "2"] - df_hvm[, "3"]
# df_hvm <- melt(df_hvm, id.vars = c("method", "seed", "tau", "nitems"),
#                measure.vars = c("dif_1", "dif_3"))
#
# df_hvm$Variable <- refactor(df_hvm$variable,
#                             c(dif_1 = "Model 2 vs 1", dif_3 = "Model 2 vs 3"))
# df_hvm$Method <- refactor(df_hvm$method, method_labels)
# df_hvm$Tau <- refactor(df_hvm$tau, tau_labels)
# df_hvm$Nitems <- refactor(df_hvm$nitems, nitems_labels)
#
# p1 <- ggplot(subset(df_hvm, nitems == 32 & method == "twostage")) +
#   aes(value, fill = Variable, color = Variable) +
#   geom_vline(xintercept = 0) +
#   geom_density(alpha = .2) +
#   facet_grid(.~Tau) +
#   xlim(c(-30, 30)) +
#   xlab("Difference in deviance") + ylab("Density") +
#   labs(color = NULL, fill = NULL)
# p2 <- ggplot(subset(df_hvm, tau == .5 & method == "twostage")) +
#   aes(value, fill = Variable, color = Variable) +
#   geom_vline(xintercept = 0) +
#   geom_density(alpha = .2) +
#   facet_grid(. ~ Nitems) +
#   xlim(c(-30, 30)) +
#   xlab("Difference in deviance") + ylab("Density") +
#   labs(color = NULL, fill = NULL)
# png("../figs/2-nohvpred-twostage.png", family = "serif", width = 6, height = 5,
#     units = "in", res = 200, pointsize = 9)
# multiplot(p1, p2, cols=1)
# dev.off()


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
ggsave("../figs/2-rsq-vs-tau.png", family = "serif", width = 3, height = 3,
       units = "in", pointsize = 9)


# Table of item covariates -----------------------------------------------------

X <- subset(ex, person == 1, c(paste0("x", 2:4)))
X <- cbind(x1 = 1, X)
colnames(X) <- paste0("$x_", 1:4, "$")
rownames(X) <- paste("Item", 1:nrow(X))
xtab <- xtable(X, align = "lrrrr", digits = c(0, 0, 2, 0, 10))
print(xtab, file = "../figs/2-table-X.tex", floating = FALSE,
      sanitize.colnames.function = function(x) x)

var(ex$xB) # upsilon^2
