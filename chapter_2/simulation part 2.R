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

# Prints multiple plots showing legend only once
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1,
                                       position = c("bottom", "right")) {

  library(gridExtra)
  library(grid)

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                            legend,
                            ncol = 1,
                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                           legend,
                           ncol = 2,
                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

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

# Bias plots
bias1 <- ggplot(subset(df_bias, nitems == 32 & model == 2)) +
  aes(x = as.factor(tau), y = mean, ymin = lower, ymax = upper,
      color = factor(par)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  facet_grid(. ~ Method) +
  labs(list(x = expression(paste("Residual item SD (", tau, ")")),
            y = "Bias", color = NULL)) +
  scale_colour_manual(values = par_colors, labels = par_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
bias2 <- ggplot(subset(df_bias, tau == .5 & model == 2)) +
  aes(x = as.factor(nitems), y = mean, ymin = lower, ymax = upper,
      color = factor(par)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  facet_grid(. ~ Method) +
  labs(list(x = expression(paste("Number of items (", italic(I), ")")),
            y = "Bias", color = NULL)) +
  scale_colour_manual(values = par_colors, labels = par_labels) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
png("../figs/2-bias.png", family = "serif", width = 6, height = 7,
    units = "in", res = 200, pointsize = 9)
multiplot(bias1, bias2, cols=1)
dev.off()

# Q-Q plots for standard errors
qq1 <- ggplot(subset(df_recov, par == "beta6" & nitems == 32)) +
  aes(sample = z) +
  geom_point(stat = "qq", size = .5) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Normal theoretical quantiles") + ylab("Observed quantiles") +
  facet_grid(tau + nitems ~ Method, labeller =
             label_bquote(rows = list(tau == .(tau), italic(I) == .(nitems)) ))
qq2 <- ggplot(subset(df_recov, par == "beta6" & tau == .5)) +
  aes(sample = z) +
  geom_point(stat = "qq", size = .5) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Normal theoretical quantiles") + ylab("") +
  facet_grid(tau + nitems ~ Method, labeller =
             label_bquote(rows = list(tau == .(tau), italic(I) == .(nitems)) ))
png("../figs/2-qqplot.png", family = "serif", width = 6, height = 4,
    units = "in", res = 200, pointsize = 9)
multiplot(qq1, qq2, cols=2)
dev.off()


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

# Selection plot for varying tau
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

# Selection plot varying I
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
df_penstack <- subset(df_penstack, type == "hv")

# Plot penalties
expand_y_overtau <- range(subset(df_penstack,
                                 method == "twostage" & nitems == 32,
                                 select = c("mean", "lower", "upper")),
                          na.rm = TRUE)
pen1 <- ggplot(subset(df_penstack, nitems == 32)) +
  aes(x = tau, y = mean, ymin = lower, ymax = upper, shape = as.factor(model),
      color = as.factor(model), fill = as.factor(model), label = text) +
  geom_ribbon(alpha = .25, show.legend = FALSE) +
  geom_point(size=2, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_line(aes(y = aic_penalty), size = .5, linetype = 2, show.legend = FALSE) +
  geom_text(hjust = "left", nudge_x = (1-.2)*.015, show.legend = FALSE) +
  expand_limits(x = 1 + (1-.2)*.1, y = expand_y_overtau) +
  facet_grid(Type~side, scale = "free_y") +
  xlab(expression(paste("Residual item SD (", tau, ")"))) + ylab("Penalty") +
  scale_x_continuous(breaks = unique(df_penstack$tau)) +
  theme(panel.grid.minor = element_blank())
expand_y_overnitems <- range(subset(df_penstack,
                                    method == "twostage" & tau == .5,
                                    select = c("mean", "lower", "upper")),
                             na.rm = TRUE)
pen2 <- ggplot(subset(df_penstack, tau == .5)) +
  aes(x = nitems, y = mean, ymin = lower, ymax = upper, shape = as.factor(model),
      color = as.factor(model), fill = as.factor(model), label = text) +
  geom_ribbon(alpha = .25, show.legend = FALSE) +
  geom_point(size=2, show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  geom_line(aes(y = aic_penalty), size = .5, linetype = 2, show.legend = FALSE) +
  geom_text(hjust = "left", nudge_x = (64-16)*.015, show.legend = FALSE) +
  expand_limits(x = 64 + (64-16)*.1, y = expand_y_overnitems) +
  facet_grid(Type~side, scale = "free_y") +
  xlab(expression(paste("Number of items (", italic(I), ")"))) + ylab("Penalty") +
  scale_x_continuous(breaks = unique(df_penstack$nitems)) +
  theme(panel.grid.minor = element_blank())
png("../figs/2-penalty.png", family = "serif", width = 6, height = 7,
    units = "in", res = 200, pointsize = 9)
multiplot(pen1, pen2, cols=1)
dev.off()


# Plots of prediction ----------------------------------------------------------

# Make data frame of relative deviances
df_rel <- subset(df, variable == "dev_test")
df_rel$variable <- NULL
df_mean <- aggregate(value ~ method + seed + tau + nitems, df_rel, mean)
names(df_mean)[names(df_mean) == "value"] <- "mean"
df_rel <- merge(df_rel, df_mean)
df_rel$relative <- df_rel$value - df_rel$mean
df_rel[, c("value", "mean")] <- list(NULL)

# Merge with selection data frame
df_pred <- merge(df_sel, df_rel)
df_pred <- subset(df_pred, selected)
df_predmean <- aggregate(relative ~ selector + tau + nitems + method,
                         df_pred, mean)
df_predmean$Selector <- refactor(df_predmean$selector, selector_labels)
df_predmean$Method <- refactor(df_predmean$method, method_labels)

# Plot prediction means for best models
pred1 <- ggplot(subset(df_predmean, nitems == 32)) +
  aes(x = tau, y = relative, color = Selector, shape = Selector) +
  geom_point(size=2) + geom_line() +
  scale_x_continuous(breaks = unique(df_predmean$tau)) +
  xlab(expression(paste("Residual item SD (", tau, ")"))) +
  ylab("Mean relative holdout deviance") +
  labs(color = NULL, shape = NULL) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  theme(panel.grid.minor = element_blank())
pred2 <- ggplot(subset(df_predmean, tau == .5)) +
  aes(x = nitems, y = relative, color = Selector, shape = Selector) +
  geom_point(size=2) + geom_line() +
  scale_x_continuous(breaks = unique(df_predmean$nitems)) +
  xlab(expression(paste("Number of items (", italic(I), ")"))) +
  ylab("") +
  labs(color = NULL, shape = NULL) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  theme(panel.grid.minor = element_blank())
png("../figs/2-prediction.png", family = "serif", width = 6, height = 7,
    units = "in", res = 200, pointsize = 9)
grid_arrange_shared_legend(pred1, pred2, ncol = 2, nrow = 1, position = "bottom")
dev.off()


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


#  ------

# LR test "penalty"
x <- seq(from = 2, to = 4, by = .05)
plot(x, pchisq(x, 1, lower.tail = FALSE))
abline(h = .05)

# BIC penalties
log(500*32)
log(32)

