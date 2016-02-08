library(reshape2)
library(ggplot2)
theme_set(theme_bw())

# Adapted from: http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
ggcolors = gg_color_hue(3)



# Plot of r-square against tau -------------------------------------------------

# Get variance of fixed part of items
example <- read.csv("example_sim_data.csv")
example.sub <- unique(example[, c("item", "xB")])
# Using population variance:
upsilon.sq <- (nrow(example.sub)-1) / nrow(example.sub) * var(example.sub$xB)

r.sq <- function(tau, upsilon.sq) upsilon.sq / (tau^2 + upsilon.sq)

tau <- c(0, .1, .3, .5)
pts <- data.frame(tau = tau, r.sq = r.sq(tau, upsilon.sq))

ggplot() +
  coord_cartesian() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(data = data.frame(pts),
             mapping = aes(x = tau, y = r.sq),
             size = 2) +
  stat_function(data = data.frame(x = c(0, .75)),
                mapping = aes(x),
                fun = r.sq,
                args = list(upsilon.sq = upsilon.sq)) +
  xlab(expression(paste("Residual item SD (", tau, ")"))) +
  ylab(expression(paste("Explained variance (", R^2, ")")))
ggsave("figs/rsq_vs_tau.pdf", family = "Times", width = 3, height = 2,
       units = "in", pointsize = 9)


# Set up data frame of simulation results --------------------------------------

# Assemble data frame from individual files
raw.files <- dir(pattern = "^sim1_result_m[0-9]*_[0-9]*.csv$")
raw <- NULL
for(i in 1:length(raw.files)) {
  raw <- rbind(raw, read.csv(raw.files[i]))
}
melted <- melt(raw, id.vars = c("seed", "condition", "model", "tau", "nitems", "npersons"),
               measure.vars = c("insample", "aic", "bic", "newpersons", "newitems", "newboth"))



# Add column for k (n eff parameters) for each fit and "variable"
insample <- subset(melted, variable == "insample",
                   select = c("seed", "model", "value"))
names(insample)[names(insample) == "value"] <- "insample"
melted <- merge(melted, insample)
melted$k <- (melted$value - melted$insample) / 2
melted$insample <- NULL


# Add column for which model was selected
minimums <- aggregate(melted$value,
                      by = list(seed = melted$seed,
                                variable = melted$variable),
                      FUN = min)
names(minimums)[names(minimums) == "x"] <- "min"
melted <- merge(melted, minimums)
melted$selected <- melted$value == melted$min
melted$min <- NULL

# Add columns for indicators for simulation strand
melted$overtau <- melted$condition %in% 1:4
melted$overnitems <- melted$tau == .3

# Add rows for selection via likelihood ratio test
insample.cast <- dcast(subset(melted, variable == "insample"),
                       seed + condition + tau + nitems + npersons + overtau + overnitems ~ model)
insample.cast$p2v1 <- 1 - pchisq(insample.cast[, "1"] - insample.cast[, "2"],
                                 df = 1)
insample.cast$p3v2 <- 1 - pchisq(insample.cast[, "2"] - insample.cast[, "3"],
                                 df = 1)
chisq.cast <- insample.cast[, !(names(insample.cast) %in% as.character(1:3)) ]
chisq.cast[, "1"] = insample.cast$p2v1 >= .05
chisq.cast[, "2"] = insample.cast$p2v1 < .05 & insample.cast$p3v2 >= .05
chisq.cast[, "3"] = insample.cast$p2v1 < .05 & insample.cast$p3v2 < .05
chisq.melt <- melt(chisq.cast,
                   id.vars = c("seed", "condition", "tau", "nitems", "npersons",
                               "overtau", "overnitems"),
                   measure.vars = as.character(1:3),
                   variable.name = "model",
                   value.name = "selected")
chisq.melt$variable <- "lrtest"
chisq.melt$value <- NA
chisq.melt$k <- NA
chisq.melt <- chisq.melt[, names(melted)]
melted <- rbind(melted, chisq.melt)

# Add column for factor variable with long names for each "variable"
labels <- c("lrtest" = "Likelihood ratio test",
            "insample" = "In-sample deviance",
            "aic" = "AIC",
            "bic" = "BIC",
            "newpersons" = "CV over persons",
            "newitems" = "CV over items",
            "newboth" = "CV over persons and items")
melted$labelled.variable <- factor(labels[as.character(melted$variable)],
                                   levels = (labels))

# List of grouping variables for aggregate()
grouping.list <- with(melted, list(tau = tau,
                                   nitems = nitems,
                                   npersons = npersons,
                                   model = model,
                                   variable = variable,
                                   overtau = overtau,
                                   overnitems = overnitems,
                                   labelled.variable = labelled.variable))


# Plots of selection -----------------------------------------------------------

selected <- aggregate(melted$selected, by = grouping.list, FUN = mean)
selected <- subset(selected, variable %in% c("lrtest", "aic", "bic",
                                             "newpersons", "newitems"))

ggplot(subset(selected, overtau)) +
  aes(factor(tau), x, fill = factor(model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(list(x = expression(paste("Residual item SD (", tau, ")")),
            y = "Proportion of times selected", fill = "Model")) +
  facet_wrap(~labelled.variable, nrow = 2)
ggsave("figs/select_overtau.pdf", family = "Times", width = 6, height = 4,
       units = "in", pointsize = 9)

ggplot(subset(selected, overnitems)) +
  aes(factor(nitems), x, fill = factor(model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(list(fill = "Model", x = "Number of items",
            y = "Proportion of times selected")) +
  facet_wrap(~labelled.variable, nrow = 2)
ggsave("figs/select_overnitems.pdf", family = "Times", width = 6, height = 4,
       units = "in", pointsize = 9)


# Plots of k -------------------------------------------------------------------

penalty <- aggregate(melted$k, by = grouping.list, FUN = mean)

ggplot(subset(penalty, variable %in% c("newpersons", "newitems") & overtau)) +
  aes(tau, x, color = factor(model)) +
  geom_point() + geom_line() +
  xlab(expression(paste("Residual item SD (", tau, ")"))) +
  ylab("Eff. N parameters (k)") +
  labs(color = "Model") +
  facet_wrap(~labelled.variable, scales = "free_y")
ggsave("figs/k_overtau.pdf", family = "Times", width = 5, height = 2,
       units = "in", pointsize = 9)

ggplot(subset(penalty, variable %in% c("newpersons", "newitems") & overnitems)) +
  aes(nitems, x, group = factor(model), color = factor(model)) +
  geom_point() + geom_line() +
  xlab("Number of items") +
  ylab("Eff. N parameters (k)") +
  labs(color = "Model") +
  facet_wrap(~labelled.variable, scales = "free_y")
ggsave("figs/k_overnitems.pdf", family = "Times", width = 5, height = 2,
       units = "in", pointsize = 9)


#  -------------------------------------------------------------------

diffs <- dcast(melted[, ! names(melted) %in% c("k", "selected")],
               ... ~ model,
               value.var = "value")
diffs <- subset(diffs, overtau & variable %in% c("aic", "bic",
                                                 "newpersons", "newitems"))
diffs$diff1 <- diffs[, "1"] - diffs[, "2"]
diffs$diff3 <- diffs[, "3"] - diffs[, "2"]

percent <- function(x, gt) {
  number <- mean( (x > 0) == gt) * 100
  string <- paste0(sprintf("%0.1f", number), "%")
  return(string)
}
wins <- aggregate(cbind(diff1, diff3) ~ variable + condition,
                  diffs, percent, gt = FALSE)
names(wins) <- sub("diff", "win", names(wins))
losses <- aggregate(cbind(diff1, diff3) ~ variable + condition,
                    diffs, percent, gt = TRUE)
names(losses) <- sub("diff", "lose", names(losses))
win_lose <- merge(wins, losses)
diffs <- merge(diffs, win_lose)

p1 <- ggplot(diffs) +
  aes(diff1) + #xlim(-200, 200) +
  facet_grid(tau ~ labelled.variable) +
  geom_histogram(binwidth = 20, fill = ggcolors[1]) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  ylab("Count") + xlab("Difference")
p1_build <- ggplot_build(p1)
x1 <- .95 * p1_build$panel$ranges[[1]]["x.range"][[1]]
y1 <- .95 * p1_build$panel$ranges[[1]]["y.range"][[1]]
p1 +
  geom_text(aes(label = win1, x = x1[1], y = y1[2], fontface = "bold"),
            stat = "unique", hjust = "left", vjust = "top", size = 3) +
  geom_text(aes(label = lose1, x = x1[2], y = y1[2], fontface = "bold"),
            stat = "unique", hjust = "right", vjust = "top", size = 3) +
  theme(axis.text=element_text(size=6))
ggsave("figs/twoway_model1.pdf", family = "Times", width = 5, height = 5, units = "in")


p3 <- ggplot(diffs) +
  aes(diff3) + #xlim(-200, 200) +
  facet_grid(tau ~ labelled.variable) +
  geom_histogram(binwidth = 20, fill = ggcolors[3]) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  ylab("Count") + xlab("Difference")
p3_build <- ggplot_build(p3)
x3 <- .95 * p3_build$panel$ranges[[1]]["x.range"][[1]]
y3 <- .95 * p3_build$panel$ranges[[1]]["y.range"][[1]]
p3 +
  geom_text(aes(label = win3, x = x3[1], y = y3[2], fontface = "bold"),
            stat = "unique", hjust = "left", vjust = "top", size = 3) +
  geom_text(aes(label = lose3, x = x3[2], y = y3[2], fontface = "bold"),
            stat = "unique", hjust = "right", vjust = "top", size = 3) +
  theme(axis.text=element_text(size=9))
ggsave("figs/twoway_model3.pdf", family = "Times", width = 5, height = 5, units = "in")




# Write latex macros -----------------------------------------------------------

mycat <- function(name, values, fmt = "%.2f", file = "figs/macros.tex") {
  fmt.values <- sprintf(fmt, sort(values))
  contents <- paste(fmt.values, collapse = ", ")
  if(length(values) > 1) {
    option <- "[1][]"
    contents <- sub(", ([0-9\\.]*)$", ", \\{#1} \\1", contents)
  } else {
    option <- ""
    contents <- contents
  }
  command <- paste0("\\newcommand{\\", name, "}", option, "{", contents, "}")
  cat(command) # show command to be written
  cat(command, "\n", file = file, append = TRUE)
}

if(file.exists("figs/macros.tex")) file.remove("figs/macros.tex")
mycat("genupsilonsq", upsilon.sq)
mycat("gentau", tau)
mycat("gentausq", tau^2)
mycat("rsq", r.sq(tau, upsilon.sq))
mycat("aic", subset(penalty, variable == "aic" & tau == 0)[, "x"], "%.0f")
mycat("bic", subset(penalty, variable == "bic" & tau == 0)[, "x"], "%.2f")
mycat("nreps", length(unique(melted$seed[melted$condition==1])), "%.0f")
mycat("nitems", unique(melted$nitems[melted$overtau]), "%.0f")
mycat("npersons", unique(melted$npersons[melted$overtau]), "%.0f")
mycat("nitemsoveritems", unique(melted$nitems[melted$overnitems]), "%.0f")
mycat("gentauoveritems", unique(melted$tau[melted$overnitems]))



