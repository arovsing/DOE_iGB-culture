## Rovsing et al. 2025
## data and scripts for generating linear model and figures for DOE 6 IgM

library(pid)
library(ggplot2)
library(dplyr)

## create variables for each factor with DOE coding in standard order
BAFF <- rep(c(1, 1, 1, 1, -1, -1), 5)
CD40L_exp <-rep(c(-1, +1, -1, +1, -1, +1), 5)
IL4 <- rep(c(-1, +1, +1, -1, +1, +1), 5)
TIME = c(rep(-0.57, 6), rep(-0.29, 6), rep(0, 6), rep(0.43, 6), rep(1, 6))

## create variables for each readout in the same DOE standard order
# IgM
DOE6_IgM <- c(
    c(2, 1, 0, 0, 1, 1), # Day 3
    c(299, 7, 43, 65, 50, 8), # Day 5
    c(10937, 73, 461, 2979, 535, 64), # Day 7
    c(2432, 125, 445, 642, 903, 189), # Day 10
    c(17682, 22272, 17848, 12000, 26248, 37866) # Day 14
)
DOE6_IgM.log10 = log10(DOE6_IgM + 1)

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE6_IgM.log10 ~ BAFF * CD40L_exp * IL4 * TIME)
model = lm(DOE6_IgM.log10 ~ BAFF + CD40L_exp + IL4 + TIME)
## final models
model_DOE6_IgM = lm(DOE6_IgM.log10 ~ IL4 + TIME)

## print the summary of the model fit
summary(model_DOE6_IgM)

## plot a pareto plot for the final model
paretoPlot(model_DOE6_IgM)  ## fig S7.A

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 6, height = 5, units = "in")

## plot model diagnostics
par(mfrow=c(2,4)) ## plots all four plots in one window
plot(model_DOE6_IgM)

## plot raw data for each variable
data <- data.frame(
    DOE6_IgM.log10 = DOE6_IgM.log10,
    TIME            = TIME,
    IL4             = IL4
)

## compute group means for IL4 = –1 and IL4 = 1
mean_IL4_neg <- data %>%
    filter(IL4 == -1) %>%
    summarize(mean_value = mean(DOE6_IgM.log10))

mean_IL4_pos <- data %>%
    filter(IL4 == 1) %>%
    summarize(mean_value = mean(DOE6_IgM.log10))

## fig S7.E
ggplot(data, aes(x = TIME, y = DOE6_IgM.log10, color = factor(IL4))) +
    geom_point(size = 3) +
    scale_color_manual(
        values = c("-1" = "#f9a042", "1" = "#1a3f71"),
        labels = c("-1" = "0",  "1" = "10"),
        name   = "IL4"
    ) +
    # mean lines
    geom_hline(yintercept = mean_IL4_neg$mean_value, color = "#f9a042", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = mean_IL4_pos$mean_value, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
    # custom x-axis breaks → days
    scale_x_continuous(
        breaks = c(-0.57, -0.29, 0,     0.43, 1),
        labels = c("3",  "5", "7", "10", "14")
    ) +
    labs(
        x = "Time (days)",
        y = "IgM (ng/mL, log10)"
    ) +
    theme_classic() +
    theme(
        text         = element_text(size = 8, color = "#221e20"),
        axis.title   = element_text(size = 8, color = "#221e20"),
        axis.text    = element_text(size = 8, color = "#221e20"),
        legend.title = element_text(size = 8, color = "#221e20"),
        legend.text  = element_text(size = 8, color = "#221e20")
    )

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 3, height = 2.5, units = "in")
