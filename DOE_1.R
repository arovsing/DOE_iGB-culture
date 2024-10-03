## Rovsing et al. 2024
## data and scripts for generating linear model and figures for DOE 1

library(pid)
library(ggplot2)
library(dplyr)

## create variables for each factor with DOE coding in standard order
## using 'expand.grid' for standard DOE full factorial design
A = B = C = D = c(-1, 1)
design = expand.grid(A=A, B=B, C=C, D=D)
A = design$A
B = design$B
C = design$C
D = design$D

BAFF = c(A, 1, -1, -1, 1)
CD40L = c(B, 1, -1, 1, -1)
IL4 = c(C, 1, -1, 1, -1)
IL21 = c(D, 1, -1, 1, -1)

## create variables for each readout in the same DOE standard order
DOE1_proliferation = c(33550, 64200, 160850, 160750, 30750, 31000, 87600, 75300,
                       49250, 58950, 142200, 172750, 23100, 22750, 108950, 
                       76600, 118000, 61700, 88150, 50950)

DOE1_viability = c(90.2, 90.4, 94.9, 96.8, 74.1, 79.8, 84.8, 85.2, 91.6, 90.1, 
                   93.1, 95.6, 78.4, 75.6, 78.0, 79.8, 80.0, 82.5, 82.1, 86.4)

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE1_proliferation ~ BAFF * CD40L * IL4 * IL21)
model = lm(DOE1_proliferation ~ BAFF + CD40L + IL4 + IL21)
## final model
model_DOE1_proliferation = lm(DOE1_proliferation ~ CD40L*IL4)

model = lm(DOE1_viability ~ BAFF * CD40L * IL4 * IL21)
model = lm(DOE1_viability ~ BAFF + CD40L + IL4 + IL21)
## final model
model_DOE1_viability = lm(DOE1_viability ~ CD40L + IL4)

## print the summary of the model fit
summary(model_DOE1_proliferation)
summary(model_DOE1_viability)

## plot pareto plots for the final model
paretoPlot(model_DOE1_proliferation)  ## fig 2.E
paretoPlot(model_DOE1_viability)  ## fig 2.D

## plot a contourplot to visualise interactions
contourPlot(model_DOE1_proliferation, 'IL4', 'CD40L') ## fig 2.G
contourPlot(model_DOE1_viability, 'IL4', 'CD40L') ## fig 2.F

## plot model diagnostics
par(mfrow=c(2,4)) ## plots all four plots in one window
plot(model)

## create a function enabling plotting of raw data
plot_data <- function(DOE_readout_x, DOE_readout_y, variable, legendLabel, xlabel, ylabel) {
    
    data = data.frame(DOE_readout_x, DOE_readout_y, variable)
    
    ## calculate the mean for each readout variable 
    mean_data1 = data %>% 
        filter(variable == -1) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data2 = data %>% 
        filter(variable == 1) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    ggplot(data, aes(x = DOE_readout_x, y = DOE_readout_y, color = factor(variable))) +
        geom_point(size = 3) +
        scale_color_manual(values = c("-1" = "#f9a042", "1" = "#1a3f71"), labels = c("-1", "1"), name = legendLabel) +
        geom_vline(xintercept = mean_data1$mean_x, color = "#f9a042", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data1$mean_y, color = "#f9a042", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data2$mean_x, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data2$mean_y, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
        labs(x = xlabel, y = ylabel) +
        theme_classic() +
        theme(
            text = element_text(size = 8, color = "#221e20"),
            axis.title = element_text(size = 8, color = "#221e20"),
            axis.text = element_text(size = 8, color = "#221e20"),
            legend.title = element_text(size = 8, color = "#221e20"),
            legend.text = element_text(size = 8, color = "#221e20")
        )
}

## plot raw data for each variable
## fig 2.M
plot_data(DOE1_proliferation, DOE1_viability, BAFF, 'BAFF', 'Number of cells', 'Cell viability (%)')
## fig 2.J
plot_data(DOE1_proliferation, DOE1_viability, IL4, 'IL4', 'Number of cells', 'Cell viability (%)')
## fig 2.L
plot_data(DOE1_proliferation, DOE1_viability, IL21, 'IL21', 'Number of cells', 'Cell viability (%)')
## fig 2.K
plot_data(DOE1_proliferation, DOE1_viability, CD40L, 'CD40L', 'Number of cells', 'Cell viability (%)')

## save plot
#ggsave(filename, plot = last_plot(), width = 2.5, height = 2.5, units = "in")