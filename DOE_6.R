## Rovsing et al. 2024
## data and scripts for generating linear model and figures for DOE 6

library(pid)
library(ggplot2)
library(dplyr)

## create variables for each factor with DOE coding in standard order
BAFF <- rep(c(1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1), 5)
CD40L_exp <-rep(c(-1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1), 5)
IL4 <- rep(c(-1, -1, +1, +1, +1, +1, -1, -1, +1, +1, +1, +1), 5)
TIME = c(rep(-0.57, 12), rep(-0.29, 12), rep(0, 12), rep(0.43, 12), rep(1, 12))

## create variables for each readout in the same DOE standard order
# IgA
DOE6_IgA <- c(
    rep(0, 12), # Day 3
    rep(0, 12), # Day 5
    c(112, 110,	8,	8,	25,	26,	39,	37,	19,	19,	3,	4), # Day 7
    c(47, 46, 2, 2,	7,	7,	1,	1,	12,	12,	7,	4), # Day 10
    c(383, 388, 60, 61, 165, 165, 40, 41, 326, 328, 295, 336) # Day 14
)
DOE6_IgA.log10 = log10(DOE6_IgA + 1)

# IgE
DOE6_IgE <- c(
    rep(0, 12), # Day 3
    rep(0, 12), # Day 5
    c(0, 0,	0,	0,	10,	8,	0,	0,	7,	7,	0,	0), # Day 7
    c(0, 0,	0,	0,	27,	26,	0,	0,	47,	48,	0,	0), # Day 10
    c(0, 0,	141, 142, 385, 387, 0, 0, 431, 440, 498, 504) # Day 14 
)
DOE6_IgE.log10 = log10(DOE6_IgE + 1)

# IgG1 1:100
DOE6_IgG1 <- c(
    c(10.7266,11.6541,13.7585,6.38411,4.18327,2.51576,1.92445,1.80428,1.83073,15.1546,1.93871,2.60247), # Day 3
    c(15.4502,6.9959,4.08868,3.3141,3.43506,3.5738,2.27526,2.31936,3.15223,3.20009,1.85462,2.00443), # Day 5
    c(100.13,93.0143,3.55612,3.10475,31.2895,30.9403,14.9253,14.708,34.0971,33.7463,3.16519,3.28706), # Day 7
    c(23.6634,25.4505,7.86727,8.34257,96.489,95.5573,2.98214,3.02079,183.215,180.147,8.48203,8.44857), # Day 10
    c(358.446,349.427,597.544,591.628,5146.74,4804.42,57.5783,60.9535,9179.09,10000,1216.71,1297.78) # Day 14 
)
DOE6_IgG1.log10 = log10(DOE6_IgG1 + 1)

# IgG3 1:100
DOE6_IgG3 <- c(
    c(0.359351,0.485441,0.436717,0.379196,0.369971,0.286118,0.207308,0.403052,0.395918,0.384882,0.306301,
      0.612845), # Day 3
    c(3.80482,3.72187,0.552133,0.484354,1.62796,1.38634,0.546641,0.856553,1.63077,1.82726,0.743616,0.966518), # Day 5
    c(294.993,290.31,2.53558,2.18229,27.4484,27.0725,50.2394,48.6797,30.0361,30.4222,5.77153,5.82593), # Day 7
    c(170.11,169.84,5.92103,5.8978,65.8348,64.7179,19.7918,19.8272,125.685,128.896,13.831,14.3236), # Day 10
    c(6847.09,10428.7,2676.07,2954.39,34601.9,34601.9,895.032,903.746,38605.4,38605.4,9026.89,9026.89) # Day 14 
)
DOE6_IgG3.log10 = log10(DOE6_IgG3 + 1)

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE6_IgA ~ BAFF * CD40L_exp * IL4 * TIME)
model = lm(DOE6_IgA ~ BAFF + CD40L_exp + IL4 + TIME)
## final models
model_DOE6_IgA = lm(DOE6_IgA.log10 ~ CD40L_exp + TIME)
model_DOE6_IgE = lm(DOE6_IgE.log10 ~ CD40L_exp*IL4 + IL4*TIME)
model_DOE6_IgG1 = lm(DOE6_IgG1.log10 ~ CD40L_exp*TIME + IL4*TIME)
model_DOE6_IgG3 = lm(DOE6_IgG3.log10 ~ CD40L_exp*TIME)

## print the summary of the model fit
summary(model_DOE6_IgA)
summary(model_DOE6_IgE)
summary(model_DOE6_IgG1)
summary(model_DOE6_IgG3)

## plot a pareto plot for the final model
paretoPlot(model_DOE6_IgA)  ## fig S7.B
paretoPlot(model_DOE6_IgE)  ## fig S7.A
paretoPlot(model_DOE6_IgG1)  ## fig S7.C
paretoPlot(model_DOE6_IgG3)  ## fig S7.D

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 6, height = 5, units = "in")

## plot a contourplot to visualise interactions
contourPlot(model_DOE6_IgE, 'IL4', 'CD40L_exp')

## plot model diagnostics
par(mfrow=c(2,4)) ## plots all four plots in one window
plot(model_DOE6_IgA)

## create a function enabling plotting of raw data
plot_data_2levels <- function(DOE_readout_x, DOE_readout_y, variable, legendLabel, xlabel, ylabel) {
    
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
        scale_color_manual(values = c("-1" = "#f9a042", "1" = "#1a3f71"), labels = c("-1", "1"), 
                           name = legendLabel) +
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
plot_data_TIME <- function(DOE_readout_x, DOE_readout_y, variable, legendLabel, xlabel, ylabel) {
    
    data = data.frame(DOE_readout_x, DOE_readout_y, variable)
    
    ## calculate the mean for each readout variable 
    mean_data1 = data %>% 
        filter(variable == -0.57) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data2 = data %>% 
        filter(variable == -0.29) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data3 = data %>% 
        filter(variable == 0) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data4 = data %>% 
        filter(variable == 0.43) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data5 = data %>% 
        filter(variable == 1) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    ggplot(data, aes(x = DOE_readout_x, y = DOE_readout_y, color = factor(variable))) +
        geom_point(size = 3) +
        scale_color_manual(values = c("-0.57" = "#f6e8a7", "-0.29" = "#f9a042", "0" = "#a24131", "0.43" = "#5f85ae",
                                      "1" = "#1a3f71"), 
                           labels = c("3", "5", "7", "10", "14"), name = legendLabel) +
        geom_vline(xintercept = mean_data1$mean_x, color = "#f6e8a7", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data1$mean_y, color = "#f6e8a7", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data2$mean_x, color = "#f9a042", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data2$mean_y, color = "#f9a042", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data3$mean_x, color = "#a24131", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data3$mean_y, color = "#a24131", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data4$mean_x, color = "#5f85ae", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data4$mean_y, color = "#5f85ae", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data5$mean_x, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data5$mean_y, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
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
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgG3.log10, BAFF, 'BAFF', 'IgG1 (log10)', 'IgG3 (log10)')
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgA.log10, BAFF, 'BAFF', 'IgG1 (log10)', 'IgA (log10)')
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgE.log10, BAFF, 'BAFF', 'IgG1 (log10)', 'IgE (log10)')

## fig S7.F
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgG3.log10, IL4, 'IL4', 'IgG1 (log10)', 'IgG3 (log10)')
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgA.log10, IL4, 'IL4', 'IgG1 (log10)', 'IgA (log10)')
## fig 6.E
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgE.log10, IL4, 'IL4', 'IgG1 (log10)', 'IgE (log10)')

## fig S7.G
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgG3.log10, CD40L_exp, 'CD40L_exp', 'IgG1 (log10)', 'IgG3 (log10)')
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgA.log10, CD40L_exp, 'CD40L_exp', 'IgG1 (log10)', 'IgA (log10)')
## fig 6.F
plot_data_2levels(DOE6_IgG1.log10, DOE6_IgE.log10, CD40L_exp, 'CD40L_exp', 'IgG1 (log10)', 'IgE (log10)')

## fig S7.E
plot_data_TIME(DOE6_IgG1.log10, DOE6_IgG3.log10, TIME, 'Time', 'IgG1 (log10)', 'IgG3 (log10)')
## fig S7.H
plot_data_TIME(DOE6_IgG1.log10, DOE6_IgA.log10, TIME, 'Time', 'IgG1 (log10)', 'IgA (log10)')
## fig 6.D
plot_data_TIME(DOE6_IgG1.log10, DOE6_IgE.log10, TIME, 'Time', 'IgG1 (log10)', 'IgE (log10)')

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 3, height = 2.5, units = "in")
