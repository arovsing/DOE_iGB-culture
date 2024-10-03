## Rovsing et al. 2024
## data and scripts for generating linear model and figures for DOE 4

library(pid)
library(ggplot2)
library(dplyr)

## create variables for each factor with DOE coding in standard order
## here the design is full factorial repeated with either IL21 or CM being constant
BAFF = c(-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1)
CD40L_exp = c(-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1)
IL4_pulse = c(-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1)
IL21 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1)
CM = c(-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

## create variables for each readout in the same DOE standard order
DOE4_proliferation = c(39812.5,69492.5,86800,59187.5,40530,24300,2225,3050,32337.5,34662.5,49375,51812.5,13250,17050,
                       1387.5,1050,40900,33950,52225,35950,10362.5,8237.5,1362.5,2500,49012.5,48400,35312.5,43812.5,
                       12050,10562.5,2537.5,1450)

DOE4_viability = c(67.2,76.6,77.6,68.6,61.0,51.2,36.1,36.1,62.2,70.6,77.2,71.8,41.8,43.6,31.0,29.3,67.4,59.5,73.0,63.7,
               39.3,34.9,28.9,34.1,74.7,68.9,73.7,67.7,33.5,38.7,33.7,30.2)

## Europium IgG ng/mL
DOE4_IgG = c(7.00249,9.22255,3.3762,6.64173,4.38132,6.89412,0.773492,0.810687,6.6547,8.36994,4.18787,4.4943,5.26207,
               5.08436,1.07477,0.769382,5.0045,8.01301,3.57891,6.11264,4.11957,4.93582,0.559887,0.652444,5.32629,
               7.42495,3.24068,3.40281,5.28322,5.04366,0.470296,0.614788)

## Europium IgE mU
DOE4_IgE = c(0,0,0,0,9.29197,24.7627,0,0,0,0,0,0,15.0322,32.2263,0,0,0,0,0,0,16.7669,29.0631,0,0,0,0.667206,0.600013,0,
             22.5999,26.4129,0,0)

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE4_proliferation ~ BAFF * CD40L_exp * IL4_pulse * IL21 * CM)
model = lm(DOE4_proliferation ~ BAFF + CD40L_exp + IL4_pulse + IL21 + CM)
## final model
model_DOE4_proliferation = lm(DOE4_proliferation ~ BAFF*CD40L_exp*IL4_pulse*CM)

## final model
model_DOE4_viability = lm(DOE4_viability ~ CD40L_exp*IL4_pulse + CM)

## final model
model_DOE4_IgE = lm(DOE4_IgE ~ CM + BAFF*IL4_pulse + IL4_pulse*CD40L_exp)

## final model
model_DOE4_IgG = lm(DOE4_IgG ~ BAFF*CD40L_exp*IL4_pulse + CD40L_exp*IL4_pulse*CM + BAFF*IL4_pulse)

## print the summary of the model fit
summary(model_DOE4_proliferation)
summary(model_DOE4_viability)
summary(model_DOE4_IgG)
summary(model_DOE4_IgE)

## plot a pareto plot for the final model
paretoPlot(model_DOE4_proliferation)  ## fig S6.A
paretoPlot(model_DOE4_viability)  ## fig 5.B
paretoPlot(model_DOE4_IgG)  ## fig S6.C
paretoPlot(model_DOE4_IgE)  ## fig 5.H

## plot a contourplot to visualise interactions
contourPlot(model_DOE4_proliferation, 'IL4_pulse', 'CD40L_exp')

## plot model diagnostics
par(mfrow=c(2,4)) ## plots all four plots in one window
plot(model_DOE4_proliferation)

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
plot_data_3levels <- function(DOE_readout_x, DOE_readout_y, variable, legendLabel, xlabel, ylabel) {
    
    data = data.frame(DOE_readout_x, DOE_readout_y, variable)
    
    ## calculate the mean for each readout variable 
    mean_data1 = data %>% 
        filter(variable == -1) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data2 = data %>% 
        filter(variable == 1) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    mean_data3 = data %>% 
        filter(variable == 0) %>% 
        summarize(mean_x = mean(DOE_readout_x), mean_y = mean(DOE_readout_y))
    
    ggplot(data, aes(x = DOE_readout_x, y = DOE_readout_y, color = factor(variable))) +
        geom_point(size = 3) +
        scale_color_manual(values = c("-1" = "#f9a042", "0" = "#f6e8a7", "1" = "#1a3f71"), 
                           labels = c("-1", "0", "1"), name = legendLabel) +
        geom_vline(xintercept = mean_data1$mean_x, color = "#f9a042", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data1$mean_y, color = "#f9a042", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data2$mean_x, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data2$mean_y, color = "#1a3f71", linetype = "dashed", linewidth = 1) +
        geom_vline(xintercept = mean_data3$mean_x, color = "#f6e8a7", linetype = "dashed", linewidth = 1) +
        geom_hline(yintercept = mean_data3$mean_y, color = "#f6e8a7", linetype = "dashed", linewidth = 1) +
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
plot_data_2levels(DOE4_proliferation, DOE4_viability, BAFF, 'BAFF', 'Number of cells', 'Cell viability (%)')
## fig 5.C
plot_data_2levels(DOE4_proliferation, DOE4_viability, IL4_pulse, 
                  'Pulse of IL4', 'Number of cells', 'Cell viability (%)')
## fig S6.E
plot_data_3levels(DOE4_proliferation, DOE4_viability, IL21, 'IL21', 'Number of cells', 'Cell viability (%)')
## fig 5.D
plot_data_2levels(DOE4_proliferation, DOE4_viability, CD40L_exp, 'CD40L_exp', 'Number of cells', 'Cell viability (%)')
plot_data_2levels(DOE4_proliferation, DOE4_viability, CM, 'Change medium', 'Number of cells', 'Cell viability (%)')

plot_data_2levels(DOE4_IgG, DOE4_IgE, BAFF, 'BAFF', 'IgG (ng/mL)', 'IgE (mU)')
## fig 5.I
plot_data_2levels(DOE4_IgG, DOE4_IgE, IL4_pulse, 'Pulse of IL4', 'IgG (ng/mL)', 'IgE (mU)')
## fig S6.F
plot_data_3levels(DOE4_IgG, DOE4_IgE, IL21, 'IL21', 'IgG (ng/mL)', 'IgE (mU)')
## fig 5.J
plot_data_2levels(DOE4_IgG, DOE4_IgE, CD40L_exp, 'CD40L_exp', 'IgG (ng/mL)', 'IgE (mU)')
plot_data_2levels(DOE4_IgG, DOE4_IgE, CM, 'Change medium', 'IgG (ng/mL)', 'IgE (mU)')

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 3, height = 2.5, units = "in")

## find maximum
int = seq(-1,1,by=0.05)
BAFF = CD40L_exp = IL4_pulse = CM = int
design = expand.grid(BAFF = BAFF, CD40L_exp = CD40L_exp, IL4_pulse = IL4_pulse, CM = CM)
df.find_max = predict(model_DOE4_proliferation, design)
max(df.find_max)
design[which.max(df.find_max), ]
