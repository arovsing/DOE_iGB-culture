## Rovsing et al. 2024
## data and scripts for generating linear model and figures for DOE 5

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
DOE5_proliferation = c(79450,81742.5,116830,119105,58607.5,49752.5,5302.5,5250,46550,37047.5,125282.5,37310,18970,32585,
                       1382.5,1400,48020,62195,58187.5,23625,14595,18427.5,3972.5,3115,76002.5,77367.5,71592.5,50242.5,
                       32095,16712.5,3027.5,3307.5)

DOE5_viability = c(90.9,92.2,93.2,89.1,88.1,81.7,55.3,53.4,82.8,84.5,91.5,93.6,72.9,71.3,34.8,46.5,86.8,91.1,93.6,91.1,
                   74.0,71.2,50.8,56.7,91.0,91.4,88.4,88.8,68.8,67.9,41.6,46.3)

## Europium IgG ng/mL
DOE5_IgG = c(3.36835,3.20493,2.20902,2.23878,4.89121,3.33943,0.506452,0.561402,2.75817,4.51315,2.72034,3.4874,2.89818,
             2.85445,0.708451,1.23398,2.62348,3.63928,1.26982,2.02476,2.49858,2.40015,0.392038,0.180397,3.93873,4.29771,
             2.24462,2.64566,3.34202,2.86562,0.691452,0.94529)

## Europium IgE mU
DOE5_IgE = c(0,0.49875,0.00803734,0,12.5029,15.7466,0,0,0,0,0,0.0922662,10.9091,13.3056,0,0,0,0,0,0,9.64682,13.4406,
             0.107534,0,0.565743,0.203281,0,0,21.4853,11.158,0.657414,0)

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE5_proliferation ~ BAFF * CD40L_exp * IL4_pulse * IL21 * CM)
model = lm(DOE5_proliferation ~ BAFF + CD40L_exp + IL4_pulse + IL21 + CM)
## final model
model_DOE5_proliferation = lm(DOE5_proliferation ~ CD40L_exp*IL4_pulse*CM)

## final model
model_DOE5_viability = lm(DOE5_viability ~ CD40L_exp*IL4_pulse + IL4_pulse*CM)

## final model
model_DOE5_IgG = lm(DOE5_IgG ~ IL21 + BAFF*CM + CD40L_exp*IL4_pulse*CM)

## final model
model_DOE5_IgE = lm(DOE5_IgE ~ BAFF*CD40L_exp*IL4_pulse*IL21)

## print the summary of the model fit
summary(model_DOE5_proliferation)
summary(model_DOE5_viability)
summary(model_DOE5_IgG)
summary(model_DOE5_IgE)

## plot a pareto plot for the final model
paretoPlot(model_DOE5_proliferation)  ## fig S6.B
paretoPlot(model_DOE5_viability)  ## fig 5.E
paretoPlot(model_DOE5_IgG)  ## fig S6.D
paretoPlot(model_DOE5_IgE)  ## fig 5.K

## plot a contourplot to visualise interactions
contourPlot(model_DOE5_IgE, 'IL4_pulse', 'CD40L_exp')

## plot model diagnostics
par(mfrow=c(2,4)) ## plots all four plots in one window
plot(model_DOE5_proliferation)

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
plot_data_2levels(DOE5_proliferation, DOE5_viability, BAFF, 'BAFF', 'Number of cells', 'Cell viability (%)')
## fig 5.F
plot_data_2levels(DOE5_proliferation, DOE5_viability, IL4_pulse, 'Pulse of IL4', 'Number of cells', 'Cell viability (%)')
plot_data_3levels(DOE5_proliferation, DOE5_viability, IL21, 'IL21', 'Number of cells', 'Cell viability (%)')
## fig 5.G
plot_data_2levels(DOE5_proliferation, DOE5_viability, CD40L_exp, 'CD40L_exp', 'Number of cells', 'Cell viability (%)')
plot_data_2levels(DOE5_proliferation, DOE5_viability, CM, 'Change medium', 'Number of cells', 'Cell viability (%)')

plot_data_2levels(DOE5_IgG, DOE5_IgE, BAFF, 'BAFF', 'IgG (ng/mL)', 'IgE (mU)')
## fig 5.L
plot_data_2levels(DOE5_IgG, DOE5_IgE, IL4_pulse, 'Pulse of IL4', 'IgG (ng/mL)', 'IgE (mU)')
plot_data_3levels(DOE5_IgG, DOE5_IgE, IL21, 'IL21', 'IgG (ng/mL)', 'IgE (mU)')
## fig 5.M
plot_data_2levels(DOE5_IgG, DOE5_IgE, CD40L_exp, 'CD40L_exp', 'IgG (ng/mL)', 'IgE (mU)')
plot_data_2levels(DOE5_IgG, DOE5_IgE, CM, 'Change medium', 'IgG (ng/mL)', 'IgE (mU)')

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 3, height = 2.5, units = "in")

## find maximum
int = seq(-1,1,by=0.05)
BAFF = CD40L_exp = IL4_pulse = CM = int
design = expand.grid(BAFF = BAFF, CD40L_exp = CD40L_exp, IL4_pulse = IL4_pulse, CM = CM)
df.find_max = predict(model_DOE5_proliferation, design)
max(df.find_max)
design[which.max(df.find_max), ]
