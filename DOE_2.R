## Rovsing et al. 2024
## data and scripts for generating linear model and figures for DOE 2

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

## design repeated with two levels of CD40L expression
BAFF = c(BAFF, BAFF)
CD40L_no = c(CD40L, CD40L)
IL4 = c(IL4, IL4)
IL21 = c(IL21, IL21)
CD40L_exp = rep(c(-1,1), each = 20)

## create variables for each readout in the same DOE standard order
DOE2_proliferation = c(45250, 64950, 197000, 198850, 31150, 43100, 85700, 73850, 56900, 74850, 142650, 158450, 37800, 
                       43000, 83650, 62500, 92500, 47000, 55900, 90750, 61600, 97250, 136100, 131650, 28100, 21450, 
                       26800, 21900, 111400, 85600, 84800, 109250, 14650, 26450, 24550, 24100, 34800, 88250, 21300, 
                       80300)

DOE2_viability = c(94.4, 92.0, 92.4 , 90.3 , 76.9 , 75.2 , 77.1 , 57.8 , 88.9 , 91.7 , 91.5 , 93.7 , 72.2 , 68.6 , 
                   47.6 , 36.1 , 45.4 , 86.6 , 55.4 , 87.1 , 90.4 , 91.5 , 90.3 , 88.6 , 35.1 , 30.3 , 33.2 , 30.4 , 
                   88.2 , 87.7 , 84.6 , 87.0 , 46.1 , 31.9 , 31.6 , 26.8 , 32.8 , 82.5 , 34.0 , 85.9 )

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE2_proliferation ~ BAFF * CD40L_no * IL4 * IL21 * CD40L_exp)
model = lm(DOE2_proliferation ~ BAFF + CD40L_no + IL4 + IL21 + CD40L_exp)
## final model
model_DOE2_proliferation = lm(DOE2_proliferation ~ CD40L_no + IL4 + IL21 + CD40L_exp + CD40L_no*IL4*IL21 + 
                                  CD40L_no*IL4*CD40L_exp)

## final model
model_DOE2_viability = lm(DOE2_viability ~ BAFF*IL4 + CD40L_no*IL4 + CD40L_no*IL21 + CD40L_no*IL4*CD40L_exp + 
                              IL4*IL21*CD40L_exp)

## print the summary of the model fit
summary(model_DOE2_proliferation)
summary(model_DOE2_viability)

## plot a pareto plot for the final model
paretoPlot(model_DOE2_proliferation)  ## fig 3.C
paretoPlot(model_DOE2_viability)  ## fig 3.B

## plot a contourplot to visualise interactions
contourPlot(model_DOE2_proliferation, 'CD40L_exp', 'CD40L_no') ## fig 3.I
contourPlot(model_DOE2_viability, 'CD40L_exp', 'CD40L_no') ## fig 3.H

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
## fig 3.K
plot_data(DOE2_proliferation, DOE2_viability, BAFF, 'BAFF', 'Number of cells', 'Cell viability (%)')
## fig 3.E
plot_data(DOE2_proliferation, DOE2_viability, IL4, 'IL4', 'Number of cells', 'Cell viability (%)')
## fig 3.J
plot_data(DOE2_proliferation, DOE2_viability, IL21, 'IL21', 'Number of cells', 'Cell viability (%)')
## fig 3.G
plot_data(DOE2_proliferation, DOE2_viability, CD40L_no, 'CD40L_no', 'Number of cells', 'Cell viability (%)')
## fig 3.F
plot_data(DOE2_proliferation, DOE2_viability, CD40L_exp, 'CD40L_exp', 'Number of cells', 'Cell viability (%)')

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 3, height = 2.5, units = "in")

