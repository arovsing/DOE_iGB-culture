## Rovsing et al. 2024
## data and scripts for generating linear model and figures for DOE 3

library(pid)
library(ggplot2)
library(dplyr)

## create variables for each factor with DOE coding in standard order
## here the design is a Box-Behnken design repeated twice
BAFF = c(-1,1,-1,1,0,0,0,0,-1,1,-1,1,0,0,0,0,-1,1,-1,1,0,0,0,0,0,0,0,-1,0,1,-1,1,-1,1,0,0,0,0,-1,1,-1,1,0,0,0,0,-1,1,
         -1,1,0,0,0,0,0,0,0,-1,0,1)
CD40L_no = c(-1,-1,1,1,0,0,0,0,0,0,0,0,-1,1,-1,1,0,0,0,0,-1,1,-1,1,0,0,0,0,1,-1,-1,-1,1,1,0,0,0,0,0,0,0,0,-1,1,-1,1,0,0,
           0,0,-1,1,-1,1,0,0,0,0,1,-1)
IL4 = c(0,0,0,0,-1,1,-1,1,0,0,0,0,-1,-1,1,1,-1,-1,1,1,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,-1,1,-1,1,0,0,0,0,-1,-1,1,1,-1,-1,1,
        1,0,0,0,0,0,0,0,1,0,-1)
IL21 = c(0,0,0,0,-1,-1,1,1,-1,-1,1,1,0,0,0,0,0,0,0,0,-1,-1,1,1,0,0,0,0,-1,1,0,0,0,0,-1,-1,1,1,-1,-1,1,1,0,0,0,0,0,0,0,
         0,-1,-1,1,1,0,0,0,0,-1,1)
CM = c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,
          1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

BAFF = c(BAFF,BAFF)
CD40L_no = c(CD40L_no,CD40L_no)
IL4 = c(IL4,IL4)
IL21 = c(IL21, IL21)
CM = c(CM,CM)
CD40L_exp = rep(c(1,-1),each = 60)

## create variables for each readout in the same DOE standard order
DOE3_proliferation = c(18500,51400,57700,55700,50400,37500,332900,94900,15800,62300,77400,55100,182000,216300,59800,
                       140800,101500,149800,71000,85800,69200,66800,40700,94800,70700,43600,72800,111000,72200,311800,
                       288800,370400,281600,239700,154000,102800,241700,97100,168900,82900,229800,416500,194200,282000,
                       60000,120600,472900,317300,100400,268200,68500,88400,301600,247900,203300,356700,210800,97700,
                       65600,121100,20900,20200,11000,8200,15400,14900,101200,10300,15200,21200,10300,21000,71700,58900,
                       7200,17800,86000,49000,7500,22200,8100,9600,16900,33600,8500,14800,9900,13500,5400,46000,187800,
                       270600,178000,160600,50000,14400,186400,14200,50800,53100,83900,134800,330000,343000,5600,16900,
                       379600,216700,15800,7400,49000,34600,271800,223500,103200,149600,160500,8500,134300,178300)

DOE3_viability = c(37.8,46.2,43.6,49.9,26.7,46.2,42.0,40.0,51.8,43.8,36.5,40.3,38.2,40.4,42.0,42.7,36.1,40.2,38.6,46.2,
                   46.3,42.4,39.3,41.1,36.9,43.3,37.5,36.0,42.8,38.2,32.1,37.2,28.7,28.1,23.6,31.5,30.1,31.0,21.8,21.8,
                   32.2,31.7,24.8,34.4,31.9,29.0,36.4,33.0,39.6,38.4,18.2,19.0,35.2,27.0,29.6,31.3,33.1,31.2,15.0,34.6,
                   55.7,44.5,57.9,61.7,58.6,68.3,71.4,67.8,67.3,78.5,41.9,66.2,81.3,77.8,62.1,50.7,67.2,69.2,72.8,61.0,
                   73.0,66.2,55.6,58.2,50.0,43.9,63.5,48.2,62.1,78.0,88.0,87.1,87.0,87.0,72.7,59.5,84.8,64.3,72.4,63.7,
                   77.0,84.0,86.3,88.3,51.4,58.3,88.2,82.7,62.2,61.2,69.7,61.9,87.6,85.7,83.0,85.8,86.1,59.0,77.2,69.6)

## Europium IgG ng/mL
DOE3_IgG = c(8.150029817,12.94901701,6.345286417,6.536793374,4.067824433,8.327778559,20.88921497,11.75331065,
             6.172748825,8.183318978,7.231370109,12.07789641,17.6413463,13.46397747,17.68184378,11.20115127,
             13.93698087,23.52249753,12.90797084,12.84644484,8.531159206,4.865904195,5.156891573,6.51868097,
             13.41537873,9.816476231,13.37750012,12.91149065,6.484373895,28.98471561,29.01435253,39.51035412,
             22.29086721,46.48168699,12.08951726,13.62351781,22.07412831,13.68827862,8.931253225,15.01328728,
             28.8335732,27.68131151,20.80587305,28.67384152,19.14128554,11.45387089,30.63721944,34.37679891,
             15.93113845,22.31140787,7.496543817,8.235621752,20.30186317,23.64756656,25.55960375,38.47188712,
             28.82592832,17.24932511,12.42822256,38.35318608,3.153618533,2.344815245,1.935216788,0.937790572,
             2.296880697,1.436791957,2.741985753,3.105604755,2.760402605,1.728105501,2.202440983,1.389116604,
             3.927472936,3.089382488,3.170338506,2.861123101,3.126886452,3.030456618,2.127224574,2.285908365,
             2.999215364,1.646764104,2.149559148,2.140039908,3.216661432,1.849670063,1.117894643,1.978852869,
             2.099586006,3.076060179,2.812846194,2.30769988,5.160356986,1.669665595,4.024251648,3.61327811,
             4.623666801,2.453176969,3.423948,4.54161099,3.124011995,1.56600551,3.911960366,6.269599149,2.679560006,
             1.698532482,3.510862901,4.139175646,2.780267178,2.917149453,3.577684933,6.674479676,3.022968103,
             3.030119123,8.032069841,4.270914058,5.521039102,1.85359596,1.981574595,5.925538686)

## Europium IgE mU
DOE3_IgE = c(1881.57,3272.70,2222.78,3235.98,99.46,1027.10,3731.80,4174.02,917.62,1901.32,3478.20,7685.72,2717.22,
             2388.66,4090.68,5164.04,1630.79,3077.46,5614.22,5742.50,1255.65,1217.16,2136.54,2966.70,4104.64,2632.30,
             4884.74,5886.88,1570.49,9580.34,1875.63,3670.74,2219.48,1701.17,81.12,1186.87,1029.20,3678.74,146.88,
             239.41,2212.74,2855.16,569.05,885.09,4698.98,4280.48,319.02,1606.60,2924.84,5845.56,262.92,129.08,
             446.40,1888.97,2193.76,2492.76,1952.73,3730.08,98.50,781.64,0.00,0.00,189.11,0.00,0.00,144.30,0.00,
             0.00,2.71,34.17,0.00,0.00,38.15,0.00,0.00,0.00,16.64,10.73,0.34,35.23,1.18,2.77,37.88,16.15,0.00,0.00,
             8.36,0.00,0.00,3.74,1087.59,51.98,33.85,0.00,0.00,45.62,33.02,0.00,1.02,16.92,0.00,13.60,29.65,25.10,
             0.00,3.00,27.68,159.18,1.74,0.00,0.00,0.00,28.82,14.10,10.10,0.00,130.45,0.00,0.00,0.00)

## create linear models and use backward elimination to eliminate variables,
## which aren't statistically significant (Pr(>|t|) < 0.05)
model = lm(DOE3_proliferation ~ BAFF * CD40L_no * IL4 * IL21 * CD40L_exp * CM)
model = lm(DOE3_proliferation ~ BAFF + CD40L_no + IL4 + IL21 + CD40L_exp + CM)
## final model
model_DOE3_proliferation = lm(DOE3_proliferation ~ BAFF + IL4 + IL21 + CM + CD40L_exp + I(IL21^2) + BAFF*IL4*CM +
                                  IL21*CM)

## final model
model_DOE3_viability = lm(DOE3_viability ~ IL4 + IL21 + CM + CD40L_exp + CM*CD40L_exp + IL4*CD40L_exp + CM*IL21)
    
    
## final model
model_DOE3_IgG = lm(DOE3_IgG ~ BAFF + IL4 + IL21 + CM + CD40L_exp + I(IL21^2) + BAFF*CD40L_exp + IL21*CD40L_exp)


## final model
model_DOE3_IgE = lm(DOE3_IgE ~ BAFF + IL4 + IL21 + CM + CD40L_exp + I(IL21^2) + BAFF*CD40L_exp + IL4*CD40L_exp + 
                        IL21*CD40L_exp + CM*CD40L_exp)

## print the summary of the model fit
summary(model_DOE3_proliferation)
summary(model_DOE3_viability)
summary(model_DOE3_IgG)
summary(model_DOE3_IgE)

## plot a pareto plot for the final model
paretoPlot(model_DOE3_proliferation)  ## fig 4.C
paretoPlot(model_DOE3_viability)  ## fig 4.B
paretoPlot(model_DOE3_IgG)  ## fig 4.I
paretoPlot(model_DOE3_IgE)  ## fig 4.J

## plot a contourplot to visualise interactions
contourPlot(model_DOE3_proliferation, 'IL21', 'IL4') ## fig 4.H

## plot model diagnostics
par(mfrow=c(2,4)) ## plots all four plots in one window
plot(model)

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
## fig S5.A
plot_data_3levels(DOE3_proliferation, DOE3_viability, BAFF, 'BAFF', 'Number of cells', 'Cell viability (%)')
## fig 4.E
plot_data_3levels(DOE3_proliferation, DOE3_viability, IL4, 'IL4', 'Number of cells', 'Cell viability (%)')
## fig 4.G
plot_data_3levels(DOE3_proliferation, DOE3_viability, IL21, 'IL21', 'Number of cells', 'Cell viability (%)')
## fig S5.B
plot_data_3levels(DOE3_proliferation, DOE3_viability, CD40L_no, 'CD40L_no', 'Number of cells', 'Cell viability (%)')
## fig 4.D
plot_data_2levels(DOE3_proliferation, DOE3_viability, CD40L_exp, 'CD40L_exp', 'Number of cells', 'Cell viability (%)')
## fig 4.F
plot_data_2levels(DOE3_proliferation, DOE3_viability, CM, 'Change medium', 'Number of cells', 'Cell viability (%)')

## fig S5.G
plot_data_3levels(DOE3_IgG, DOE3_IgE, BAFF, 'BAFF', 'IgG (ng/mL)', 'IgE (mU)')
## fig S5.D
plot_data_3levels(DOE3_IgG, DOE3_IgE, IL4, 'IL4', 'IgG (ng/mL)', 'IgE (mU)')
## fig S5.C
plot_data_3levels(DOE3_IgG, DOE3_IgE, IL21, 'IL21', 'IgG (ng/mL)', 'IgE (mU)')
## fig S5.F
plot_data_3levels(DOE3_IgG, DOE3_IgE, CD40L_no, 'CD40L_no', 'IgG (ng/mL)', 'IgE (mU)')
## fig 4.K
plot_data_2levels(DOE3_IgG, DOE3_IgE, CD40L_exp, 'CD40L_exp', 'IgG (ng/mL)', 'IgE (mU)')
## fig S5.E
plot_data_2levels(DOE3_IgG, DOE3_IgE, CM, 'Change medium', 'IgG (ng/mL)', 'IgE (mU)')

## save plot
#ggsave('this.pdf', plot = last_plot(), width = 3, height = 2.5, units = "in")

## find maximum
int = seq(-1,1,by=0.1)
BAFF = CD40L_exp = IL4 = IL21 = CM = int
design = expand.grid(BAFF = BAFF, CD40L_exp = CD40L_exp, IL4 = IL4, IL21 = IL21, CM = CM)
df.find_max = predict(model_DOE3_proliferation, design)
max(df.find_max)
design[which.max(df.find_max), ]
