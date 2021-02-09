#### Code to create a new theme that looks like base and set as default

# Load libraries 
library(cowplot); library(ggplot2); library(ggthemes)

# Uses theme_base() from the ggthemes library
# Adds theme arguments to make the plot transparent and change axis spacing and sizes 
# Sets to default using theme_set()
theme_set(theme_bw()+theme(panel.grid = element_blank(),
                           plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                           strip.background = element_blank()))
