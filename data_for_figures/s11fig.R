#### S11 Fig ####
library(tidyverse)
library(patchwork)

MM_fs4_11chr_df <- read.csv("../summary_statistics/fs4_pi_theta_tajimasD_maxmiss0.5.csv")

d<- ggplot(data = MM_fs4_11chr_df, aes(tajimasD, factor(pop))) + 
  geom_boxplot(aes(fill = factor(pop)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))+
  scale_colour_manual(values = c("#984EA3","#377EB8","#EA9999","#E41A1C","#E5B110"), aesthetics = "fill")+
  labs(y = "Population", x = "Tajima's D (10kb windows)")+
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed")+
  coord_cartesian(xlim = c(-3,3))
  
t <- ggplot(data = MM_fs4_11chr_df, aes(wattersons_theta, factor(pop))) + 
  geom_boxplot(aes(fill = factor(pop)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))+
  scale_colour_manual(values = c("#984EA3","#377EB8","#EA9999","#E41A1C","#E5B110"), aesthetics = "fill")+
  labs(y = "Population", x = "Watterson's theta (10kb windows)") +
  coord_cartesian(xlim = c(0,0.0075))
  
t+d+plot_layout(axes = "collect") +plot_annotation(tag_levels = 'a') & theme(axis.title = element_text(size = 12))
