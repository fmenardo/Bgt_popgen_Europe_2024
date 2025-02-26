#### S10 Fig ####
library(tidyverse)
library(patchwork)

fst_dxy_df <- read.csv("../summary_statistics/dxy_wc-fst_fs_level4.csv")

d <- ggplot(data = subset(fst_dxy_df, fst_dxy_df$statistic == "avg_dxy"), aes(value, factor(pop_pair) )) + 
  geom_boxplot(fill = "grey",alpha=0.8,show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))+
  labs(y = "Population pair", x = "average dxy (10 kb windows)") + 
  coord_cartesian(xlim=c(0,0.005))+ stat_summary("median")
  
wc <- ggplot(data = subset(fst_dxy_df, fst_dxy_df$statistic == "avg_wc_fst"), aes(value, factor(pop_pair) )) + 
  geom_boxplot(fill = "grey",alpha=0.8, show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))+
  scale_colour_viridis_b()+
  labs(y = "Population pair", x = "average WC-Fst (10 kb windows)") +
  coord_cartesian(xlim=c(-0.25,0.5))
  
d + wc  + plot_layout(axes = "collect") +plot_annotation(tag_levels = 'a') & theme(axis.title = element_text(size = 12)) 

