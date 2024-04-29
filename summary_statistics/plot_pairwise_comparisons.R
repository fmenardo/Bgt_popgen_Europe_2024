### plot distributions of pairwise comparisons (dxy and fst) as calculated by pixy
library(tidyverse)
setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/summary_statistics/pixy/")

#### fs4 ####
# list files to be read
# fst wc
all_fst_files <- list.files("output_fs_level4_10kb", pattern = "fst.txt", full.names = TRUE)
fst_unwanted <- c("output_fs_level4_10kb/pixy_tritici_ext_eur_recent_fs_level4_Bgt_MAT_1_1_3_fst.txt",
                  "output_fs_level4_10kb/pixy_tritici_ext_eur_recent_fs_level4_LR026995.1_Un_fst.txt",
                  "output_fs_level4_10kb/pixy_tritici_ext_eur_recent_fs_level4_MT880591.1_fst.txt")
fst_11chr <- setdiff(all_fst_files, fst_unwanted)


# dxy
all_dxy_files <- list.files("output_fs_level4_10kb", pattern = "dxy.txt", full.names = TRUE)
dxy_unwanted <- c("output_fs_level4_10kb/pixy_tritici_ext_eur_recent_fs_level4_Bgt_MAT_1_1_3_dxy.txt",
                  "output_fs_level4_10kb/pixy_tritici_ext_eur_recent_fs_level4_LR026995.1_Un_dxy.txt",
                  "output_fs_level4_10kb/pixy_tritici_ext_eur_recent_fs_level4_MT880591.1_dxy.txt")
dxy_11chr <- setdiff(all_dxy_files, dxy_unwanted)

read_and_combine <- function(fst_list, dxy_list){
  df <- lapply(c(fst_list, dxy_list), read_tsv)
  pixy_df <- lapply(df, pivot_longer, -c(pop1, pop2, window_pos_1, window_pos_2, chromosome), names_to="statistic", values_to="value") 
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  }

fst_dxy_df <- read_and_combine(fst_11chr, dxy_11chr)
fst_dxy_df$pop_pair <- paste(fst_dxy_df$pop1, fst_dxy_df$pop2, sep = "/")

#### plot dxy fs4 ####
ggplot(data = subset(fst_dxy_df, fst_dxy_df$statistic == "avg_dxy"), aes(factor(pop_pair), value)) + 
  geom_violin(aes(fill = factor(pop_pair)), draw_quantiles = c(0.25,0.5,0.75), show.legend = FALSE)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_viridis_b()+
  labs(x = "Population pair", y = "average dxy (10 kb windows)")

d <- ggplot(data = subset(fst_dxy_df, fst_dxy_df$statistic == "avg_dxy"), aes(value, factor(pop_pair) )) + 
  geom_boxplot(aes(fill = factor(pop_pair)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_viridis_b()+
  labs(y = "Population pair", x = "average dxy (10 kb windows)") + 
  coord_cartesian(xlim=c(0,0.005))

#### plot Fst fs4 ####
ggplot(data = subset(fst_dxy_df, fst_dxy_df$statistic == "avg_wc_fst"), aes(factor(pop_pair), value)) + 
  geom_violin(aes(fill = factor(pop_pair)), draw_quantiles = c(0.25,0.5,0.75), show.legend = FALSE)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_viridis_b()+
  labs(x = "Population pair", y = "average Fst (10 kb windows)")

wc <- ggplot(data = subset(fst_dxy_df, fst_dxy_df$statistic == "avg_wc_fst"), aes(value, factor(pop_pair) )) + 
  geom_boxplot(aes(fill = factor(pop_pair)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_viridis_b()+
  labs(y = "Population pair", x = "average WC Fst (10 kb windows)") +
  coord_cartesian(xlim=c(-0.25,0.5))

d + wc  + plot_layout(axes = "collect") & theme(axis.title = element_text(size = 12))

