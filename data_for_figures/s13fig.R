#### S13 Fig ####
set.seed(123)
library(patchwork)
library(tidyverse)

gendist_all <- as.matrix(read.csv("../distance_matrix/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))

## Europe+_2022_2023
samples_2022_2023 <- readLines("../Datasets/tritici_2022_2023.args")
gendist_2022_2023 <- gendist_all[samples_2022_2023,samples_2022_2023]
xy <- t(combn(colnames(gendist_2022_2023), 2))
gg_2022_2023 <- data.frame(xy, dist=gendist_2022_2023[xy])

## N_EUR_2022_2023
samples_2022_2023_neur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_N_EUR.args")
gendist_2022_2023_neur <- gendist_all[samples_2022_2023_neur,samples_2022_2023_neur]
xy_neur <- t(combn(colnames(gendist_2022_2023_neur), 2))
gg_2022_2023_neur <- data.frame(xy_neur, dist=gendist_2022_2023_neur[xy_neur])

## S_EUR2_2022_2023
samples_2022_2023_seur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_S_EUR+.args")
gendist_2022_2023_seur <- gendist_all[samples_2022_2023_seur,samples_2022_2023_seur]
xy_seur <- t(combn(colnames(gendist_2022_2023_seur), 2))
gg_2022_2023_seur <- data.frame(xy_seur, dist=gendist_2022_2023_seur[xy_seur])

my_cols <- c("#000000", "#377EB8", "#E41A1C")
all_dists <- bind_rows(gg_2022_2023, gg_2022_2023_neur, gg_2022_2023_seur, .id = "ID")
a <- ggplot(all_dists, aes(dist, colour = ID, fill = ID))+geom_density(alpha = 0.4 ) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))+
  labs(y = "Density", x = "Pairwise genetic distance", colour = "Dataset", fill = "Dataset")+
  scale_color_manual(values = my_cols, labels = c("Europe+_2022_2023", "N_EUR_2022_2023", "S_EUR_2022_2023"))+
  scale_fill_manual(values = my_cols, labels = c("Europe+_2022_2023", "N_EUR_2022_2023", "S_EUR_2022_2023"))

