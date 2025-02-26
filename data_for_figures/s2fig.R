#### S2 Fig ####
library(ggplot2)
library(patchwork)
library(Polychrome)
library(RColorBrewer)
library(pals)


## World dataset ##

world_pca <- read.csv("../PCA/tritici_world.pca.csv") # output from pca_clean.R
world_eig <- read.csv("../PCA/tritici_world.eig.csv") # output from pca_clean.R
world_eig$PC <- as.numeric(rownames(world_eig))
world_eig$var_perc <- (world_eig$pca.eig/sum(world_eig$pca.eig))*100
global_regions <- read.csv("../PCA/global_regions.csv")
world_pca_reg <- merge(world_pca, global_regions, by.x = "isolate", by.y = "name")

world_pca_reg$collection_3_levels <- ifelse(world_pca_reg$year_of_collection == 2022 | world_pca_reg$year_of_collection == 2023, 
                                            world_pca_reg$year_of_collection,
                                                    "1980-2019")

world_pca_reg$global_region <- as.factor(world_pca_reg$global_region)

# percentage of variance explained by each principal component
var <- ggplot(data=world_eig[1:20,], aes(PC, var_perc))+geom_bar(stat = "identity")+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs(x="Principal Component", y="Percentage of variance explained",title = "PCA - World dataset")+
  scale_x_continuous(breaks = c(1:20))

# PC1 vs PC2
world_pc12 <- ggplot(world_pca_reg %>%
                   arrange(factor(global_region, levels = c("Argentina","United States of America","China",
                                                            "Egypt", "Israel", "Japan", "Northern Europe", "Russia/Central Asia",
                                                            "Southern Europe", "Turkey", "Australia"))), 
                     aes(PC1, PC2, colour=global_region,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(world_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(world_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Region", shape = "Year of sampling") + scale_color_manual(values = as.vector(kelly(12))[2:12]) 

# PC1 vs PC3
world_pc13 <- ggplot(world_pca_reg %>%
                   arrange(factor(global_region, levels = c("Argentina","United States of America","China",
                                                            "Egypt", "Israel", "Japan", "Northern Europe", "Russia/Central Asia",
                                                            "Southern Europe", "Turkey", "Australia"))), 
                     aes(PC1, PC3, colour=global_region,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(world_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(world_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Region", shape = "Year of sampling") + scale_color_manual(values = as.vector(kelly(12))[2:12]) 

#PC2 vs PC3
world_pc23 <- ggplot(world_pca_reg %>%
                   arrange(factor(global_region, levels = c("Argentina","United States of America","China",
                                                            "Egypt", "Israel", "Japan", "Northern Europe", "Russia/Central Asia",
                                                            "Southern Europe", "Turkey", "Australia"))), 
                     aes(PC2, PC3, colour=global_region,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(world_eig$var_perc[2], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(world_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Region", shape = "Year of sampling") + scale_color_manual(values = (kelly(12))[2:12]) 

world_all <- world_pc12 + world_pc13 + world_pc23 + var + plot_layout(guides = "collect")+plot_annotation(tag_levels = 'a')

#ggsave(world_all, filename = "tritici_world_regions_pca_plots.pdf", height = 30, width = 45, unit = "cm")
