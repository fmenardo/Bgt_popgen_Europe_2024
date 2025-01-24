library(ggplot2)
library(patchwork)
library(Polychrome)
library(RColorBrewer)
library(pals)


#### World dataset ####

world_pca <- read.csv("tritici_world.pca.csv") # output from pca_clean.R
world_eig <- read.csv("tritici_world.eig.csv") # output from pca_clean.R
world_eig$PC <- as.numeric(rownames(world_eig))
world_eig$var_perc <- (world_eig$pca.eig/sum(world_eig$pca.eig))*100
global_regions <- read.csv("global_regions.csv")
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

#### Europe+ dataset ####

eur_pca <- read.csv("tritici_europe+.pca.csv")  # output from pca_clean.R
eur_eig <- read.csv("tritici_europe+.eig.csv") # output from pca_clean.R
eur_eig$PC <- as.numeric(rownames(eur_eig))
eur_eig$var_perc <- (eur_eig$pca.eig/sum(eur_eig$pca.eig))*100

eur_pca$collection_3_levels <- ifelse(eur_pca$year_of_collection == 2022 | eur_pca$year_of_collection == 2023, 
                                            eur_pca$year_of_collection,
                                            "1980-2019")

# percentage of variance explained by each principal component
var_eur <- ggplot(data=eur_eig[1:20,], aes(PC, var_perc))+geom_bar(stat = "identity")+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs(x="Principal Component", y="Percentage of variance explained")+
  scale_x_continuous(breaks = c(1:20))

## modified as.vector(polychrome(30)) to better distinguish between russia and turkey
my_pal <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B",
            "#C4451C", "#B10DA1", "#85660D", "#1C8356", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5")

# PC1 vs PC2
eur_pc12 <- ggplot(eur_pca, aes(PC1, PC2, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal)

# PC1 vc PC3
eur_pc13 <- ggplot(eur_pca, aes(PC1, PC3, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal) 

#PC2 vs PC3
eur_pc23 <- ggplot(eur_pca, aes(PC2, PC3, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal) 

eur_all <- eur_pc12 + eur_pc13 + eur_pc23 + var_eur + plot_layout(guides = "collect")+plot_annotation(tag_levels = 'a')

ggsave(eur_all, filename = "tritici_europe+_pca_plots.pdf", height = 30, width = 45, unit = "cm")


## plots with year of sampling and coloured by population

new_metadata <- read.csv("S1_Data.csv")
pca_data <- eur_pca[,1:11]
pca_with_new_metadata <- merge(pca_data, new_metadata, by.x = "isolate", by.y = "Sample.ID")

pca_with_new_metadata$collection_3_levels <- ifelse(pca_with_new_metadata$year_of_collection == 2022 | pca_with_new_metadata$year_of_collection == 2023, 
                                                    pca_with_new_metadata$year_of_collection,
                                                    "1980-2019")

pca_with_new_metadata <- pca_with_new_metadata %>% 
  mutate_at(c('Country', 'Region', 'Collection', 'fs_level_4','collection_3_levels'), as.factor)

fs4_pal = c("#984EA3", "#377EB8", "#EA9999", "#E41A1C", "#E5B110")

# PC1 vs PC2
eur_pc12_pop <- ggplot(pca_with_new_metadata, aes(PC1, PC2, colour=fs_level_4,shape=collection_3_levels)) + geom_point(size=4) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1), axis.text = element_text(size = 18))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Population", shape = "Year of sampling") + scale_color_manual(values = fs4_pal)
#ggsave(eur_pc12_pop, filename = "tritici_europe+_pca12_pop_4_20.pdf", height = 30, width = 45, unit = "cm")


