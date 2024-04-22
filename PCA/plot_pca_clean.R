library(ggplot2)
library(patchwork)
library(Polychrome)
library(RColorBrewer)
library(pals)

setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/pca/")

#### World dataset ####

world_pca <- read.csv("tritici_world.pca.csv")
world_eig <- read.csv("tritici_world.eig.csv")
world_eig$PC <- as.numeric(rownames(world_eig))
world_eig$var_perc <- (world_eig$pca.eig/sum(world_eig$pca.eig))*100
global_regions <- read.csv("../../data/global_regions.csv")
world_pca_reg <- merge(world_pca, global_regions, by.x = "isolate", by.y = "name")


# percentage of variance explained by each principal component
ggplot(data=world_eig[1:20,], aes(PC, var_perc))+geom_bar(stat = "identity")+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs(x="Principal Component", y="Percentage of variance explained",title = "PCA - World dataset")+
  scale_x_continuous(breaks = c(1:20))


# PC1 vs PC2
world_pc12 <- ggplot(world_pca_reg, aes(PC1, PC2, colour=global_region,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(world_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(world_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Region", shape = "Collection") + scale_color_manual(values = as.vector(kelly(12))[2:12]) 

# PC1 vs PC3
world_pc13 <- ggplot(world_pca_reg, aes(PC1, PC3, colour=global_region,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(world_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(world_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Region", shape = "Collection") + scale_color_manual(values = as.vector(kelly(12))[2:12]) 

#PC2 vs PC3
world_pc23 <- ggplot(world_pca_reg, aes(PC2, PC3, colour=global_region,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(world_eig$var_perc[2], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(world_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Region", shape = "Collection") + scale_color_manual(values = (kelly(12))[2:12]) 

world_all <- world_pc12 + world_pc13 + world_pc23 + guide_area() + plot_layout(guides = "collect") & guides(shape=guide_legend(order = 1, ncol = 2 ), colour = guide_legend(order = 2))

ggsave(world_all, filename = "tritici_world_regions_pca_plots.pdf", height = 30, width = 45, unit = "cm")

#### Europe+ dataset ####

eur_pca <- read.csv("tritici_europe+.pca.csv")
eur_eig <- read.csv("tritici_europe+.eig.csv")
eur_eig$PC <- as.numeric(rownames(eur_eig))
eur_eig$var_perc <- (eur_eig$pca.eig/sum(eur_eig$pca.eig))*100

# percentage of variance explained by each principal component
ggplot(data=eur_eig[1:20,], aes(PC, var_perc))+geom_bar(stat = "identity")+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs(x="Principal Component", y="Percentage of variance explained",title = "PCA - Europe+ dataset")+
  scale_x_continuous(breaks = c(1:20))

# modified as.vector(polychrome(30)) to better distinguish between russia and turkey
my_pal <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B",
            "#C4451C", "#B10DA1", "#85660D", "#1C8356", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5")

# PC1 vs PC2
eur_pc12 <- ggplot(eur_pca, aes(PC1, PC2, colour=Country,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Collection") + scale_color_manual(values = my_pal)

# PC1 vc PC3
eur_pc13 <- ggplot(eur_pca, aes(PC1, PC3, colour=Country,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Collection") + scale_color_manual(values = my_pal) 

#PC2 vs PC3
eur_pc23 <- ggplot(eur_pca, aes(PC2, PC3, colour=Country,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Collection") + scale_color_manual(values = my_pal) 

eur_all <- eur_pc12 + eur_pc13 + eur_pc23 + guide_area() + plot_layout(guides = "collect") & guides(shape=guide_legend(order = 1, ncol = 2 ), colour = guide_legend(order = 2))

ggsave(eur_all, filename = "tritici_europe+_pca_plots.pdf", height = 30, width = 45, unit = "cm")

