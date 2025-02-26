#### S5 Fig ####
library(ggplot2)
library(patchwork)
library(Polychrome)
library(RColorBrewer)
library(pals)
library(tidyverse)

## stringent GATK filters
eur_pca <- read.csv("../PCA/rev_filt_tritici_europe+.pca.csv")
eur_eig <- read.csv("../PCA/rev_filt_tritici_europe+.eig.csv")
eur_eig$PC <- as.numeric(rownames(eur_eig))
eur_eig$var_perc <- (eur_eig$pca.eig/sum(eur_eig$pca.eig))*100
eur_pca$collection_3_levels <- ifelse(eur_pca$Collection == 2022 | eur_pca$Collection == 2023, 
                                      eur_pca$Collection,
                                      "1980-2019")

## modified as.vector(polychrome(30)) to better distinguish between russia and turkey
my_pal <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B",
            "#C4451C", "#B10DA1", "#85660D", "#1C8356", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5")

# percentage of variance explained by each principal component
var_eur <- ggplot(data=eur_eig[1:20,], aes(PC, var_perc))+geom_bar(stat = "identity")+
  theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs(x="Principal Component", y="% variance explained")+
  scale_x_continuous(breaks = c(1:20))
  
# PC1 vs PC2
eur_pc12 <- ggplot(eur_pca, aes(PC1, PC2, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal)
# PC1 vc PC3
eur_pc13 <- ggplot(eur_pca, aes(PC1, PC3, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal) 

#PC2 vs PC3
eur_pc23 <- ggplot(eur_pca, aes(PC2, PC3, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal) 
        

## usual GATK filters
old_eur_pca <- read.csv("../PCA/tritici_europe+.pca.csv")
old_eur_eig <- read.csv("../PCA/tritici_europe+.eig.csv")
old_eur_eig$PC <- as.numeric(rownames(old_eur_eig))
old_eur_eig$var_perc <- (old_eur_eig$pca.eig/sum(old_eur_eig$pca.eig))*100
old_eur_pca$collection_3_levels <- ifelse(old_eur_pca$Collection == 2022 | old_eur_pca$Collection == 2023, 
                                      old_eur_pca$Collection,
                                      "1980-2019")

old_var_eur <- ggplot(data=old_eur_eig[1:20,], aes(PC, var_perc))+geom_bar(stat = "identity")+
  theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs(x="Principal Component", y="% variance explained")+
  scale_x_continuous(breaks = c(1:20))

# PC1 vs PC2
old_eur_pc12 <- ggplot(old_eur_pca, aes(PC1, PC2, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(old_eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(old_eur_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal)
# PC1 vc PC3
old_eur_pc13 <- ggplot(old_eur_pca, aes(PC1, PC3, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(old_eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(old_eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal) 

#PC2 vs PC3
old_eur_pc23 <- ggplot(old_eur_pca, aes(PC2, PC3, colour=Country,shape=collection_3_levels)) + geom_point(size=3) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(old_eur_eig$var_perc[2], 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(old_eur_eig$var_perc[3], 3), "%",")",  sep = ""),
        colour = "Country", shape = "Year of sampling") + scale_color_manual(values = my_pal) 
a <- (old_eur_pc12 + old_eur_pc13 + old_eur_pc23 + old_var_eur) / (eur_pc12 + eur_pc13 + eur_pc23 + var_eur) + plot_layout(guides = "collect")+plot_annotation(tag_levels = 'a')

ggsave(a, filename = "S5Fig.pdf", height = 45, width = 40, unit = "cm")
