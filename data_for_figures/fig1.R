#### Fig 1 ####
library(ape)
library(patchwork)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dendextend)

## Fig 1b dendrogram
my_tree <- read.nexus("../fineStructure/fs_tree.nex")
hc <- as.hclust(my_tree)
dendro <- as.dendrogram(hc)
dend_col <- color_branches(dendro, k=5)
dend_col <- set(dend_col, "labels", NULL)

dendro_plot <- as.ggdend(dend_col) %>% ggplot(linewidth = 0.5)+theme_void()+
  scale_color_manual(values = c("#E41A1C", "#EA9999", "#E5B110", "#377EB8", "#984EA3", "#000000"))+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),axis.title.y = element_blank(),
    axis.title = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank(),axis.text.x=element_blank()
  )+coord_cartesian(expand = FALSE)
  
write.csv(my_tree$tip.label,"tip_order", row.names = FALSE, quote = FALSE) ## save tip order for coancestry matrix and admixture barplot


## Fig 1c fs coancestry matrix

tritici_order <- read.delim("tip_order", header = TRUE)
chunk_mat <- as.matrix(read.table("../fineStructure/Europe_large_linked_hap.chunkcounts.out", row.names=1, header = T, skip = 1))
datamatrix <- chunk_mat[tritici_order$x, tritici_order$x]
datamatrix_log <- log(datamatrix)
coancestry_melt <- melt(as.matrix(datamatrix_log))

coancestry_plot <- ggplot(coancestry_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(aes(color = after_scale(fill))) +
  scale_fill_gradientn(colours = c("#FBD786","#FF3030","#68228B", "#000000"))+
  theme_void() +
  labs(fill = "fineSTRUCTURE chunk count (log)")+
  theme(axis.text.x=element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), legend.position = "none",axis.title.y = element_blank())

        
## Fig 1d admixture barplot

anc_df <- read.delim("../ADMIXTURE/r10_k9_admx_prop_with_names.Q", sep = " ", header = FALSE)
anc_ext_eur <- merge(tritici_order, anc_df, by.x = "x", by.y = "V1")
rownames(anc_ext_eur) <- anc_ext_eur[,1]
anc_ext_eur[,1] <- NULL
anc_df_ord <- anc_ext_eur[match(tritici_order$x, rownames(anc_ext_eur)),]
mypal = c("#440154FF", "#472D7BFF", "#2C728EFF", "#27AD81FF","#21908CFF","#5DC863FF", "#FDE725FF", "#3B528BFF", "#AADC32FF" )
admixture_df <- data.frame(individual=1:nrow(anc_df_ord), anc_df_ord)
admixture_df_melt <- melt(admixture_df, id.vars="individual")

admixture_plot <- ggplot(admixture_df_melt, aes(x=individual, y=value, fill=variable, colour = variable)) +
  geom_bar(stat="identity", width = 0.9) +
  theme_void() +
  #labs(x="Individuals", y="Ancestry Proportion") +
  scale_fill_manual(values = mypal)+scale_color_manual(values = mypal)+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(),legend.position = "none")+coord_cartesian(expand = FALSE)


## Fig 1e pi

MM_fs4_11chr_df <- read.csv("../summary_statistics/fs4_pi_theta_tajimasD_maxmiss0.5.csv")

pi <- ggplot(data = MM_fs4_11chr_df, aes(avg_pi, factor(pop))) + 
  geom_boxplot(aes(fill = factor(pop)), show.legend = FALSE, notch = TRUE, outlier.shape = NA, width = 0.5)+
  theme_classic(base_size = 22)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1), axis.text = element_text(size = 18))+
  scale_colour_manual(values = c("#984EA3","#377EB8","#EA9999","#E41A1C","#E5B110"), aesthetics = "fill")+
  labs(y = "Population", x = "Average pi (per-site, 10kb windows)") +
  coord_cartesian(xlim = c(0,0.004))

        
## Fig 1f PCA

eur_pca <- read.csv("../PCA/tritici_europe+.pca.csv")
eur_eig <- read.csv("../PCA/tritici_europe+.eig.csv")
eur_eig$PC <- as.numeric(rownames(eur_eig))
eur_eig$var_perc <- (eur_eig$pca.eig/sum(eur_eig$pca.eig))*100
new_metadata <- read.csv("../Datasets/S1_Data.csv")
pca_data <- eur_pca[,1:11]
pca_with_new_metadata <- merge(pca_data, new_metadata, by.x = "isolate", by.y = "Sample.ID")
pca_with_new_metadata$collection_3_levels <- ifelse(pca_with_new_metadata$year_of_collection == 2022 | pca_with_new_metadata$year_of_collection == 2023, 
                                                    pca_with_new_metadata$year_of_collection,
                                                    "1990-2018")

pca_with_new_metadata <- pca_with_new_metadata %>% 
  mutate_at(c('Country', 'Region', 'Collection', 'fs_level_4','collection_3_levels'), as.factor)
  
fs4_pal = c("#984EA3", "#377EB8", "#EA9999", "#E41A1C", "#E5B110")
     
# PC1 vs PC2
eur_pc12_pop <- ggplot(pca_with_new_metadata, aes(PC1, PC2, colour=fs_level_4,shape=collection_3_levels)) + geom_point(size=4) + theme_classic(base_size = 20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1), axis.text = element_text(size = 18))+
  labs (x = paste("PC1", " ", "(", signif(eur_eig$var_perc[1], 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(eur_eig$var_perc[2], 3), "%",")",  sep = ""),
        colour = "Population", shape = "Year of sampling") + scale_color_manual(values = fs4_pal)   



