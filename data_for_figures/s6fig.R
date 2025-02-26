#### S6 Fig ####
## code same as Fig1, just changed colour palettes
library(ape)
library(patchwork)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dendextend)

## S6a dendrogram
my_tree <- read.nexus("../fineStructure/fs_tree.nex")
hc <- as.hclust(my_tree)
dendro <- as.dendrogram(hc)

dendro_p <- as.ggdend(dendro) %>% ggplot(linewidth = 0.1)+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),axis.title.y = element_blank(),
    axis.title = element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank(),axis.text.x=element_blank()
  )+coord_cartesian(expand = FALSE)

## S6b fs coancestry matrix
tritici_order <- my_tree$tip.label  ## order of samples
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

## S6c admixture barplot
anc_df <- read.delim("../ADMIXTURE/r10_k9_admx_prop_with_names.Q", sep = " ", header = FALSE)
anc_ext_eur <- merge(tritici_order, anc_df, by.x = "x", by.y = "V1")
rownames(anc_ext_eur) <- anc_ext_eur[,1]
anc_ext_eur[,1] <- NULL
anc_df_ord <- anc_ext_eur[match(tritici_order$x, rownames(anc_ext_eur)),]
mypal = c("#E41A1C","#984EA3","#4DAF4A","#87CEEB","#FF7F00","#377EB8","#FFFF33","#F781BF","#999999")
admixture_df <- data.frame(individual=1:nrow(anc_df_ord), anc_df_ord)
admixture_df_melt <- melt(admixture_df, id.vars="individual")

admixture_plot <- ggplot(admixture_df_melt, aes(x=individual, y=value, fill=variable, colour = variable)) +
  geom_bar(stat="identity", width = 0.9) +
  theme_void() +
  #labs(x="Individuals", y="Ancestry Proportion") +
  scale_fill_manual(values = mypal)+scale_color_manual(values = mypal)+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(),legend.position = "none")+coord_cartesian(expand = FALSE)
        
## combine the plots
comb_pl <- plot_grid(dendro_p, coancestry_plot, admixture_plot, ncol = 1, align = 'v', rel_heights = c(2,6,1), scale = 0.97)+theme_map()

