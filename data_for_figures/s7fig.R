#### S7 Fig ####
library(ape)
library(patchwork)
library(tidyverse)
library(reshape2)
library(ggplot2)

my_tree <- read.nexus("../fineStructure/fs_tree.nex")
tritici_order <- my_tree$tip.label
chunk_mat <- as.matrix(read.table("../fineStructure/Europe_large_linked_hap.chunkcounts.out", row.names=1, header = T, skip = 1))
datamatrix <- chunk_mat[tritici_order$x, tritici_order$x]
datamatrix_log <- log(datamatrix)

## subset to include only samples belonging to S_EUR2 population
es_subs <- datamatrix_log[134:202,134:202]

coancestry_melt_es <- melt(as.matrix(es_subs))

coancestry_plot_es <- ggplot(coancestry_melt_es, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(aes(color = after_scale(fill))) +
  scale_fill_gradientn(colours = c("#FBD786","#FF3030","#68228B", "#000000"))+
  labs(fill = paste("fineSTRUCTURE","chunk count (log)", sep = "\n"))+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), 
       axis.title.x = element_blank(), legend.position = "right",axis.title.y = element_blank())

