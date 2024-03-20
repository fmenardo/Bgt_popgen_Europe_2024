library(tidyverse)
library(RColorBrewer)
library(patchwork)
setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/admixture/2022+before2022+2023+ncsu/tritici_ALL_25kb_0.1_LDp_admixture/")

#### cv error ####

log <- read.table("tritici_ALL_25kb_0.1_LDp_admixture_K_r1.log")[,c(3:4)]
log$V3<-gsub("\\(K=", "", log$V3)
log$V3<-gsub("):", "", log$V3)
#interpret K values as numerical
colnames(log)<-c("k","cv")
log$k <- as.numeric(log$k)
log$cv <- as.numeric(log$cv)

p <- ggplot(log, aes(k,cv))+geom_line()+geom_point()+
  ylab("Cross-validation error")+
  xlab("K")+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

p

#### bar plots ####

# read .fam file with sample names
samps<-read.table("tritici_ALL_25kb_0.1_LDp.fam")[,1]

# desired order of samples in plot
tritici_order <- read.delim("../tritici_sorted_updated", header = FALSE)

admx_ord <- function(anc_prop){
  df <- read.delim(anc_prop, sep = " ", header = FALSE)
  rownames(df) <- samps
  df_ord <- df[match(tritici_order$V1,rownames(df)),]
  df_ord_t <- t(as.matrix(df_ord))
  return(df_ord_t)
}

list_files <- c("tritici_ALL_25kb_0.1_LDp.k4.r1.Q", "tritici_ALL_25kb_0.1_LDp.k5.r1.Q", "tritici_ALL_25kb_0.1_LDp.k6.r1.Q",
                "tritici_ALL_25kb_0.1_LDp.k7.r1.Q", "tritici_ALL_25kb_0.1_LDp.k8.r1.Q")

tritici_admx <- lapply(list_files, admx_ord)

# set up plot
par(mfrow=c(5,1))

## K=4 ##
my_pal_4 <- c( "#4DAF4A", "#999999", "#377EB8", "#E41A1C")
bp4 <- barplot(tritici_admx[[1]], col = my_pal_4 , border = my_pal_4 , las = 2, cex.names = 0.5, ylab = "Ancestry", main = "K=4", adj = 0.05)

## K=5 ## 
my_pal_5 <- c( "#E41A1C", "#FFFF33","#4DAF4A","#377EB8","#999999")
bp5 <- barplot(tritici_admx[[2]],col = my_pal_5 , border = my_pal_5 , las = 2, cex.names = 0.5, ylab = "Ancestry", main = "K=5", adj = 0.05)

## K=6 ##
my_pal_6 <- c("#E41A1C", "#999999","#FFFF33","#FF7F00","#4DAF4A","#377EB8")
bp6 <- barplot(tritici_admx[[3]],col = my_pal_6 , border = my_pal_6 , las = 2, cex.names = 0.5, ylab = "Ancestry", main = "K=6", adj = 0.05)

## K=7 ##
my_pal_7 <- c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FFFF33", "#999999" , "#FF7F00")
bp7 <- barplot(tritici_admx[[4]],col = my_pal_7 , border = my_pal_7 , las = 2, cex.names = 0.5, ylab = "Ancestry", main = "K=7", adj = 0.05)

## K=8 ##
my_pal_8 <- c("#4DAF4A","#E41A1C","#FFFF33","#999999", "#F781BF", "#FF7F00" ,"#377EB8","#984EA3")
bp8 <- barplot(tritici_admx[[5]],col = my_pal_8 , border = my_pal_8 , las = 2, cex.names = 0.5, ylab = "Ancestry", main = "K=8", adj = 0.05) 


