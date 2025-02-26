#### S4 Fig ####
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(data.table)

# read .fam file with sample names
samps<-read.table("../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.fam")[,1]

# desired order of samples in plot
tritici_order <- read.delim("../ADMIXTURE/tritici_sorted_25042024", header = FALSE)

admx_ord <- function(anc_prop){
  df <- read.delim(anc_prop, sep = " ", header = FALSE)
  rownames(df) <- samps
  df_ord <- df[match(tritici_order$V1,rownames(df)),]
  df_ord_t <- t(as.matrix(df_ord))
  return(df_ord_t)
}

## list of ADMIXTURE output files with least CV error for each value of K between 4-9. 
list_files <- c("../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.k4.r10.Q", 
                "../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.k5.r8.Q", 
                "../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.k6.r3.Q", 
                "../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.k7.r8.Q", 
                "../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.k8.r2.Q", 
                "../ADMIXTURE/tritici_ALL_25kb_0.1_LDp.k9.r10.Q")

tritici_admx <- lapply(list_files, admx_ord)

# dummy barplot to get label segment coordinates for subsequent plots
bp <- barplot(tritici_admx[[1]])
# regions : Northern/Central Europe, Southern Europe, Russia/Central Asia, N. Turkey, S. Turkey, Israel, Egypt, China, Japan, USA/Australia, Argentina
start_samples <- c(bp[1],bp[206],bp[277],bp[313],bp[370],bp[381],bp[409],bp[444],bp[506],bp[516],bp[561])
end_samples <- c(bp[205],bp[276],bp[312],bp[369],bp[380],bp[408],bp[443],bp[505],bp[515],bp[560],bp[568])
start_segments <- start_samples + 0.5
end_segments <- end_samples - 0.5
text_pos<- start_segments+((end_segments - start_segments)/2)

col_pals <- list(c("#999999","#377EB8","#E41A1C","#4DAF4A"),
                 c("#999999","#FFFF33","#E41A1C","#377EB8","#4DAF4A"),
                 c("#FF7F00","#377EB8","#E41A1C","#4DAF4A","#FFFF33","#999999"),
                 c("#FF7F00","#4DAF4A","#FFFF33","#F781BF","#E41A1C","#999999","#377EB8"),
                 c("#377EB8","#F781BF","#FF7F00","#FFFF33","#999999","#E41A1C","#4DAF4A","#984EA3"),
                 c("#E41A1C","#984EA3","#4DAF4A","#87CEEB","#FF7F00","#377EB8","#FFFF33","#F781BF","#999999"))


six_plots <- function(list_of_six, colpal){
  par(mfrow = c(6,1), mar = c(1,1,0.2,0), oma=c(1.5,0,0.3,0))
  for (i in 1:6){
    if (i<6){
      barplot(list_of_six[[i]],col=colpal[[i]], border=colpal[[i]], las = 2, cex.names = 0.3, ylab = "",
              ylim = c(-0.1,1), axes = FALSE,xaxt = "n")
      segments(x0=start_segments,
               x1=end_segments,
               y0=rep(-0.04,11),y1=rep(-0.04,11),lwd = 2, col = "black")
      title(ylab = paste0("K=",i+3), line = -4.5, cex.lab = 1.8)
      } else{
      barplot(list_of_six[[i]],col=colpal[[i]], border=colpal[[i]], las = 2, cex.names = 0.3, ylab = "",
              ylim = c(-0.1,1), axes = FALSE,xaxt = "n")
      segments(x0=start_segments,
               x1=end_segments,
               y0=rep(-0.04,11),y1=rep(-0.04,11),lwd = 2, col = "black")
      title(ylab = paste0("K=",i+3), line = -4.5, cex.lab = 1.8)
      text(x = text_pos, y = -0.17,
           labels = c(paste("North/Central", "Europe", sep = "\n"),paste("South", "Europe", sep = "\n"), 
                      "RUS+", paste("North", "Turkey", sep = "\n"), paste("S", "TUR",sep = "\n"), 
                      "ISR", "EGY", "China", "JPN", paste("USA/","AUS",sep = "\n"), "ARG"),
           xpd = NA, srt = 0, adj = 0.45, cex = 1.5 )
    }
  }
}


six_plots(tritici_admx,col_pals)


#### data to plot S4b can be found in S1 data.

