# MH plot for zoom in a region iHS
library(ggplot2)
library(cowplot)
library(dplyr)

setwd("~/projects/nikos/selection_scans/iHS/2_output/")

# Open and read the file list including all the variable names
f <- list.files(full.names = T, pattern = "csv")
f_short1 <- list.files(full.names = F, pattern = "csv")
f_short <- str_sub(f_short1, end = -5)

for (i in 1:length(f)) {
  gwscan1 <- read.csv(f[i], header = TRUE)
  gwscan <- gwscan1[complete.cases(gwscan1), ]
  head(gwscan)
  
  BC <- -log10(0.05/nrow(gwscan))
  iHS_top1percent <- quantile(gwscan$LOGPVALUE, probs = c(.99), na.rm = T)
  iHS_top01percent <- quantile(gwscan$LOGPVALUE, probs = c(.999), na.rm = T)
  n <- length(gwscan$LOGPVALUE)
  
  # Sort the negative log10(p-values) from largest to smallest.
  y <- rev(sort(gwscan$LOGPVALUE))
  
  gwscan <- cbind(gwscan, marker = 1:n)
  
  ## Replace chromosome names with numbers
  unique(gwscan$CHR)
  gwscan$CHR[gwscan$CHR == "LR026984.1_chr1"] <- 1
  gwscan$CHR[gwscan$CHR == "LR026985.1_chr2"] <- 2
  gwscan$CHR[gwscan$CHR == "LR026986.1_chr3"] <- 3
  gwscan$CHR[gwscan$CHR == "LR026987.1_chr4"] <- 4
  gwscan$CHR[gwscan$CHR == "LR026988.1_chr5"] <- 5
  gwscan$CHR[gwscan$CHR == "LR026989.1_chr6"] <- 6
  gwscan$CHR[gwscan$CHR == "LR026990.1_chr7"] <- 7
  gwscan$CHR[gwscan$CHR == "LR026991.1_chr8"] <- 8
  gwscan$CHR[gwscan$CHR == "LR026992.1_chr9"] <- 9
  gwscan$CHR[gwscan$CHR == "LR026993.1_chr10"] <- 10
  gwscan$CHR[gwscan$CHR == "LR026994.1_chr11"] <- 11
  
  gwscan$CHR <- as.numeric(gwscan$CHR)
  
  # Add column "odd.CHR" to the table, and find the positions of the chromosomes along the x-axis
  gwscan <- transform(gwscan, odd.CHR = (CHR %% 2) == 1)
  x.CHR  <- tapply(gwscan$marker, gwscan$CHR, mean)
  
  plot.region.pvalues <- function(gwscan, chr_start, chr_end) {
    
    # Convert the positions to the Megabase (Mb) scale.
    gwscan <- transform(gwscan,POSITION = POSITION/1e6)
    
    # Create a Manhattan plot.
    return(ggplot(gwscan,aes(x = POSITION,y = LOGPVALUE)) +
             geom_point(color = "azure4",size = 4, shape = 20) +
             labs(x = "base-pair position (Mb)",y = "-log10 p-value") +
             geom_hline(yintercept = iHS_top1percent, linetype = "dotted") +
             geom_hline(yintercept = iHS_top01percent) +
             geom_vline(xintercept = chr_start, color = "purple", linewidth = 0.4, linetype = "dotted") +
             geom_vline(xintercept = chr_end, color = "purple", linewidth = 0.4, linetype = "dotted") +
             theme_cowplot(font_size = 10) +
             theme(axis.line.y = element_line(lineend = "butt"),
                   axis.line = element_blank(),
                   axis.ticks.x = element_line(),
                   plot.title = element_text(hjust = 0.5)))
  }
  
  # zoom in chromosome 8 in the region 5.5 millionth base to 6 millionth base
  
  AvrPm3e <- subset(gwscan, CHR == 4 & POSITION > 3e6 & POSITION < 3.25e6)
  MH_plot_zoom <- plot.region.pvalues(AvrPm3e, 3.134197, 3.139303)
  
  #MH_plot_zoom

  # create plots directory if it doesn't exist
  
  output_dir <- "./plots/regions/"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  #pdf(paste("2_output/MH_isoRelate_maf0.02_cm4_Avr_Svr.pdf"), width = 12, height = 4)
  #print(MH_Avr_Svr)
  #dev.off()
  saveRDS(MH_plot_zoom, paste0("./plots/regions/MH_",f_short[i],"_AvrPm3e.RDS"))
  png(paste0("./plots/regions/MH_",f_short[i],"_AvrPm3e.png"), width = 12, height = 4, units = "in", res = 150)
  print(MH_plot_zoom)
  dev.off()
  # for png: The units in which height and width are given. Can be px (pixels, the default), in (inches), cm or mm.
}
