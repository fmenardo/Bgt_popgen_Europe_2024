# MH plot for zoom in a region isoRelate
library(ggplot2)
library(cowplot)
library(dplyr)

setwd("~/projects/nikos/selection_scans/isoRelate/0_data/")
# Open and read the file list including all the variable names
f <- list.files(full.names = T, pattern = "txt")
f_short1 <- list.files(full.names = F, pattern = "txt")
f_short <- str_sub(f_short1, end = -5)

for (i in 1:length(f)) {
  
  gwscan <- read.table(f[i], header = TRUE)
  
  head(gwscan)

  BC <- -log10(0.05/nrow(gwscan))

  n <- length(gwscan$new_p)

  # Sort the negative log10(p-values) from largest to smallest.
  y <- rev(sort(gwscan$new_p))

  gwscan <- cbind(gwscan, marker = 1:n)

  # Add column "odd.chr" to the table, and find the positions of the chromosomes along the x-axis
  gwscan <- transform(gwscan, odd.chr = (chr %% 2) == 1)
  
  plot.region.pvalues <- function(gwscan, chr_start, chr_end) {
    
    # Convert the positions to the Megabase (Mb) scale.
    gwscan <- transform(gwscan, pos_bp = pos_bp/1e6)
    
    # Create a Manhattan plot.
    return(ggplot(gwscan, aes(x = pos_bp, y = new_p)) +
             geom_point(color = "azure4", size = 4, shape = 20) +
             labs(x = "base-pair position (Mb)", y = "-log10 p-value") +
             geom_hline(yintercept = BC) +
             theme_cowplot(font_size = 10) +
             geom_vline(xintercept = chr_start, color = "purple", linewidth = 0.4, linetype = "dotted") +
             geom_vline(xintercept = chr_end, color = "purple", linewidth = 0.4, linetype = "dotted") +
             theme(axis.line.y = element_line(lineend = "butt"),
                   axis.line = element_blank(),
                   axis.ticks.x = element_line(),
                   plot.title = element_text(hjust = 0.5)))
  }
  
  # zoom in chromosome 8 in the region 5.5 millionth base to 6 millionth base

  AvrPm3e <- subset(gwscan, chr == 4 & pos_bp > 3e6 & pos_bp < 3.25e6)
  MH_plot_zoom <- plot.region.pvalues(AvrPm3e, 3.134197, 3.139303)
  #MH_plot_zoom_cyp51
  
  # create plots directory if it doesn't exist
  
  output_dir <- "../2_output/plots/no_cap/regions/"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  #pdf(paste("2_output/MH_isoRelate_maf0.02_cm4_Avr_Svr.pdf"), width = 12, height = 4)
  #print(MH_Avr_Svr)
  #dev.off()
  saveRDS(MH_plot_zoom, paste0("../2_output/plots/no_cap/regions/MH_",f_short[i],"_AvrPm3e.RDS"))
  png(paste0("../2_output/plots/no_cap/regions/MH_",f_short[i],"_AvrPm3e.png"), width = 12, height = 4, units = "in", res = 150)
  print(MH_plot_zoom)
  dev.off()
  # for png: The units in which height and width are given. Can be px (pixels, the default), in (inches), cm or mm.
}
