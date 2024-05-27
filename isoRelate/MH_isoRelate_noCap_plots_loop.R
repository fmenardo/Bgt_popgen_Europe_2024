library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

setwd("~/projects/nikos/selection_scans/isoRelate/0_data/")

# Open and read the file list including all the variable names
f <- list.files(full.names = T, pattern = "txt")
f_short1 <- list.files(full.names = F, pattern = "txt")
f_short <- str_sub(f_short1, end = -5)

for (i in 1:length(f)) {
  
  gwscan <- read.table(f[i], header = TRUE)
  
  head(gwscan)
  
  BC <- -log10(0.05/nrow(gwscan))
  FDR <- function(pvals, FDR){
    pvalss <- sort(pvals, decreasing = F)
    m = length(pvalss)
    cutoffs <- ((1:m)/m)*FDR
    logicvec <- pvalss <= cutoffs
    postrue <- which(logicvec)
    k <- max(c(postrue,0))
    cutoff <- (((0:m)/m)*FDR)[k + 1]
    return(cutoff)
  }
  
  fdrlog10 <- FDR(gwscan$new_p, 0.05)
  
  n <- length(gwscan$new_p)
  
  # Sort the negative log10(p-values) from largest to smallest.
  y <- rev(sort(gwscan$new_p))
  
  gwscan <- cbind(gwscan, marker = 1:n)
  
  # Add column "odd.chr" to the table, and find the positions of the chromosomes along the x-axis
  gwscan <- transform(gwscan, odd.chr = (chr %% 2) == 1)
  x.chr  <- tapply(gwscan$marker, gwscan$chr, mean)
  
  # Create the genome-wide scan (no highlights)
  MH <- ggplot(gwscan, aes(x = marker, y = new_p, color = odd.chr)) +
    geom_point(size = 2, shape = 20) +
    scale_x_continuous(breaks = x.chr, labels = 1:11) +
    scale_y_continuous(limits = c(0, max(y) + 1), expand = c(0,0)) +
    scale_color_manual(values = c("black", "azure4"), guide = "none") +
    labs(x = "",y = "-log10 p-value") +
    geom_hline(yintercept = BC) +
    geom_hline(yintercept = fdrlog10, linetype = "dotted") +
    theme_cowplot() +
    theme(axis.line.y = element_line(lineend = "butt"),
          axis.line = element_blank(),
          axis.ticks.x = element_line(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, size = 10))
  
  # Create a data-set to highlight. Using my gene middle point for Avrs and Svr. ()
  chr <- c("9", "9", "7", "6", "4", "5", "6", "8", "1", "1", "11", "4")
  pos_bp <- c("347003", "369134", "935062", "2135001", "10109438", "18863147", "4863234", "10725805", "4365208", "4370275", "2660889", "3136750")
  gene_of_int_Avr_Svr_mean <- data.frame(chr, pos_bp)
  
  gwscan_chr1 <- gwscan %>% filter(chr==1)
  gwscan_chr2 <- gwscan %>% filter(chr==2)
  gwscan_chr3 <- gwscan %>% filter(chr==3)
  gwscan_chr4 <- gwscan %>% filter(chr==4)
  gwscan_chr5 <- gwscan %>% filter(chr==5)
  gwscan_chr6 <- gwscan %>% filter(chr==6)
  gwscan_chr7 <- gwscan %>% filter(chr==7)
  gwscan_chr8 <- gwscan %>% filter(chr==8)
  gwscan_chr9 <- gwscan %>% filter(chr==9)
  gwscan_chr10 <- gwscan %>% filter(chr==10)
  gwscan_chr11 <- gwscan %>% filter(chr==11)
  
  # Avr/Svr
  gene_of_int_Avr_Svr_marker <- 0
  gene_of_int_Avr_Svr_marker[1] <- gwscan_chr9[which(abs(gwscan_chr9$pos_bp-347003)==min(abs(gwscan_chr9$pos_bp-347003))), "marker"]
  gene_of_int_Avr_Svr_marker[2] <- gwscan_chr9[which(abs(gwscan_chr9$pos_bp-369134)==min(abs(gwscan_chr9$pos_bp-369134))), "marker"]
  gene_of_int_Avr_Svr_marker[3] <- gwscan_chr7[which(abs(gwscan_chr7$pos_bp-935062)==min(abs(gwscan_chr7$pos_bp-935062))), "marker"]
  gene_of_int_Avr_Svr_marker[4] <- gwscan_chr6[which(abs(gwscan_chr6$pos_bp-2135001)==min(abs(gwscan_chr6$pos_bp-2135001))), "marker"]
  gene_of_int_Avr_Svr_marker[5] <- gwscan_chr4[which(abs(gwscan_chr4$pos_bp-10109438)==min(abs(gwscan_chr4$pos_bp-10109438))), "marker"]
  gene_of_int_Avr_Svr_marker[6] <- gwscan_chr5[which(abs(gwscan_chr5$pos_bp-18863147)==min(abs(gwscan_chr5$pos_bp-18863147))), "marker"]
  gene_of_int_Avr_Svr_marker[7] <- gwscan_chr6[which(abs(gwscan_chr6$pos_bp-4863234)==min(abs(gwscan_chr6$pos_bp-4863234))), "marker"]
  gene_of_int_Avr_Svr_marker[8] <- gwscan_chr8[which(abs(gwscan_chr8$pos_bp-10725805)==min(abs(gwscan_chr8$pos_bp-10725805))), "marker"]
  gene_of_int_Avr_Svr_marker[9] <- gwscan_chr1[which(abs(gwscan_chr1$pos_bp-4365208)==min(abs(gwscan_chr1$pos_bp-4365208))), "marker"]
  gene_of_int_Avr_Svr_marker[10] <- gwscan_chr1[which(abs(gwscan_chr1$pos_bp-4370275)==min(abs(gwscan_chr1$pos_bp-4370275))), "marker"]
  gene_of_int_Avr_Svr_marker[11] <- gwscan_chr11[which(abs(gwscan_chr11$pos_bp-2660889)==min(abs(gwscan_chr11$pos_bp-2660889))), "marker"]
  gene_of_int_Avr_Svr_marker[12] <- gwscan_chr4[which(abs(gwscan_chr4$pos_bp-3136750)==min(abs(gwscan_chr4$pos_bp-3136750))), "marker"]
  
  # fungicide resistance
  # get means
  mean(c(7766015,7767656))
  mean(c(5728014,5729706))
  mean(c(5629949,5630812))
  mean(c(3466484,3467217))
  mean(c(13068412,13069018))
  
  
  gene_of_int_fungicide <- 0
  gene_of_int_fungicide[1] <- gwscan_chr10[which(abs(gwscan_chr10$pos_bp-7766836)==min(abs(gwscan_chr10$pos_bp-7766836))), "marker"]
  gene_of_int_fungicide[2] <- gwscan_chr8[which(abs(gwscan_chr8$pos_bp-5728860)==min(abs(gwscan_chr8$pos_bp-5728860))), "marker"]
  gene_of_int_fungicide[3] <- gwscan_chr9[which(abs(gwscan_chr9$pos_bp-5630380)==min(abs(gwscan_chr9$pos_bp-5630380))), "marker"]
  gene_of_int_fungicide[4] <- gwscan_chr3[which(abs(gwscan_chr3$pos_bp-3466850)==min(abs(gwscan_chr3$pos_bp-3466850))), "marker"]
  gene_of_int_fungicide[5] <- gwscan_chr2[which(abs(gwscan_chr2$pos_bp-13068715)==min(abs(gwscan_chr2$pos_bp-13068715))), "marker"]
  
  # Create the genome-wide scan with ALL Avr and Svr
  MH_Avr_Svr <- ggplot(gwscan, aes(x = marker, y = new_p, color = odd.chr)) +
    geom_point(size = 1.5, shape = 20) +
    geom_vline(xintercept = gene_of_int_Avr_Svr_marker, color = "red", linewidth = 0.2) +
    geom_vline(xintercept = gene_of_int_fungicide, color = "dodgerblue3", linewidth = 0.4, linetype = "dotted") +
    geom_text(aes(x=gene_of_int_fungicide[1], label="Btub", y=max(y) - 20), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[2], label="cyp51", y=max(y) - 20), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[3], label="sdhB", y=max(y) - 20), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[4], label="sdhC", y=max(y) - 20), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[5], label="sdhD", y=max(y) - 20), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[1], label="AvrPm3d3", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[3], label="AvrPm2", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[4], label="AvrPm3a2/f2", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[5], label="SvrPm3a1/f1", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[6], label="AvrPm3b2/c2", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[7], label="AvrPm1a.1", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[8], label="AvrPm1a.2", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[9], label="AvrPm17", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[11], label="AvrPm8", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[12], label="AvrPm3e", y=max(y) - 20), colour="red3", angle=60, size = 2.5) +
    scale_x_continuous(breaks = x.chr, labels = 1:11) +
    scale_y_continuous(limits = c(0, max(y) + 1), expand = c(0,0)) +
    scale_color_manual(values = c("black", "azure4"), guide = "none") +
    labs(x = "", y = "-log10 p-value") +
    geom_hline(yintercept = BC) +
    #geom_hline(yintercept = fdrlog10, linetype = "dotted") +
    theme_cowplot() +
    theme(axis.line.y = element_line(lineend = "butt"),
          axis.line = element_blank(),
          axis.ticks.x = element_line(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, size = 10))
  
  
  # create plots directory if it doesn't exist
  
  output_dir <- "../2_output/plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  #pdf(paste("2_output/MH_isoRelate_maf0.02_cm4_Avr_Svr.pdf"), width = 12, height = 4)
  #print(MH_Avr_Svr)
  #dev.off()
  saveRDS(MH_Avr_Svr, paste0("../2_output/plots/MH_",f_short[i],".RDS"))
  png(paste0("../2_output/plots/MH_",f_short[i],".png"), width = 12, height = 4, units = "in", res = 150)
  print(MH_Avr_Svr)
  dev.off()
# for png: The units in which height and width are given. Can be px (pixels, the default), in (inches), cm or mm.
}
