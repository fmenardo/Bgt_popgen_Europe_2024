library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

setwd("~/projects/nikos/selection_scans/iHS/2_output/")

# Open and read the file list including all the variable names
f <- list.files(full.names = T, pattern = "csv")
f_short1 <- list.files(full.names = F, pattern = "csv")
f_short <- str_sub(f_short1, end = -5)

for (i in 1:length(f)) {
  
  gwscan1 <- read.csv(f[i], header = TRUE)
  
  # exclude rows with NaN
  gwscan <- gwscan1[complete.cases(gwscan1), ]
  head(gwscan)
  
  BC <- -log10(0.05/nrow(gwscan))
  iHS_top1percent <- quantile(gwscan$LOGPVALUE, probs = c(.99), na.rm = T)
  iHS_top01percent <- quantile(gwscan$LOGPVALUE, probs = c(.999), na.rm = T)
  n <- length(gwscan$LOGPVALUE)
  
  # Sort the negative log10(p-values) from largest to smallest.
  y <- rev(sort(gwscan$LOGPVALUE))
  
  gwscan <- cbind(gwscan, marker = 1:n)
  
  ## Replace CHRomosome names with numbers
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
  
  # Add column "odd.CHR" to the table, and find the positions of the CHRomosomes along the x-axis
  gwscan <- transform(gwscan, odd.CHR = (CHR %% 2) == 1)
  x.CHR  <- tapply(gwscan$marker, gwscan$CHR, mean)
  
  # Create a data-set to highlight. Using my gene middle point for Avrs and Svr. ()
  CHR <- c("9", "9", "7", "6", "4", "5", "6", "8", "1", "1", "11", "4")
  pos_bp <- c("347003", "369134", "935062", "2135001", "10109438", "18863147", "4863234", "10725805", "4365208", "4370275", "2660889", "3136750")
  gene_of_int_Avr_Svr_mean <- data.frame(CHR, pos_bp)
  
  gwscan_CHR1 <- gwscan %>% filter(CHR==1)
  gwscan_CHR2 <- gwscan %>% filter(CHR==2)
  gwscan_CHR3 <- gwscan %>% filter(CHR==3)
  gwscan_CHR4 <- gwscan %>% filter(CHR==4)
  gwscan_CHR5 <- gwscan %>% filter(CHR==5)
  gwscan_CHR6 <- gwscan %>% filter(CHR==6)
  gwscan_CHR7 <- gwscan %>% filter(CHR==7)
  gwscan_CHR8 <- gwscan %>% filter(CHR==8)
  gwscan_CHR9 <- gwscan %>% filter(CHR==9)
  gwscan_CHR10 <- gwscan %>% filter(CHR==10)
  gwscan_CHR11 <- gwscan %>% filter(CHR==11)
  
  gene_of_int_Avr_Svr_marker <- 0
  gene_of_int_Avr_Svr_marker[1] <- gwscan_CHR9[which(abs(gwscan_CHR9$POSITION-347003)==min(abs(gwscan_CHR9$POSITION-347003))), "marker"]
  gene_of_int_Avr_Svr_marker[2] <- gwscan_CHR9[which(abs(gwscan_CHR9$POSITION-369134)==min(abs(gwscan_CHR9$POSITION-369134))), "marker"]
  gene_of_int_Avr_Svr_marker[3] <- gwscan_CHR7[which(abs(gwscan_CHR7$POSITION-935062)==min(abs(gwscan_CHR7$POSITION-935062))), "marker"]
  gene_of_int_Avr_Svr_marker[4] <- gwscan_CHR6[which(abs(gwscan_CHR6$POSITION-2135001)==min(abs(gwscan_CHR6$POSITION-2135001))), "marker"]
  gene_of_int_Avr_Svr_marker[5] <- gwscan_CHR4[which(abs(gwscan_CHR4$POSITION-10109438)==min(abs(gwscan_CHR4$POSITION-10109438))), "marker"]
  gene_of_int_Avr_Svr_marker[6] <- gwscan_CHR5[which(abs(gwscan_CHR5$POSITION-18863147)==min(abs(gwscan_CHR5$POSITION-18863147))), "marker"]
  gene_of_int_Avr_Svr_marker[7] <- gwscan_CHR6[which(abs(gwscan_CHR6$POSITION-4863234)==min(abs(gwscan_CHR6$POSITION-4863234))), "marker"]
  gene_of_int_Avr_Svr_marker[8] <- gwscan_CHR8[which(abs(gwscan_CHR8$POSITION-10725805)==min(abs(gwscan_CHR8$POSITION-10725805))), "marker"]
  gene_of_int_Avr_Svr_marker[9] <- gwscan_CHR1[which(abs(gwscan_CHR1$POSITION-4365208)==min(abs(gwscan_CHR1$POSITION-4365208))), "marker"]
  gene_of_int_Avr_Svr_marker[10] <- gwscan_CHR1[which(abs(gwscan_CHR1$POSITION-4370275)==min(abs(gwscan_CHR1$POSITION-4370275))), "marker"]
  gene_of_int_Avr_Svr_marker[11] <- gwscan_CHR11[which(abs(gwscan_CHR11$POSITION-2660889)==min(abs(gwscan_CHR11$POSITION-2660889))), "marker"]
  gene_of_int_Avr_Svr_marker[12] <- gwscan_CHR4[which(abs(gwscan_CHR4$POSITION-3136750)==min(abs(gwscan_CHR4$POSITION-3136750))), "marker"]
  
  # fungicide resistance

  gene_of_int_fungicide <- 0
  gene_of_int_fungicide[1] <- gwscan_CHR10[which(abs(gwscan_CHR10$POSITION-7766836)==min(abs(gwscan_CHR10$POSITION-7766836))), "marker"]
  gene_of_int_fungicide[2] <- gwscan_CHR8[which(abs(gwscan_CHR8$POSITION-5728860)==min(abs(gwscan_CHR8$POSITION-5728860))), "marker"]
  gene_of_int_fungicide[3] <- gwscan_CHR9[which(abs(gwscan_CHR9$POSITION-5630380)==min(abs(gwscan_CHR9$POSITION-5630380))), "marker"]
  gene_of_int_fungicide[4] <- gwscan_CHR3[which(abs(gwscan_CHR3$POSITION-3466850)==min(abs(gwscan_CHR3$POSITION-3466850))), "marker"]
  gene_of_int_fungicide[5] <- gwscan_CHR2[which(abs(gwscan_CHR2$POSITION-13068715)==min(abs(gwscan_CHR2$POSITION-13068715))), "marker"]
  
  # Create the genome-wide scan with ALL Avr and Svr
  MH_Avr_Svr <- ggplot(gwscan, aes(x = marker, y = LOGPVALUE, color = odd.CHR)) +
    geom_point(size = 1.5, shape = 20) +
    geom_vline(xintercept = gene_of_int_Avr_Svr_marker, color = "red3", linewidth = 0.4, linetype = "dotted") +
    geom_vline(xintercept = gene_of_int_fungicide, color = "dodgerblue3", linewidth = 0.4, linetype = "dotted") +
    geom_text(aes(x=gene_of_int_fungicide[1], label="Btub", y=max(y) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[2], label="cyp51", y=max(y) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[3], label="sdhB", y=max(y) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[4], label="sdhC", y=max(y) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[5], label="sdhD", y=max(y) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[1], label="AvrPm3d3", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[1], label="AvrPm3d3", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[3], label="AvrPm2", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[4], label="AvrPm3a2/f2", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[5], label="SvrPm3a1/f1", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[6], label="AvrPm3b2/c2", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[7], label="AvrPm1a.1", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[8], label="AvrPm1a.2", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[9], label="AvrPm17", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[11], label="AvrPm8", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[12], label="AvrPm3e", y=max(y) - 1), colour="red3", angle=60, size = 2.5) +
    scale_x_continuous(breaks = x.CHR, labels = 1:11) +
    scale_y_continuous(limits = c(0, max(y) + 1), expand = c(0,0)) +
    scale_color_manual(values = c("black", "azure4"), guide = "none") +
    labs(x = "", y = "-log10 p-value") +
    geom_hline(yintercept = iHS_top01percent) +
    geom_hline(yintercept = iHS_top1percent, linetype = "dotted") +
    theme_cowplot() +
    ggtitle(f_short[i]) +
    theme(axis.line.y = element_line(lineend = "butt"),
          axis.line = element_blank(),
          axis.ticks.x = element_line(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, size = 10))
  
  # create plots directory if it doesn't exist
  
  output_dir <- "./plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  #pdf(paste("2_output/MH_isoRelate_maf0.02_cm4_Avr_Svr.pdf"), width = 12, height = 4)
  #print(MH_Avr_Svr)
  #dev.off()
  saveRDS(MH_Avr_Svr, paste0("./plots/MH_",f_short[i],".RDS"))
  png(paste0("./plots/MH_",f_short[i],".png"), width = 12, height = 4, units = "in", res = 150)
  print(MH_Avr_Svr)
  dev.off()
# for png: The units in which height and width are given. Can be px (pixels, the default), in (inches), cm or mm.
}