library(tidyverse)
library(rehh)
library(vcfR)

setwd("~/projects/nikos/selection_scans/xp-EHH/")

# xp-ehh

## Sabeti et al. 2007:
## We developed, evaluated and applied a new test,
## Cross Population Extended Haplotype Homozogysity (XP-EHH),
## to detect selective sweeps in which the selected allele has approached
## or achieved fixation in one population but remains polymorphic in the human population as a whole.
population_scans_ME <- readRDS("2_output/population_tritici_recent_extended_europe_fs_level4_ME.vcf.gz_scans.RDS")
population_scans_TUR <- readRDS("2_output/population_tritici_recent_extended_europe_fs_level4_TUR.vcf.gz_scans.RDS")
population_scans_S_EUR1 <- readRDS("2_output/population_tritici_recent_extended_europe_fs_level4_S_EUR1.vcf.gz_scans.RDS")
population_scans_S_EUR2 <- readRDS("2_output/population_tritici_recent_extended_europe_fs_level4_S_EUR2.vcf.gz_scans.RDS")
population_scans_N_EUR <- readRDS("2_output/population_tritici_recent_extended_europe_fs_level4_N_EUR.vcf.gz_scans.RDS")

population_scans_ME[12] <- "ME"
population_scans_TUR[12] <- "TUR"
population_scans_S_EUR1[12] <- "S_EUR1"
population_scans_S_EUR2[12] <- "S_EUR2"
population_scans_N_EUR[12] <- "N_EUR"

population_scans <- list(population_scans_ME,
                              population_scans_TUR,
                              population_scans_S_EUR1,
                              population_scans_S_EUR2,
                              population_scans_N_EUR)

for (i in 1:4) {
  if (i==1){
    secondary = 2
  } else if (i==2) {
    secondary = 3
  } else if (i==3) {
    secondary = 4
  } else if (i==4) {
    secondary = 5
  }
  for (j in secondary:5) {
    xpehh_chr1 <- ies2xpehh(as.data.frame(population_scans[[i]][1]), 
                           as.data.frame(population_scans[[j]][1]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr2 <- ies2xpehh(as.data.frame(population_scans[[i]][2]), 
                           as.data.frame(population_scans[[j]][2]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr3 <- ies2xpehh(as.data.frame(population_scans[[i]][3]), 
                           as.data.frame(population_scans[[j]][3]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr4 <- ies2xpehh(as.data.frame(population_scans[[i]][4]), 
                           as.data.frame(population_scans[[j]][4]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr5 <- ies2xpehh(as.data.frame(population_scans[[i]][5]), 
                           as.data.frame(population_scans[[j]][5]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr6 <- ies2xpehh(as.data.frame(population_scans[[i]][6]), 
                           as.data.frame(population_scans[[j]][6]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr7 <- ies2xpehh(as.data.frame(population_scans[[i]][7]), 
                           as.data.frame(population_scans[[j]][7]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr8 <- ies2xpehh(as.data.frame(population_scans[[i]][8]), 
                           as.data.frame(population_scans[[j]][8]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr9 <- ies2xpehh(as.data.frame(population_scans[[i]][9]), 
                           as.data.frame(population_scans[[j]][9]),
                           as.data.frame(population_scans[[i]][12]),
                           as.data.frame(population_scans[[j]][12]))
    xpehh_chr10 <- ies2xpehh(as.data.frame(population_scans[[i]][10]), 
                            as.data.frame(population_scans[[j]][10]),
                            as.data.frame(population_scans[[i]][12]),
                            as.data.frame(population_scans[[j]][12]))
    xpehh_chr11 <- ies2xpehh(as.data.frame(population_scans[[i]][11]), 
                            as.data.frame(population_scans[[j]][11]),
                            as.data.frame(population_scans[[i]][12]),
                            as.data.frame(population_scans[[j]][12]))
    xpehh_all <- rbind(xpehh_chr1, xpehh_chr2, xpehh_chr3, xpehh_chr4, xpehh_chr5, xpehh_chr6, xpehh_chr7, xpehh_chr8, xpehh_chr9, xpehh_chr10, xpehh_chr11)
    write.csv(xpehh_all, paste0("2_output/xpehh_",population_scans[[i]][12],"_",population_scans[[j]][12],".csv"))
  }
}

# xpehh MH

#
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)

setwd("~/projects/nikos/selection_scans/xp-EHH/2_output/")

# Open and read the file list including all the variable names
f <- list.files(full.names = T, pattern = "csv")
f_short1 <- list.files(full.names = F, pattern = "csv")
f_short <- str_sub(f_short1, end = -5)

for (i in 1:length(f)) {
  
  gwscan1 <- read.csv(f[i], header = TRUE)
  
  # exclude rows with NaN
  gwscan <- gwscan1[complete.cases(gwscan1), ]
  head(gwscan)
  
  names(gwscan)[4] <- "XPEHH"
  
  BC <- -log10(0.05/nrow(gwscan))
  iHS_top1percent <- quantile(gwscan$XPEHH, probs = c(.99), na.rm = T)
  iHS_top01percent <- quantile(gwscan$XPEHH, probs = c(.999), na.rm = T)
  iHS_bot01percent <- quantile(gwscan$XPEHH, probs = c(.001), na.rm = T)
  n <- length(gwscan$XPEHH)
  
  # Sort the negative log10(p-values) from largest to smallest.
  #y <- rev(sort(gwscan$LOGPVALUE))
  
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
  
  # Add column "odd.CHR" to the table, and find the positions of the chromosomes along the x-axis
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
  MH_Avr_Svr <- ggplot(gwscan, aes(x = marker, y = XPEHH, color = odd.CHR)) +
    geom_point(size = 1.5, shape = 20) +
    geom_vline(xintercept = gene_of_int_Avr_Svr_marker, color = "red3", linewidth = 0.4, linetype = "dotted") +
    geom_vline(xintercept = gene_of_int_fungicide, color = "dodgerblue3", linewidth = 0.4, linetype = "dotted") +
    geom_text(aes(x=gene_of_int_fungicide[1], label="Btub", y=max(gwscan$XPEHH) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[2], label="cyp51", y=max(gwscan$XPEHH) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[3], label="sdhB", y=max(gwscan$XPEHH) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[4], label="sdhC", y=max(gwscan$XPEHH) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_fungicide[5], label="sdhD", y=max(gwscan$XPEHH) - 1), colour="dodgerblue3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[1], label="AvrPm3d3", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[1], label="AvrPm3d3", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[3], label="AvrPm2", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[4], label="AvrPm3a2/f2", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[5], label="SvrPm3a1/f1", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[6], label="AvrPm3b2/c2", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[7], label="AvrPm1a.1", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[8], label="AvrPm1a.2", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[9], label="AvrPm17", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[11], label="AvrPm8", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    geom_text(aes(x=gene_of_int_Avr_Svr_marker[12], label="AvrPm3e", y=max(gwscan$XPEHH) - 1), colour="red3", angle=60, size = 2.5) +
    scale_x_continuous(breaks = x.CHR, labels = 1:11) +
    scale_y_continuous(limits = c(min(gwscan$XPEHH) - 1, max(gwscan$XPEHH) + 1), expand = c(0,0)) +
    scale_color_manual(values = c("black", "azure4"), guide = "none") +
    labs(x = "", y = "xp-EHH") +
    geom_hline(yintercept = iHS_bot01percent, linetype = "dotted") +
    geom_hline(yintercept = iHS_top01percent, linetype = "dotted") +
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



