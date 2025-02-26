#### S18 Fig ####
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(ggpubr)

# Open and read the file list including all the variable names
f <- list.files(path = "../isoRelate/", full.names = T, pattern = "iR_table.txt")

names <- gsub("\\./BgtE\\+r_","",as.character(f))
names <- gsub("_2cM_iR_table.txt","",names)
names <- gsub("\\+","2",names)

###### merge tables for all populations into one with all snps
merged_data <- read.table(f[1], header = TRUE)
BC <- -log10(0.05/nrow(merged_data))
merged_data <- subset(merged_data, select = -c(pos_M,snp_id,pop,subpop,iR,log10_pvalue))
colnames(merged_data)<-c("chr","pos_bp",paste0(names[1],"_p_value"))

for (i in 2:length(f)) {
  print(i)
  gwscan <- read.table(f[i], header = TRUE)
  
  
  BC <- c(BC,-log10(0.05/nrow(gwscan)))
  
  gwscan <- subset(gwscan, select = -c(pos_M,snp_id,pop,subpop,iR,log10_pvalue))
  colnames(gwscan)<-c("chr","pos_bp",paste0(names[i],"_p_value"))
  merged_data <- merge(merged_data,gwscan,by="pos_bp",all=TRUE)
  for (r in 1:nrow(merged_data)){
    if (is.na(merged_data$chr.x[r])){merged_data$chr.x[r]=merged_data$chr.y[r]}
    if (is.na(merged_data$chr.y[r])){merged_data$chr.y[r]=merged_data$chr.x[r]}
  }
  merged_data <- subset(merged_data, select = -c(chr.y))
  names(merged_data)[names(merged_data) == "chr.x"] <- "chr"
  merged_data <- merged_data[order(merged_data$chr, merged_data$pos_bp), ]
}

n <- nrow(merged_data)
merged_data <- cbind(merged_data, marker = 1:n)

# Add column "odd.chr" to the table, and find the positions of the chromosomes along the x-axis
merged_data <- transform(merged_data, odd.chr = (chr %% 2) == 1)

#### find positions of avr and fungidies targets
# Create a data-set to highlight. Using my gene middle point for Avrs and Svr. ()
# avrpm3d,avrpm3d,avrpm2m,avrpm3a2/f2,Svrpm3a1/f1,Avrpm3b2/c2,avrpm1.1.avrpm1.2,avrpm17,avrpm17,avrpm8


gwscan_chr1 <- merged_data %>% filter(chr==1)
gwscan_chr2 <- merged_data %>% filter(chr==2)
gwscan_chr3 <- merged_data %>% filter(chr==3)
gwscan_chr4 <- merged_data %>% filter(chr==4)
gwscan_chr5 <- merged_data %>% filter(chr==5)
gwscan_chr6 <- merged_data %>% filter(chr==6)
gwscan_chr7 <- merged_data %>% filter(chr==7)
gwscan_chr8 <- merged_data %>% filter(chr==8)
gwscan_chr9 <- merged_data %>% filter(chr==9)
gwscan_chr10 <- merged_data %>% filter(chr==10)
gwscan_chr11 <- merged_data %>% filter(chr==11)

# Avr/Svr
gene_of_int_Avr_Svr_marker <- 0
gene_of_int_Avr_Svr_marker[1] <- gwscan_chr9[which(abs(gwscan_chr9$pos_bp-347003)==min(abs(gwscan_chr9$pos_bp-347003))), "marker"]
gene_of_int_Avr_Svr_marker[2] <- gwscan_chr9[which(abs(gwscan_chr9$pos_bp-369134)==min(abs(gwscan_chr9$pos_bp-369134))), "marker"]
gene_of_int_Avr_Svr_marker[3] <- gwscan_chr7[which(abs(gwscan_chr7$pos_bp-935062)==min(abs(gwscan_chr7$pos_bp-935062))), "marker"]
gene_of_int_Avr_Svr_marker[4] <- gwscan_chr6[which(abs(gwscan_chr6$pos_bp-2135001)==min(abs(gwscan_chr6$pos_bp-2135001))), "marker"]
gene_of_int_Avr_Svr_marker[5] <- gwscan_chr4[which(abs(gwscan_chr4$pos_bp-10109438)==min(abs(gwscan_chr4$pos_bp-10109438))), "marker"]
gene_of_int_Avr_Svr_marker[6] <- gwscan_chr1[which(abs(gwscan_chr1$pos_bp-4365208)==min(abs(gwscan_chr1$pos_bp-4365208))), "marker"]
gene_of_int_Avr_Svr_marker[7] <- gwscan_chr1[which(abs(gwscan_chr1$pos_bp-4370275)==min(abs(gwscan_chr1$pos_bp-4370275))), "marker"]

# fungicide resistance
# get means
mean(c(5728014,5729706))
gene_of_int_fungicide <- 0
gene_of_int_fungicide[1] <- gwscan_chr8[which(abs(gwscan_chr8$pos_bp-5728860)==min(abs(gwscan_chr8$pos_bp-5728860))), "marker"]

#### make plot 5 populations
col_names <- c("N_EUR_p_value","S_EUR2_p_value","TUR_p_value","S_EUR1_p_value","ME_p_value")#colnames(merged_data)[3:7]


merged_data1 <- merged_data


colours=c("#377EB8","#E41A1C","#E5B110","#EA9999","#984EA3")

list_plots <- list()

for (i in 1:length(col_names)){
  print(i)
  col <- col_names[i]
  
  merged_data1[[col]][merged_data[[col]] < 0.4] <- NA  
  merged_data1[[col]][merged_data[[col]] > 50] <- 50  
  y <- rev(sort(merged_data1[[col]]))
  
  
  
  list_plots[[col]] <- ggplot(merged_data1, aes(x = marker, y = .data[[col]], color = odd.chr))+  #coord_cartesian(ylim = c(0, 80)) +
    geom_point(size = 0.5, shape = 20) +
    scale_x_continuous(breaks = x.chr, labels = 1:11) +
    scale_y_continuous(limits = c(0, 52), expand = c(0,0)) +
    scale_color_manual(values = c("black", "azure4"), guide = "none") +
    labs(x = "",y = "-log10 p-value") +
    geom_hline(yintercept = BC[1]) +
    geom_vline(xintercept = gene_of_int_Avr_Svr_marker, color = "black",
               linewidth = 0.4, linetype = "dotted") +
    geom_vline(xintercept = gene_of_int_fungicide, color = "black",
               linewidth = 0.4, linetype = "dotted")+
    theme_cowplot() +
    theme(axis.line.y = element_line(lineend = "butt"),
          axis.line = element_blank(),
          axis.ticks.x = element_line(),
          axis.text = element_text(size = 5),
          axis.title = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 5))
  
  
  list_plots[[col]] <- list_plots[[col]] +  
    annotate("text", x = gene_of_int_fungicide[1]-15000, y = 50, label = "cyp51",size = 2)+
    annotate("text", x = gene_of_int_Avr_Svr_marker[1]+5000, y = 50, label = "AvrPm3d3",size = 2)+
    annotate("text", x = gene_of_int_Avr_Svr_marker[3], y = 50, label = "AvrPm2",size = 2)+
    annotate("text", x = gene_of_int_Avr_Svr_marker[4], y = 50, label = "AvrPm3a2/f2",size = 2)+
    annotate("text", x = gene_of_int_Avr_Svr_marker[5], y = 50, label = "SvrPm3a1/f1",size = 2)+
    annotate("text", x = gene_of_int_Avr_Svr_marker[6], y = 50, label = "AvrPm17",size = 2)+
    annotate("text", x = x.chr[11], y = 50, label = gsub("_p_value","",col),size = 3,color=colours[i],fontface="bold")
  
}

ggarrange(list_plots[[1]],list_plots[[2]],list_plots[[3]],list_plots[[4]],list_plots[[5]],nrow=5,align="v")

ggsave("Manhattan_plots.pdf", width = 20, height = 30, unit="cm")
