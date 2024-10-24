### calculate tajima's d using pi from pixy output and S from vcftools --SNPdensity
### calculations based on https://ocw.mit.edu/courses/hst-508-quantitative-genomics-fall-2005/0900020a633e85338e9510495c2e01a6_tajimad1.pdf 
### Ref: Tajima 1989  doi: 10.1093/genetics/123.3.585
### Watterson's theta = S/a1 (https://en.wikipedia.org/wiki/Watterson_estimator) 
### for theta, we use an estimate of S normalised by the actual number of sites in each window to make the distributions across windows comparable

library(tidyverse)
library(patchwork)
setwd("/summary_statistics/tajimas_d/")

## set up functions to calculate Tajima's D 

a1 <- function(n){i <- 1:(n-1) ; return(sum(1/i))}
a2 <- function(n){i <- 1:(n-1) ; return(sum(1/i^2))}
b1 <- function(n){return((n+1)/(3*(n-1)))}
b2 <- function(n){return(2*(n^2+n+3)/(9*n*(n-1)))}
c1 <- function(n){return(b1(n) - 1/a1(n))}
c2 <- function(n){return(b2(n) - ((n+2)/(a1(n)*n)) + (a2(n)/a1(n)^2))}
e1 <- function(n){return(c1(n)/a1(n))}
e2 <- function(n){return(c2(n)/(a1(n)^2+a2(n)))}

tajimas_d <- function(pi,S,n){return((pi-S/a1(n))/sqrt(e1(n)*S + e2(n)*S*(S-1)))}

#### FS LEVEL 4 ####
## read vcftools and pixy output to extract S and pi

MM_S_ME <- read_tsv("tritici_fs4_ME_maxmiss0.5.snpden")
MM_S_ME$BIN_START_PLUS_ONE <- MM_S_ME$BIN_START + 1

MM_S_N_EUR <- read_tsv("tritici_fs4_N_EUR_maxmiss0.5.snpden")
MM_S_N_EUR$BIN_START_PLUS_ONE <- MM_S_N_EUR$BIN_START + 1

MM_S_S_EUR <- read_tsv("tritici_fs4_S_EUR+_maxmiss0.5.snpden")
MM_S_S_EUR$BIN_START_PLUS_ONE <- MM_S_S_EUR$BIN_START + 1

MM_S_S_EUR1 <- read_tsv("tritici_fs4_S_EUR1_maxmiss0.5.snpden")
MM_S_S_EUR1$BIN_START_PLUS_ONE <- MM_S_S_EUR1$BIN_START + 1

MM_S_TUR <- read_tsv("tritici_fs4_TUR_maxmiss0.5.snpden")
MM_S_TUR$BIN_START_PLUS_ONE <- MM_S_TUR$BIN_START + 1

pi_all_files <- list.files("../pixy/output_fs_level4_10kb_per_pop", pattern = "pi.txt", full.names = TRUE)
pi_Un_files <- c(list.files("../pixy/output_fs_level4_10kb_per_pop", pattern = "Un_pi.txt", full.names = TRUE),
                 list.files("../pixy/output_fs_level4_10kb_per_pop", pattern = "MT880591.1_pi.txt", full.names = TRUE),
                 list.files("../pixy/output_fs_level4_10kb_per_pop", pattern = "MAT_1_1_3_pi.txt", full.names = TRUE))
pi_files <- setdiff(pi_all_files,pi_Un_files)  # keep only 11 chromosomes

read_pi <- function(pi_files){
  pi_chr <- lapply(pi_files, read_tsv)
  pi_WG <- do.call(rbind, pi_chr)
  return (pi_WG)
}
MM_pi_WG <- read_pi(pi_files)
MM_pi_WG$pi_in_window <- MM_pi_WG$avg_pi*MM_pi_WG$no_sites

## fs4_ME (n = 42) ##
MM_fs4_ME_df <- merge(subset(MM_pi_WG, MM_pi_WG$pop == "ME"), MM_S_ME, by.x = c("chromosome","window_pos_1"), by.y = c("CHROM", "BIN_START_PLUS_ONE"))
MM_fs4_ME_df$tajimasD <- tajimas_d(MM_fs4_ME_df$pi_in_window, MM_fs4_ME_df$SNP_COUNT, 42)
MM_fs4_ME_df$wattersons_theta <- MM_fs4_ME_df$SNP_COUNT / (MM_fs4_ME_df$no_sites*a1(42))

## fs4_N_EUR (n = 188) ##
MM_fs4_N_EUR_df <- merge(subset(MM_pi_WG, MM_pi_WG$pop == "N_EUR"), MM_S_N_EUR, by.x = c("chromosome","window_pos_1"), by.y = c("CHROM", "BIN_START_PLUS_ONE"))
MM_fs4_N_EUR_df$tajimasD <- tajimas_d(MM_fs4_N_EUR_df$pi_in_window, MM_fs4_N_EUR_df$SNP_COUNT, 188)
MM_fs4_N_EUR_df$wattersons_theta <- MM_fs4_N_EUR_df$SNP_COUNT / (MM_fs4_N_EUR_df$no_sites*a1(188))

## fs4_S_EUR+ (n = 69) ##
MM_fs4_S_EUR_df <- merge(subset(MM_pi_WG, MM_pi_WG$pop == "S_EUR+"), MM_S_S_EUR, by.x = c("chromosome","window_pos_1"), by.y = c("CHROM", "BIN_START_PLUS_ONE"))
MM_fs4_S_EUR_df$tajimasD <- tajimas_d(MM_fs4_S_EUR_df$pi_in_window, MM_fs4_S_EUR_df$SNP_COUNT, 69)
MM_fs4_S_EUR_df$wattersons_theta <- MM_fs4_S_EUR_df$SNP_COUNT / (MM_fs4_S_EUR_df$no_sites*a1(69))

## fs4_S_EUR1 (n = 12) ##
MM_fs4_S_EUR1_df <- merge(subset(MM_pi_WG, MM_pi_WG$pop == "S_EUR1"), MM_S_S_EUR1, by.x = c("chromosome","window_pos_1"), by.y = c("CHROM", "BIN_START_PLUS_ONE"))
MM_fs4_S_EUR1_df$tajimasD <- tajimas_d(MM_fs4_S_EUR1_df$pi_in_window, MM_fs4_S_EUR1_df$SNP_COUNT, 12)
MM_fs4_S_EUR1_df$wattersons_theta <- MM_fs4_S_EUR1_df$SNP_COUNT / (MM_fs4_S_EUR1_df$no_sites*a1(12))

## fs4_TUR (n = 57) ##
MM_fs4_TUR_df <- merge(subset(MM_pi_WG, MM_pi_WG$pop == "TUR"), MM_S_TUR, by.x = c("chromosome","window_pos_1"), by.y = c("CHROM", "BIN_START_PLUS_ONE"))
MM_fs4_TUR_df$tajimasD <- tajimas_d(MM_fs4_TUR_df$pi_in_window, MM_fs4_TUR_df$SNP_COUNT, 57)
MM_fs4_TUR_df$wattersons_theta <- MM_fs4_TUR_df$SNP_COUNT / (MM_fs4_TUR_df$no_sites*a1(57))


#### PLOTS FS4 ####

MM_fs4_11chr_df <- rbind(MM_fs4_ME_df, MM_fs4_N_EUR_df, MM_fs4_S_EUR_df, MM_fs4_S_EUR1_df, MM_fs4_TUR_df)
write_csv(MM_fs4_11chr_df, "fs4_pi_theta_tajimasD_maxmiss0.5.csv")

# tajimas D
ggplot(data = MM_fs4_11chr_df, aes(tajimasD,factor(pop))) + 
  geom_violin(aes(fill = factor(pop)), draw_quantiles = c(0.25,0.5,0.75), show.legend = FALSE)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_manual(values = c("#FF7F00","#377EB8","#E41A1C","#EA9999","#FFFF33"), aesthetics = "fill")+
  labs(y = "Population", x = "Tajima's D (10kb windows)")

d<- ggplot(data = MM_fs4_11chr_df, aes(tajimasD, factor(pop))) + 
  geom_boxplot(aes(fill = factor(pop)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_manual(values = c("#FF7F00","#377EB8","#E41A1C","#EA9999","#FFFF33"), aesthetics = "fill")+
  labs(y = "Population", x = "Tajima's D (10kb windows)")+
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed")+
  coord_cartesian(xlim = c(-3,3))

# watterson's theta
ggplot(data = MM_fs4_11chr_df, aes(wattersons_theta, factor(pop),)) + 
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), aes(fill = factor(pop)), show.legend = FALSE)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_manual(values = c("#FF7F00","#377EB8","#E41A1C","#EA9999","#FFFF33"), aesthetics = "fill")+
  labs(y = "Population", x = "Watterson's theta (10kb windows)")

t <- ggplot(data = MM_fs4_11chr_df, aes(wattersons_theta, factor(pop))) + 
  geom_boxplot(aes(fill = factor(pop)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_manual(values = c("#FF7F00","#377EB8","#E41A1C","#EA9999","#FFFF33"), aesthetics = "fill")+
  labs(y = "Population", x = "Watterson's theta (10kb windows, normalised)") +
  coord_cartesian(xlim = c(0,0.0075))

# pi
ggplot(data = MM_fs4_11chr_df, aes(avg_pi, factor(pop))) + 
  geom_violin(aes(fill = factor(pop)), draw_quantiles = c(0.25,0.5,0.75), show.legend = FALSE)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_manual(values = c("#FF7F00","#377EB8","#E41A1C","#EA9999","#FFFF33"), aesthetics = "fill")+
  labs(y = "Population", x = "Average pi (10kb windows)")

p <- ggplot(data = MM_fs4_11chr_df, aes(avg_pi, factor(pop))) + 
  geom_boxplot(aes(fill = factor(pop)), show.legend = FALSE, notch = TRUE, outlier.shape = NA)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_colour_manual(values = c("#FF7F00","#377EB8","#E41A1C","#EA9999","#FFFF33"), aesthetics = "fill")+
  labs(y = "Population", x = "Average pi (per-site, 10kb windows)") +
  coord_cartesian(xlim = c(0,0.005))

p + t + d + plot_layout(axes = "collect") & theme(axis.title = element_text(size = 12))


