#####  1.3.24: something is off with adegenet?? (https://github.com/thibautjombart/adegenet/issues/307)
#### TL;DR: use as.vector(dist objects) for functions like lm, kde2d. used to work without it, not anymore.

library(adegenet)
library(MASS)
library(fields)
library(tidyverse)

setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/isolation_by_distance/")
gendist_all <- as.matrix(read.csv("../dist_mat/gw_dist_mat_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
geo_loc_all <- read.csv("../../../vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK7_19032024.csv")

#### Tritici Europe+_recent ####

samples_eur <- readLines("../../../vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args")
sample_df_eur = as.data.frame(samples_eur)

# extract samples from distance matrix
tritici_2022_2023_gendist_eur <- gendist_all[samples_eur,samples_eur]
gendist_matr_eur <- as.dist(tritici_2022_2023_gendist_eur)

# extract coordinates from metadata file
tritici_2022_2023_geoloc_eur <- merge(geo_loc_all,sample_df_eur,by.x = "Sample.Name", by.y = "samples_eur")
coord_s_eur <- tritici_2022_2023_geoloc_eur[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_eur <- rdist.earth(coord_s_eur[,2:3], miles = FALSE, R = NULL)
geodist_matr_eur <- as.dist(dist_kms_eur)

ibd_eur <- mantel.randtest(gendist_matr_eur,geodist_matr_eur)
ibd_eur

# histogram
plot(ibd_eur, main = "Europe+_recent
     corr = 0.03788523, p = 0.055")

# density plot
dens <- kde2d(as.vector(geodist_matr_eur),as.vector(gendist_matr_eur), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(geodist_matr_eur, gendist_matr_eur, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = "Europe+_recent 
corr = 0.03788523, p = 0.055")
title(ylab="Genetic distance (no. of SNPs)", line=2, cex.lab=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_eur)~as.vector(geodist_matr_eur)),col = "red", lty=1, lwd=2)


#### Tritici europe+_2022_2023
samples_all <- readLines("../../../vcf_project_tritici/tritici_2022_2023.args")
sample_df_all = as.data.frame(samples_all)

tritici_2022_2023_gendist_all <- gendist_all[samples_all,samples_all]
gendist_matr_all <- as.dist(tritici_2022_2023_gendist_all)

tritici_2022_2023_geoloc_all <- merge(geo_loc_all,sample_df_all,by.x = "Sample.Name", by.y = "samples_all")
coord_s <- tritici_2022_2023_geoloc_all[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms <- rdist.earth(coord_s[,2:3], miles = FALSE, R = NULL)
geodist_matr_all <- as.dist(dist_kms)

#geodist_matr_all <- geodist(tritici_2022_geoloc_list_all, measure = "geodesic")

ibd_all <- mantel.randtest(gendist_matr_all,geodist_matr_all)
ibd_all

# histogram
plot(ibd_all, main = "2022+2023 tritici all
  corr = 0.4750431, p = 0.001")

# density plot
dens <- kde2d(as.vector(geodist_matr_all),as.vector(gendist_matr_all), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(geodist_matr_all, gendist_matr_all, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = "Europe+_2022_2023 
corr = 0.4750431, p = 0.001")
title(ylab="Genetic distance (no. of SNPs)", line=2, cex.lab=1)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_all)~as.vector(geodist_matr_all)),col = "red", lty=1, lwd=2)

