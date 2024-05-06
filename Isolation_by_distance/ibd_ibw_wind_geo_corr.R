## ISOLATION BY DISTANCE, ISOLATION BY WIND
#####  1.3.24: something is off with adegenet?? (https://github.com/thibautjombart/adegenet/issues/307)
#### TL;DR: use as.vector(dist objects) for functions like lm, kde2d. used to work without it, not anymore.

library(adegenet)
library(MASS)
library(fields)
library(tidyverse)

set.seed(123)

setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/isolation_by_distance/")

# read in genetic dist matrix, wind distance matrix and metadata file for all isolates

gendist_all <- as.matrix(read.csv("../dist_mat/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
geo_loc_all <- read.csv("../../../vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK9_02052024.csv")
wind_all <- as.matrix(read.csv("../dist_mat/windscape_2012-2021.wind_distance_sym.csv", header = TRUE, row.names = 1))

##par(mfrow=c(2,2))
#### TRITICI EUR+_RECENT ####

samples_eur <- readLines("../../../vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args")
sample_df_eur = as.data.frame(samples_eur)

# extract samples from distance matrix
tritici_gendist_eur <- gendist_all[samples_eur,samples_eur]
gendist_matr_eur <- as.dist(tritici_gendist_eur)

# extract coordinates from metadata file
tritici_geoloc_eur <- merge(geo_loc_all,sample_df_eur,by.x = "Sample.Name", by.y = "samples_eur")
coord_s_eur <- tritici_geoloc_eur[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_eur <- rdist.earth(coord_s_eur[,2:3], miles = FALSE, R = NULL)
geodist_matr_eur <- as.dist(dist_kms_eur)

# extract samples from wind-distance matrix
tritici_windist_eur <- wind_all[samples_eur,samples_eur]
windist_matr_eur <- as.dist(tritici_windist_eur)

## ISOLATION BY DISTANCE
ibd_eur <- mantel.randtest(gendist_matr_eur,geodist_matr_eur)
# histogram
plot(ibd_eur, main = paste("Europe+_recent \ncorr = ", signif(ibd_eur$obs, 6),", p = ",ibd_eur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_eur <- mantel.randtest(gendist_matr_eur,windist_matr_eur)
# histogram
plot(ibw_eur, main = paste("Europe+_recent \ncorr = ", signif(ibw_eur$obs, 6),", p = ",ibw_eur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND-GEO CORR
corr_geo_win_eur <- mantel.randtest(windist_matr_eur,geodist_matr_eur)
# histogram
plot(corr_geo_win_eur, main = paste("Europe+_recent \ncorr = ", signif(corr_geo_win_eur$obs, 6),", p = ",corr_geo_win_eur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

'''
#### TRITICI 2022_2023 ####

samples_2022_2023 <- readLines("../../../vcf_project_tritici/tritici_2022_2023.args")
sample_df_2022_2023 = as.data.frame(samples_2022_2023)

# extract samples from distance matrix
tritici_gendist_2022_2023 <- gendist_all[samples_2022_2023,samples_2022_2023]
gendist_matr_2022_2023 <- as.dist(tritici_gendist_2022_2023)

# extract coordinates from metadata file
tritici_geoloc_2022_2023 <- merge(geo_loc_all,sample_df_2022_2023,by.x = "Sample.Name", by.y = "samples_2022_2023")
coord_s_2022_2023 <- tritici_geoloc_2022_2023[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_2022_2023 <- rdist.earth(coord_s_2022_2023[,2:3], miles = FALSE, R = NULL)
geodist_matr_2022_2023 <- as.dist(dist_kms_2022_2023)

# extract samples from wind-distance matrix
tritici_windist_2022_2023 <- wind_all[samples_2022_2023,samples_2022_2023]
windist_matr_2022_2023 <- as.dist(tritici_windist_2022_2023)

## ISOLATION BY DISTANCE
ibd_2022_2023 <- mantel.randtest(gendist_matr_2022_2023,geodist_matr_2022_2023)
# histogram
plot(ibd_2022_2023, main = paste("Europe+_2022_2023 \ncorr = ", signif(ibd_2022_2023$obs, 6),", p = ",ibd_2022_2023$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_2022_2023 <- mantel.randtest(gendist_matr_2022_2023,windist_matr_2022_2023)
# histogram
plot(ibw_2022_2023, main = paste("Europe+_2022_2023 \ncorr = ",signif(ibw_2022_2023$obs,6),", p = ",ibw_2022_2023$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND-GEO CORR
corr_geo_win_2022_2023 <- mantel.randtest(windist_matr_2022_2023,geodist_matr_2022_2023)
# histogram
plot(corr_geo_win_2022_2023, main = paste("Europe+_2022_2023 \ncorr = ",signif(corr_geo_win_2022_2023$obs,6),", p = ",corr_geo_win_2022_2023$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

'''
#### TRITICI 2022_2023_NO_IL_TR ####

samples_2022_2023_no_il <- readLines("tritici_2022_2023_no_il.args")
sample_df_2022_2023_no_il = as.data.frame(samples_2022_2023_no_il)

# extract samples from distance matrix
tritici_gendist_2022_2023_no_il <- gendist_all[samples_2022_2023_no_il,samples_2022_2023_no_il]
gendist_matr_2022_2023_no_il <- as.dist(tritici_gendist_2022_2023_no_il)

# extract coordinates from metadata file
tritici_geoloc_2022_2023_no_il <- merge(geo_loc_all,sample_df_2022_2023_no_il,by.x = "Sample.Name", by.y = "samples_2022_2023_no_il")
coord_s_2022_2023_no_il <- tritici_geoloc_2022_2023_no_il[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_2022_2023_no_il <- rdist.earth(coord_s_2022_2023_no_il[,2:3], miles = FALSE, R = NULL)
geodist_matr_2022_2023_no_il <- as.dist(dist_kms_2022_2023_no_il)

# extract samples from wind-distance matrix
tritici_windist_2022_2023_no_il <- wind_all[samples_2022_2023_no_il,samples_2022_2023_no_il]
windist_matr_2022_2023_no_il <- as.dist(tritici_windist_2022_2023_no_il)

## ISOLATION BY DISTANCE
ibd_2022_2023_no_il <- mantel.randtest(gendist_matr_2022_2023_no_il,geodist_matr_2022_2023_no_il)
# histogram
plot(ibd_2022_2023_no_il, main = paste("Europe_2022_2023_NO_IL_TR \ncorr = ",signif(ibd_2022_2023_no_il$obs,6),", p = ",ibd_2022_2023_no_il$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_2022_2023_no_il <- mantel.randtest(gendist_matr_2022_2023_no_il,windist_matr_2022_2023_no_il)
# histogram
plot(ibw_2022_2023_no_il, main = paste("Europe_2022_2023_NO_IL_TR \ncorr = ",signif(ibw_2022_2023_no_il$obs,6),", p = ",ibw_2022_2023_no_il$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND GEO CORR
corr_geo_win_2022_2023_no_il <- mantel.randtest(windist_matr_2022_2023_no_il,geodist_matr_2022_2023_no_il)
# histogram
plot(corr_geo_win_2022_2023_no_il, main = paste("Europe_2022_2023_NO_IL_TR \ncorr = ",signif(corr_geo_win_2022_2023_no_il$obs,6),", p = ",corr_geo_win_2022_2023_no_il$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### fs level4 N_EUR 2022-2023 ####

samples_2022_2023_neur <- readLines("tritici_2022+2023_fs_level4_N_EUR.args")
sample_df_2022_2023_neur = as.data.frame(samples_2022_2023_neur)

# extract samples from distance matrix
tritici_gendist_2022_2023_neur <- gendist_all[samples_2022_2023_neur,samples_2022_2023_neur]
gendist_matr_2022_2023_neur <- as.dist(tritici_gendist_2022_2023_neur)

# extract coordinates from metadata file
tritici_geoloc_2022_2023_neur <- merge(geo_loc_all,sample_df_2022_2023_neur,by.x = "Sample.Name", by.y = "samples_2022_2023_neur")
coord_s_2022_2023_neur <- tritici_geoloc_2022_2023_neur[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_2022_2023_neur <- rdist.earth(coord_s_2022_2023_neur[,2:3], miles = FALSE, R = NULL)
geodist_matr_2022_2023_neur <- as.dist(dist_kms_2022_2023_neur)

# extract samples from wind-distance matrix
tritici_windist_2022_2023_neur <- wind_all[samples_2022_2023_neur,samples_2022_2023_neur]
windist_matr_2022_2023_neur <- as.dist(tritici_windist_2022_2023_neur)

## ISOLATION BY DISTANCE
ibd_2022_2023_neur <- mantel.randtest(gendist_matr_2022_2023_neur,geodist_matr_2022_2023_neur)
# histogram
plot(ibd_2022_2023_neur, main = paste("2022_2023_N_EUR \ncorr = ",signif(ibd_2022_2023_neur$obs,6),", p = ",ibd_2022_2023_neur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_2022_2023_neur <- mantel.randtest(gendist_matr_2022_2023_neur,windist_matr_2022_2023_neur)
# histogram
plot(ibw_2022_2023_neur, main = paste("2022_2023_N_EUR \ncorr = ",signif(ibw_2022_2023_neur$obs,6),", p = ",ibw_2022_2023_neur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND GEO CORR
corr_geo_win_2022_2023_neur <- mantel.randtest(windist_matr_2022_2023_neur,geodist_matr_2022_2023_neur)
# histogram
plot(corr_geo_win_2022_2023_neur, main = paste("2022_2023_N_EUR \ncorr = ",signif(corr_geo_win_2022_2023_neur$obs,6),", p = ",corr_geo_win_2022_2023_neur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### fs level4 S_EUR 2022-2023 ####
samples_2022_2023_seur <- readLines("tritici_2022+2023_fs_level4_S_EUR+.args")
sample_df_2022_2023_seur = as.data.frame(samples_2022_2023_seur)

# extract samples from distance matrix
tritici_gendist_2022_2023_seur <- gendist_all[samples_2022_2023_seur,samples_2022_2023_seur]
gendist_matr_2022_2023_seur <- as.dist(tritici_gendist_2022_2023_seur)

# extract coordinates from metadata file
tritici_geoloc_2022_2023_seur <- merge(geo_loc_all,sample_df_2022_2023_seur,by.x = "Sample.Name", by.y = "samples_2022_2023_seur")
coord_s_2022_2023_seur <- tritici_geoloc_2022_2023_seur[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_2022_2023_seur <- rdist.earth(coord_s_2022_2023_seur[,2:3], miles = FALSE, R = NULL)
geodist_matr_2022_2023_seur <- as.dist(dist_kms_2022_2023_seur)

# extract samples from wind-distance matrix
tritici_windist_2022_2023_seur <- wind_all[samples_2022_2023_seur,samples_2022_2023_seur]
windist_matr_2022_2023_seur <- as.dist(tritici_windist_2022_2023_seur)

## ISOLATION BY DISTANCE
ibd_2022_2023_seur <- mantel.randtest(gendist_matr_2022_2023_seur,geodist_matr_2022_2023_seur)
# histogram
plot(ibd_2022_2023_seur, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(ibd_2022_2023_seur$obs,6),", p = ",ibd_2022_2023_seur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_2022_2023_seur <- mantel.randtest(gendist_matr_2022_2023_seur,windist_matr_2022_2023_seur)
# histogram
plot(ibw_2022_2023_seur, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(ibw_2022_2023_seur$obs,6),", p = ",ibw_2022_2023_seur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND GEO CORR
corr_geo_win_2022_2023_seur <- mantel.randtest(windist_matr_2022_2023_seur,geodist_matr_2022_2023_seur)
# histogram
plot(corr_geo_win_2022_2023_seur, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(corr_geo_win_2022_2023_seur$obs,6),", p = ",corr_geo_win_2022_2023_seur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### TRITICI 2022-2023 N+S ####
samples_2022_2023_NS <- readLines("tritici_2022+2023_fs_level4_N_EUR_+_S_EUR2")
sample_df_2022_2023_NS = as.data.frame(samples_2022_2023_NS)

# extract samples from distance matrix
tritici_gendist_2022_2023_NS <- gendist_all[samples_2022_2023_NS,samples_2022_2023_NS]
gendist_matr_2022_2023_NS <- as.dist(tritici_gendist_2022_2023_NS)

# extract coordinates from metadata file
tritici_geoloc_2022_2023_NS <- merge(geo_loc_all,sample_df_2022_2023_NS,by.x = "Sample.Name", by.y = "samples_2022_2023_NS")
coord_s_2022_2023_NS <- tritici_geoloc_2022_2023_NS[, c(1,9,8)]

# calculate great circle distance in Km
dist_kms_2022_2023_NS <- rdist.earth(coord_s_2022_2023_NS[,2:3], miles = FALSE, R = NULL)
geodist_matr_2022_2023_NS <- as.dist(dist_kms_2022_2023_NS)

# extract samples from wind-distance matrix
tritici_windist_2022_2023_NS <- wind_all[samples_2022_2023_NS,samples_2022_2023_NS]
windist_matr_2022_2023_NS <- as.dist(tritici_windist_2022_2023_NS)

## ISOLATION BY DISTANCE
ibd_2022_2023_NS <- mantel.randtest(gendist_matr_2022_2023_NS,geodist_matr_2022_2023_NS)
# histogram
plot(ibd_2022_2023_neur, main = paste("2022_2023_N+S \ncorr = ",signif(ibd_2022_2023_NS$obs,6),", p = ",ibd_2022_2023_NS$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_2022_2023_neur <- mantel.randtest(gendist_matr_2022_2023_neur,windist_matr_2022_2023_neur)
# histogram
plot(ibw_2022_2023_neur, main = paste("2022_2023_N_EUR \ncorr = ",signif(ibw_2022_2023_neur$obs,6),", p = ",ibw_2022_2023_neur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND GEO CORR
corr_geo_win_2022_2023_neur <- mantel.randtest(windist_matr_2022_2023_neur,geodist_matr_2022_2023_neur)
# histogram
plot(corr_geo_win_2022_2023_neur, main = paste("2022_2023_N_EUR \ncorr = ",signif(corr_geo_win_2022_2023_neur$obs,6),", p = ",corr_geo_win_2022_2023_neur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### DENSITY PLOTS ####
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))

## ISOLATION BY DISTANCE
par(mfrow=c(2,2))

# 1. ibd density plot TRITICI EUR+_RECENT
dens_eur <- kde2d(as.vector(geodist_matr_eur),as.vector(gendist_matr_eur), n=300)
plot(geodist_matr_eur, gendist_matr_eur, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("Europe+_recent \ncorr = ",signif(ibd_eur$obs,6),", p = ",ibd_eur$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_eur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_eur)~as.vector(geodist_matr_eur)),col = "red", lty=1, lwd=2)

# 2. ibd density plot TRITICI 2022_2023 NO IL,TR
dens_2022_2023_no_il <- kde2d(as.vector(geodist_matr_2022_2023_no_il),as.vector(gendist_matr_2022_2023_no_il), n=300)
plot(geodist_matr_2022_2023_no_il, gendist_matr_2022_2023_no_il, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023_NO_IL_TR \ncorr = ",signif(ibd_2022_2023_no_il$obs,6),", p = ",ibd_2022_2023_no_il$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_2022_2023_no_il, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_no_il)~as.vector(geodist_matr_2022_2023_no_il)),col = "red", lty=1, lwd=2)

# 3.ibd density plot TRITICI 2022_2023 N_EUR
dens_2022_2023_neur <- kde2d(as.vector(geodist_matr_2022_2023_neur),as.vector(gendist_matr_2022_2023_neur), n=300)
plot(geodist_matr_2022_2023_neur, gendist_matr_2022_2023_neur, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("2022_2023_N_EUR \ncorr = ",signif(ibd_2022_2023_neur$obs,6),", p = ",ibd_2022_2023_neur$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_2022_2023_neur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_neur)~as.vector(geodist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

# 4. ibd density plot TRITICI 2022_2023 S_EUR2
dens_2022_2023_seur <- kde2d(as.vector(geodist_matr_2022_2023_seur),as.vector(gendist_matr_2022_2023_seur), n=300)
plot(geodist_matr_2022_2023_seur, gendist_matr_2022_2023_seur, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(ibd_2022_2023_seur$obs,6),", p = ",ibd_2022_2023_seur$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_2022_2023_seur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_seur)~as.vector(geodist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)


## ISOLATION BY WIND
par(mfrow=c(2,2))

# 1. ibw density plot TRITICI EUR+_RECENT
dens_wind_eur <- kde2d(as.vector(windist_matr_eur),as.vector(gendist_matr_eur), n=300)
plot(windist_matr_eur, gendist_matr_eur, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_recent \ncorr = ",signif(ibw_eur$obs,6),", p = ",ibw_eur$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_wind_eur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_eur)~as.vector(windist_matr_eur)),col = "red", lty=1, lwd=2)

# 2. ibw density plot TRITICI 2022_2023 NO IL,TR
dens_wind_2022_2023_no_il <- kde2d(as.vector(windist_matr_2022_2023_no_il),as.vector(gendist_matr_2022_2023_no_il), n=300)
plot(windist_matr_2022_2023_no_il, gendist_matr_2022_2023_no_il, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe_2022_2023_NO_IL_TR \ncorr = ",signif(ibw_2022_2023_no_il$obs,6),", p = ",ibw_2022_2023_no_il$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_wind_2022_2023_no_il, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_no_il)~as.vector(windist_matr_2022_2023_no_il)),col = "red", lty=1, lwd=2)

# 3. ibw density plot TRITICI 2022_2023 N_EUR
dens_wind_2022_2023_neur <- kde2d(as.vector(windist_matr_2022_2023_neur),as.vector(gendist_matr_2022_2023_neur), n=300)
plot(windist_matr_2022_2023_neur, gendist_matr_2022_2023_neur, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("2022_2023_N_EUR \ncorr = ",signif(ibw_2022_2023_neur$obs,6),", p = ",ibw_2022_2023_neur$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_wind_2022_2023_neur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_neur)~as.vector(windist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

# 4. ibw density plot TRITICI 2022_2023 S_EUR2
dens_wind_2022_2023_seur <- kde2d(as.vector(windist_matr_2022_2023_seur),as.vector(gendist_matr_2022_2023_seur), n=300)
plot(windist_matr_2022_2023_seur, gendist_matr_2022_2023_seur, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(ibw_2022_2023_seur$obs,6),", p = ",ibw_2022_2023_seur$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_wind_2022_2023_seur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_seur)~as.vector(windist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)


# WIND GEO CORR
par(mfrow=c(2,2))

# 1. wind geo corr density plot TRITICI EUR+_RECENT
dens_corr_eur <- kde2d(as.vector(windist_matr_eur),as.vector(geodist_matr_eur), n=300)
plot(windist_matr_eur, geodist_matr_eur, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_recent \ncorr = ",signif(corr_geo_win_eur$obs,6),", p = ",corr_geo_win_eur$pvalue))
title(ylab="Geographic distance (km)", line=2, cex.lab=1)
image(dens_corr_eur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(geodist_matr_eur)~as.vector(windist_matr_eur)),col = "red", lty=1, lwd=2)

# 2. wind geo corr density plot TRITICI 2022_2023 NO IL,TR
dens_corr_2022_2023_no_il <- kde2d(as.vector(windist_matr_2022_2023_no_il),as.vector(geodist_matr_2022_2023_no_il), n=300)
plot(windist_matr_2022_2023_no_il, geodist_matr_2022_2023_no_il, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023_NO_IL_TR \ncorr = ",signif(corr_geo_win_2022_2023_no_il$obs,6),", p = ",corr_geo_win_2022_2023_no_il$pvalue))
title(ylab="Geographic distance (km)", line=2, cex.lab=1)
image(dens_corr_2022_2023_no_il, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(geodist_matr_2022_2023_no_il)~as.vector(windist_matr_2022_2023_no_il)),col = "red", lty=1, lwd=2)

# 3. wind geo corr density plot TRITICI 2022_2023 N_EUR
dens_corr_2022_2023_neur <- kde2d(as.vector(windist_matr_2022_2023_neur),as.vector(geodist_matr_2022_2023_neur), n=300)
plot(windist_matr_2022_2023_neur, geodist_matr_2022_2023_neur, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("2022_2023_N_EUR \ncorr = ",signif(corr_geo_win_2022_2023_neur$obs,6),", p = ",corr_geo_win_2022_2023_neur$pvalue))
title(ylab="Geographic distance (km)", line=2, cex.lab=1)
image(dens_corr_2022_2023_neur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(geodist_matr_2022_2023_neur)~as.vector(windist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

# 4. wind geo corr density plot TRITICI 2022_2023 S_EUR
dens_corr_2022_2023_seur <- kde2d(as.vector(windist_matr_2022_2023_seur),as.vector(geodist_matr_2022_2023_seur), n=300)
plot(windist_matr_2022_2023_seur, geodist_matr_2022_2023_seur, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(corr_geo_win_2022_2023_seur$obs,6),", p = ",corr_geo_win_2022_2023_seur$pvalue))
title(ylab="Geographic distance (km)", line=2, cex.lab=1)
image(dens_corr_2022_2023_seur, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(geodist_matr_2022_2023_seur)~as.vector(windist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)


### density plots tritici 2022-2023 incl tr, il
# ibd density plot TRITICI 2022_2023
dens_2022_2023 <- kde2d(as.vector(geodist_matr_2022_2023),as.vector(gendist_matr_2022_2023), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(geodist_matr_2022_2023, gendist_matr_2022_2023, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023 \ncorr = ", signif(ibd_2022_2023$obs, 6),", p = ",ibd_2022_2023$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_2022_2023, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(geodist_matr_2022_2023)),col = "red", lty=1, lwd=2)

# ibw density plot TRITICI 2022_2023
dens_wind_2022_2023 <- kde2d(as.vector(windist_matr_2022_2023),as.vector(gendist_matr_2022_2023), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(windist_matr_2022_2023, gendist_matr_2022_2023, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023 \ncorr = ", signif(ibw_2022_2023$obs, 6),", p = ",ibw_2022_2023$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_wind_2022_2023, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(windist_matr_2022_2023)),col = "red", lty=1, lwd=2)

# wind geo corr density plot TRITICI 2022_2023
dens_corr_2022_2023 <- kde2d(as.vector(windist_matr_2022_2023),as.vector(geodist_matr_2022_2023), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(windist_matr_2022_2023, geodist_matr_2022_2023, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023  \ncorr = ", signif(corr_geo_win_2022_2023$obs, 6),", p = ",corr_geo_win_2022_2023$pvalue))
title(ylab="Geographic distance (km)", line=2, cex.lab=1)
image(dens_corr_2022_2023, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(geodist_matr_2022_2023)~as.vector(windist_matr_2022_2023)),col = "red", lty=1, lwd=2)


# TRITICI N+S 2022_2023
# ibd
dens_NS <- kde2d(as.vector(geodist_matr_2022_2023_N),as.vector(gendist_matr_2022_2023), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(geodist_matr_2022_2023, gendist_matr_2022_2023, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023 \ncorr = ", signif(ibd_2022_2023$obs, 6),", p = ",ibd_2022_2023$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_2022_2023, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(geodist_matr_2022_2023)),col = "red", lty=1, lwd=2)
