## ISOLATION BY DISTANCE, ISOLATION BY WIND
#####  1.3.24: something is off with adegenet?? (https://github.com/thibautjombart/adegenet/issues/307)
#### TL;DR: use as.vector(dist objects) for functions like lm, kde2d. used to work without it, not anymore.

library(adegenet)
library(MASS)
library(fields)
library(tidyverse)
library(vegan)

set.seed(123)

setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/isolation_by_distance/")

# read in genetic dist matrix, wind distance matrix and metadata file for all isolates

gendist_all <- as.matrix(read.csv("../dist_mat/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
wind_all <- as.matrix(read.csv("../dist_mat/windscape_2012-2021.wind_distance_sym.csv", header = TRUE, row.names = 1))
geodist_all <- as.matrix(read.csv("../dist_mat/geo_dist_matrix_tritici_europe_recent.csv", header = TRUE, row.names = 1))


#### TRITICI EUR+_RECENT ####

samples_eur <- readLines("../../../vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args")

# extract samples from the 3 distance matrices 
gendist_matr_eur <- as.dist(gendist_all[samples_eur,samples_eur])
geodist_matr_eur <- as.dist(geodist_all[samples_eur,samples_eur])
windist_matr_eur <- as.dist(wind_all[samples_eur,samples_eur])

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
corr_geo_win_eur <- mantel.randtest(geodist_matr_eur,windist_matr_eur)
# histogram
plot(corr_geo_win_eur, main = paste("Europe+_recent \ncorr = ", signif(corr_geo_win_eur$obs, 6),", p = ",corr_geo_win_eur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

#### TRITICI 2022_2023_NO_IL_TR ####

samples_2022_2023_no_il <- readLines("tritici_2022_2023_no_il.args")

# extract samples from distance matrices
gendist_matr_2022_2023_no_il <- as.dist(gendist_all[samples_2022_2023_no_il,samples_2022_2023_no_il])
geodist_matr_2022_2023_no_il <- as.dist(geodist_all[samples_2022_2023_no_il,samples_2022_2023_no_il])
windist_matr_2022_2023_no_il <- as.dist(wind_all[samples_2022_2023_no_il,samples_2022_2023_no_il])

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
corr_geo_win_2022_2023_no_il <- mantel.randtest(geodist_matr_2022_2023_no_il,windist_matr_2022_2023_no_il)
# histogram
plot(corr_geo_win_2022_2023_no_il, main = paste("Europe_2022_2023_NO_IL_TR \ncorr = ",signif(corr_geo_win_2022_2023_no_il$obs,6),", p = ",corr_geo_win_2022_2023_no_il$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### fs level4 N_EUR 2022-2023 ####

samples_2022_2023_neur <- readLines("tritici_2022+2023_fs_level4_N_EUR.args")

# extract samples from 3 matrices
gendist_matr_2022_2023_neur <- as.dist(gendist_all[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_2022_2023_neur <- as.dist(geodist_all[samples_2022_2023_neur,samples_2022_2023_neur])
windist_matr_2022_2023_neur <- as.dist(wind_all[samples_2022_2023_neur,samples_2022_2023_neur])

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
corr_geo_win_2022_2023_neur <- mantel.randtest(geodist_matr_2022_2023_neur,windist_matr_2022_2023_neur)
# histogram
plot(corr_geo_win_2022_2023_neur, main = paste("2022_2023_N_EUR \ncorr = ",signif(corr_geo_win_2022_2023_neur$obs,6),", p = ",corr_geo_win_2022_2023_neur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### fs level4 S_EUR 2022-2023 ####
samples_2022_2023_seur <- readLines("tritici_2022+2023_fs_level4_S_EUR+.args")

# extract samples from distance matrices
gendist_matr_2022_2023_seur <- as.dist(gendist_all[samples_2022_2023_seur, samples_2022_2023_seur])
geodist_matr_2022_2023_seur <- as.dist(geodist_all[samples_2022_2023_seur, samples_2022_2023_seur])
windist_matr_2022_2023_seur <- as.dist(wind_all[samples_2022_2023_seur, samples_2022_2023_seur])

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
corr_geo_win_2022_2023_seur <- mantel.randtest(geodist_matr_2022_2023_seur,windist_matr_2022_2023_seur)
# histogram
plot(corr_geo_win_2022_2023_seur, main = paste("2022_2023_S_EUR2 \ncorr = ",signif(corr_geo_win_2022_2023_seur$obs,6),", p = ",corr_geo_win_2022_2023_seur$pvalue), xlab = "")
title(xlab = "Permuted correlation values")


#### TRITICI 2022-2023 N+S ####
samples_2022_2023_NS <- readLines("tritici_2022+2023_fs_level4_N_EUR_+_S_EUR2")

# extract samples from distance matrices 
gendist_matr_2022_2023_NS <- as.dist(gendist_all[samples_2022_2023_NS,samples_2022_2023_NS])
geodist_matr_2022_2023_NS <- as.dist(geodist_all[samples_2022_2023_NS,samples_2022_2023_NS])
windist_matr_2022_2023_NS <- as.dist(wind_all[samples_2022_2023_NS,samples_2022_2023_NS])

## ISOLATION BY DISTANCE
ibd_2022_2023_NS <- mantel.randtest(gendist_matr_2022_2023_NS,geodist_matr_2022_2023_NS)
# histogram
plot(ibd_2022_2023_NS, main = paste("2022_2023_N+S \ncorr = ",signif(ibd_2022_2023_NS$obs,6),", p = ",ibd_2022_2023_NS$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## ISOLATION BY WIND
ibw_2022_2023_NS <- mantel.randtest(gendist_matr_2022_2023_NS,windist_matr_2022_2023_NS)
# histogram
plot(ibw_2022_2023_NS, main = paste("2022_2023_N+S \ncorr = ",signif(ibw_2022_2023_NS$obs,6),", p = ",ibw_2022_2023_NS$pvalue), xlab = "")
title(xlab = "Permuted correlation values")

## WIND GEO CORR
corr_geo_win_2022_2023_NS <- mantel.randtest(geodist_matr_2022_2023_NS,windist_matr_2022_2023_NS)
# histogram
plot(corr_geo_win_2022_2023_NS, main = paste("2022_2023_N+S \ncorr = ",signif(corr_geo_win_2022_2023_NS$obs,6),", p = ",corr_geo_win_2022_2023_NS$pvalue), xlab = "")
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
par(mfrow=c(2,2))
# ibd
dens_NS <- kde2d(as.vector(geodist_matr_2022_2023_NS),as.vector(gendist_matr_2022_2023_NS), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(geodist_matr_2022_2023_NS, gendist_matr_2022_2023_NS, pch=20,cex=.5,xlab="Geographic distance (km)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023_N+S \ncorr = ", signif(ibd_2022_2023_NS$obs, 6),", p = ",ibd_2022_2023_NS$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_NS, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_NS)~as.vector(geodist_matr_2022_2023_NS)),col = "red", lty=1, lwd=2)

dens_wind_NS <- kde2d(as.vector(windist_matr_2022_2023_NS),as.vector(gendist_matr_2022_2023_NS), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(windist_matr_2022_2023_NS, gendist_matr_2022_2023_NS, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023_N+S \ncorr = ", signif(ibw_2022_2023_NS$obs, 6),", p = ",ibw_2022_2023_NS$pvalue))
title(ylab="Genetic distance (normalised no. of SNPs)", line=2, cex.lab=1)
image(dens_wind_NS, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_NS)~as.vector(windist_matr_2022_2023_NS)),col = "red", lty=1, lwd=2)

dens_corr_NS <- kde2d(as.vector(windist_matr_2022_2023_NS),as.vector(geodist_matr_2022_2023_NS), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(windist_matr_2022_2023_NS, geodist_matr_2022_2023_NS, pch=20,cex=.5,xlab="Average wind distance (wind hours)", 
     ylab="", cex.lab = 1, main = paste("Europe+_2022_2023_N+S \ncorr = ", signif(corr_geo_win_2022_2023_NS$obs, 6),", p = ",corr_geo_win_2022_2023_NS$pvalue))
title(ylab="Geographic distance (km)", line=2, cex.lab=1)
image(dens_corr_NS, col=transp(myPal(300),.7), add=TRUE)
abline(lm(as.vector(geodist_matr_2022_2023_NS)~as.vector(windist_matr_2022_2023_NS)),col = "red", lty=1, lwd=2)


#### Partial mantel test 

# EUR
mantel.partial(gendist_matr_eur, geodist_matr_eur, windist_matr_eur, method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

mantel.partial(gendist_matr_eur, windist_matr_eur, geodist_matr_eur,  method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

# 2022_2023_no_il_tr
mantel.partial(gendist_matr_2022_2023_no_il, geodist_matr_2022_2023_no_il, windist_matr_2022_2023_no_il, method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

mantel.partial(gendist_matr_2022_2023_no_il, windist_matr_2022_2023_no_il, geodist_matr_2022_2023_no_il,  method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

# 2022_2023 N_eur
mantel.partial(gendist_matr_2022_2023_neur, geodist_matr_2022_2023_neur, windist_matr_2022_2023_neur, method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

mantel.partial(gendist_matr_2022_2023_neur, windist_matr_2022_2023_neur, geodist_matr_2022_2023_neur,  method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

# 2022_2023 S_eur2
mantel.partial(gendist_matr_2022_2023_seur, geodist_matr_2022_2023_seur, windist_matr_2022_2023_seur, method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

mantel.partial(gendist_matr_2022_2023_seur, windist_matr_2022_2023_seur, geodist_matr_2022_2023_seur,  method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

# 2022_2023 N_eur + S_eur2
mantel.partial(gendist_matr_2022_2023_NS, geodist_matr_2022_2023_NS, windist_matr_2022_2023_NS, method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

mantel.partial(gendist_matr_2022_2023_NS, windist_matr_2022_2023_NS, geodist_matr_2022_2023_NS,  method = "pearson", permutations = 999, 
               strata = NULL, na.rm = FALSE)

