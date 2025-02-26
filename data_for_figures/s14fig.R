library(tidyverse)
library(adegenet)
library(MASS)
library(fields)
library(vegan)

set.seed(123)

# read in all distance matrices for all isolates

gendist_all <- as.matrix(read.csv("../distance_matrix/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
wind_all <- as.matrix(read.csv("../windscape/windscape_2012-2021.wind_distance_sym.csv", header = TRUE, row.names = 1))
geodist_all <- as.matrix(read.csv("../Isolation_by_distance/geo_dist_matrix_tritici_europe_recent.csv", header = TRUE, row.names = 1))
clim_all <- as.matrix(read.csv("../Isolation_by_distance/clim_dist_unweighted_7pc.csv", header = TRUE, row.names = 1))

#### TRITICI 2022_2023 ####
samples_2022_2023 <- readLines("../Datasets/tritici_2022_2023.args")

# extract samples from distance matrices
gendist_matr_2022_2023 <- as.dist(gendist_all[samples_2022_2023,samples_2022_2023])
geodist_matr_2022_2023 <- as.dist(geodist_all[samples_2022_2023,samples_2022_2023])
windist_matr_2022_2023 <- as.dist(wind_all[samples_2022_2023,samples_2022_2023])
climdist_matr_2022_2023 <- as.dist(clim_all[samples_2022_2023, samples_2022_2023])

## IBD, IBW, IBC
ibd_2022_2023 <- mantel.randtest(gendist_matr_2022_2023,geodist_matr_2022_2023)
ibw_2022_2023 <- mantel.randtest(gendist_matr_2022_2023, windist_matr_2022_2023)
ibc_2022_2023 <- mantel.randtest(gendist_matr_2022_2023, climdist_matr_2022_2023)

## correlations with geography
corr_wg_2022_2023 <- mantel.randtest(geodist_matr_2022_2023, windist_matr_2022_2023)
corr_cg_2022_2023 <- mantel.randtest(geodist_matr_2022_2023, climdist_matr_2022_2023)

#### fs level4 N_EUR 2022-2023 ####

samples_2022_2023_neur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_N_EUR.args")

# extract samples from 3 matrices
gendist_matr_2022_2023_neur <- as.dist(gendist_all[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_2022_2023_neur <- as.dist(geodist_all[samples_2022_2023_neur,samples_2022_2023_neur])
windist_matr_2022_2023_neur <- as.dist(wind_all[samples_2022_2023_neur,samples_2022_2023_neur])
climdist_matr_2022_2023_neur <- as.dist(clim_all[samples_2022_2023_neur, samples_2022_2023_neur])


## IBD, IBW, IBC
ibd_2022_2023_neur <- mantel.randtest(gendist_matr_2022_2023_neur,geodist_matr_2022_2023_neur)
ibw_2022_2023_neur <- mantel.randtest(gendist_matr_2022_2023_neur,windist_matr_2022_2023_neur)
ibc_2022_2023_neur <- mantel.randtest(gendist_matr_2022_2023_neur,climdist_matr_2022_2023_neur)

## correlations with geography
corr_wg_2022_2023_neur <- mantel.randtest(geodist_matr_2022_2023_neur, windist_matr_2022_2023_neur)
corr_cg_2022_2023_neur <- mantel.randtest(geodist_matr_2022_2023_neur, climdist_matr_2022_2023_neur)

#### fs level4 S_EUR2 2022-2023 ####

samples_2022_2023_seur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_S_EUR+.args")

# extract samples from 3 matrices
gendist_matr_2022_2023_seur <- as.dist(gendist_all[samples_2022_2023_seur,samples_2022_2023_seur])
geodist_matr_2022_2023_seur <- as.dist(geodist_all[samples_2022_2023_seur,samples_2022_2023_seur])
windist_matr_2022_2023_seur <- as.dist(wind_all[samples_2022_2023_seur,samples_2022_2023_seur])
climdist_matr_2022_2023_seur <- as.dist(clim_all[samples_2022_2023_seur, samples_2022_2023_seur])


## IBD, IBW, IBC
ibd_2022_2023_seur <- mantel.randtest(gendist_matr_2022_2023_seur,geodist_matr_2022_2023_seur)
ibw_2022_2023_seur <- mantel.randtest(gendist_matr_2022_2023_seur,windist_matr_2022_2023_seur)
ibc_2022_2023_seur <- mantel.randtest(gendist_matr_2022_2023_seur,climdist_matr_2022_2023_seur)

## correlations with geography
corr_wg_2022_2023_seur <- mantel.randtest(geodist_matr_2022_2023_seur, windist_matr_2022_2023_seur)
corr_cg_2022_2023_seur <- mantel.randtest(geodist_matr_2022_2023_seur, climdist_matr_2022_2023_seur)

#### DENSITY PLOTS ####

##### only plot 1000 random data points to have manageable image sizes
indices_all <- runif(1000,1,length(gendist_matr_2022_2023))
indices_neur <- runif(1000,1,length(gendist_matr_2022_2023_neur))
indices_seur <- runif(1000,1,length(gendist_matr_2022_2023_seur))

myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))

par(mfcol = c(3,3))
par(oma=c(5,4,4,1)+0.1)

## ISOLATION BY DISTANCE

#1. IBD-2022-2023
dens_2022_2023 <- kde2d(as.vector(geodist_matr_2022_2023),as.vector(gendist_matr_2022_2023), n=200)
plot(geodist_matr_2022_2023[indices_all], gendist_matr_2022_2023[indices_all], pch=20,cex=.5,xlab="Geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,3000), ylim = c(0.0014,0.0025))
image(dens_2022_2023, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(geodist_matr_2022_2023)),col = "red", lty=1, lwd=2)

#2. IBD-2022-2023-NEUR
dens_2022_2023_neur <- kde2d(as.vector(geodist_matr_2022_2023_neur),as.vector(gendist_matr_2022_2023_neur), n=200)
plot(geodist_matr_2022_2023_neur[indices_neur], gendist_matr_2022_2023_neur[indices_neur], pch=20,cex=.5,xlab="Geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,3000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_neur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_neur)~as.vector(geodist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

#3. IBD-2022-2023-SEUR2 
dens_2022_2023_seur <- kde2d(as.vector(geodist_matr_2022_2023_seur),as.vector(gendist_matr_2022_2023_seur), n=200)
plot(geodist_matr_2022_2023_seur[indices_seur], gendist_matr_2022_2023_seur[indices_seur], pch=20,cex=.5,xlab="Geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,3000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_seur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_seur)~as.vector(geodist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)


## ISOLATION BY WIND

#4. IBW-2022-2023
dens_2022_2023_w <- kde2d(as.vector(windist_matr_2022_2023),as.vector(gendist_matr_2022_2023), n=200)
plot(windist_matr_2022_2023[indices_all], gendist_matr_2022_2023[indices_all], pch=20,cex=.5,xlab="Wind distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0,xlim = c(0,2000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_w, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(windist_matr_2022_2023)),col = "red", lty=1, lwd=2)

#5. IBW-2022-2023-NEUR
dens_2022_2023_w_neur <- kde2d(as.vector(windist_matr_2022_2023_neur),as.vector(gendist_matr_2022_2023_neur), n=200)
plot(windist_matr_2022_2023_neur[indices_neur], gendist_matr_2022_2023_neur[indices_neur], pch=20,cex=.5,xlab="Wind distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0,xlim = c(0,2000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_w_neur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_neur)~as.vector(windist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

#6. IBW-2022-2023-SEUR2
dens_2022_2023_w_seur <- kde2d(as.vector(windist_matr_2022_2023_seur),as.vector(gendist_matr_2022_2023_seur), n=200)
plot(windist_matr_2022_2023_seur[indices_seur], gendist_matr_2022_2023_seur[indices_seur], pch=20,cex=.5,xlab="Wind distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0,xlim = c(0,2000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_w_seur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_seur)~as.vector(windist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)

## ISOLATION BY CLIMATE

#7. IBC-2022-2023
dens_2022_2023_c <- kde2d(as.vector(climdist_matr_2022_2023),as.vector(gendist_matr_2022_2023), n=200)
plot(climdist_matr_2022_2023[indices_all], gendist_matr_2022_2023[indices_all], pch=20,cex=.5,xlab="Climate distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,8),ylim = c(0.0014,0.0025))
image(dens_2022_2023_c, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(climdist_matr_2022_2023)),col = "red", lty=1, lwd=2)

#8. IBC-2022-2023-NEUR
dens_2022_2023_c_neur <- kde2d(as.vector(climdist_matr_2022_2023_neur),as.vector(gendist_matr_2022_2023_neur), n=200)
plot(climdist_matr_2022_2023_neur[indices_neur], gendist_matr_2022_2023_neur[indices_neur], pch=20,cex=.5,xlab="Climate distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,8),ylim = c(0.0014,0.0025))
image(dens_2022_2023_c_neur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_neur)~as.vector(climdist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

#9. #IBC-2022-2023-SEUR2
dens_2022_2023_c_seur <- kde2d(as.vector(climdist_matr_2022_2023_seur),as.vector(gendist_matr_2022_2023_seur), n=200)
plot(climdist_matr_2022_2023_seur[indices_seur], gendist_matr_2022_2023_seur[indices_seur], pch=20,cex=.5,xlab="Climate distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,8),ylim = c(0.0014,0.0025))
image(dens_2022_2023_c_seur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_seur)~as.vector(climdist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)
