## ISOLATION BY DISTANCE, ISOLATION BY WIND, ISOLATION BY CLIMATE
#####  1.3.24: something is off with adegenet?? (https://github.com/thibautjombart/adegenet/issues/307)
#### TL;DR: use as.vector(dist objects) for functions like lm, kde2d. used to work without it, not anymore.

library(tidyverse, lib.loc = "/usr/local/lib/R/site-library")
library(adegenet)
library(MASS)
library(fields)
library(vegan)

set.seed(123)

setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/isolation_by_distance/")

# read in all distance matrices for all isolates

gendist_all <- as.matrix(read.csv("../dist_mat/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
wind_all <- as.matrix(read.csv("../dist_mat/windscape_2012-2021.wind_distance_sym.csv", header = TRUE, row.names = 1))
geodist_all <- as.matrix(read.csv("../dist_mat/geo_dist_matrix_tritici_europe_recent.csv", header = TRUE, row.names = 1))
clim_all <- as.matrix(read.csv("../dist_mat/clim_dist_unweighted_7pc.csv", header = TRUE, row.names = 1))

#### TRITICI 2022_2023 ####
samples_2022_2023 <- readLines("tritici_2022_2023.args")

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

samples_2022_2023_neur <- readLines("tritici_2022+2023_fs_level4_N_EUR.args")

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

samples_2022_2023_seur <- readLines("tritici_2022+2023_fs_level4_S_EUR+.args")

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



#### E-W / N-S IBD ####

# n_eur
geo_dist_ew <- as.matrix(read.csv("geo_dist_X_n_eur_2022-2023.csv", header = TRUE, row.names = 1))
geo_dist_ns <- as.matrix(read.csv("geo_dist_Y_n_eur_2022-2023.csv", header = TRUE, row.names = 1))
samples_2022_2023_neur <- readLines("tritici_2022+2023_fs_level4_N_EUR.args")

# s_eur2
geo_dist_seur_ew <- as.matrix(read.csv("../dist_mat/geo_dist_X_s_eur2_2022-2023.csv", header = TRUE, row.names=1))
geo_dist_seur_ns <- as.matrix(read.csv("../dist_mat/geo_dist_Y_s_eur2_2022-2023.csv", header = TRUE, row.names=1))
samples_2022_2023_seur <- readLines("tritici_2022+2023_fs_level4_S_EUR+.args")


# extract samples from matrices and get in correct order
# n_eur
gendist_matr_neur <- as.dist(gendist_all[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_ew_neur <- as.dist(geo_dist_ew[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_ns_neur <- as.dist(geo_dist_ns[samples_2022_2023_neur,samples_2022_2023_neur])

#s_eur2
gendist_matr_seur <- as.dist(gendist_all[samples_2022_2023_seur,samples_2022_2023_seur])
geodist_matr_ew_seur <- as.dist(geo_dist_seur_ew[samples_2022_2023_seur,samples_2022_2023_seur])
geodist_matr_ns_seur <- as.dist(geo_dist_seur_ns[samples_2022_2023_seur,samples_2022_2023_seur])


## ISOLATION BY DISTANCE

#n_eur

ibd_ew <- mantel.randtest(gendist_matr_neur, geodist_matr_ew_neur)
ibd_ew

ibd_ns <- mantel.randtest(gendist_matr_neur, geodist_matr_ns_neur)
ibd_ns

#s_eur2
ibd_ew_seur <- mantel.randtest(gendist_matr_seur, geodist_matr_ew_seur)
ibd_ew_seur

ibd_ns_seur <- mantel.randtest(gendist_matr_seur, geodist_matr_ns_seur)
ibd_ns_seur

myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))

ind_subs <- runif(1000,1,15400)
ind_seur <- runif(1000,1,1711)

## ISOLATION BY DISTANCE plots
par(mfrow=c(2,2))

#N_EUR
# 1. ibd density plot e-w
dens_ew <- kde2d(as.vector(geodist_matr_ew_neur),as.vector(gendist_matr_neur), n=200)
plot(geodist_matr_ew_neur[ind_subs], gendist_matr_neur[ind_subs], pch=20,cex=.5,xlab="E-W geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0)
image(dens_ew, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_neur)~as.vector(geodist_matr_ew_neur)),col = "red", lty=1, lwd=2)

# 2. ibd density plot n-s
dens_ns <- kde2d(as.vector(geodist_matr_ns_neur),as.vector(gendist_matr_neur), n=200)
plot(geodist_matr_ns_neur[ind_subs], gendist_matr_neur[ind_subs], pch=20,cex=.5,xlab="N-S geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0)
image(dens_ns, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_neur)~as.vector(geodist_matr_ns_neur)),col = "red", lty=1, lwd=2)

#S_EUR2
# 1. ibd density plot e-w
dens_ew_seur <- kde2d(as.vector(geodist_matr_ew_seur),as.vector(gendist_matr_seur), n=200)
plot(geodist_matr_ew_seur[ind_seur], gendist_matr_seur[ind_seur], pch=20,cex=.5,xlab="E-W geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0)
image(dens_ew_seur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_seur)~as.vector(geodist_matr_ew_seur)),col = "red", lty=1, lwd=2)

# 2. ibd density plot n-s
dens_ns_seur <- kde2d(as.vector(geodist_matr_ns_seur),as.vector(gendist_matr_seur), n=200)
plot(geodist_matr_ns_seur[ind_seur], gendist_matr_seur[ind_seur], pch=20,cex=.5,xlab="N-S geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0)
image(dens_ns_seur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_seur)~as.vector(geodist_matr_ns_seur)),col = "red", lty=1, lwd=2)
