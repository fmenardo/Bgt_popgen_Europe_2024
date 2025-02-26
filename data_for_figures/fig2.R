#### Fig 2 ####
library(tidyverse)
library(adegenet)
library(MASS)
library(fields)
library(vegan)

set.seed(123)

# read in all distance matrices for all isolates

gendist_all <- as.matrix(read.csv("../distance_matrix/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
geodist_all <- as.matrix(read.csv("../Isolation_by_distance/geo_dist_matrix_tritici_europe_recent.csv", header = TRUE, row.names = 1))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))

## Fig 2a - N_EUR

samples_2022_2023_neur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_N_EUR.args")
gendist_matr_2022_2023_neur <- as.dist(gendist_all[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_2022_2023_neur <- as.dist(geodist_all[samples_2022_2023_neur,samples_2022_2023_neur])
indices_neur <- runif(1000,1,length(gendist_matr_2022_2023_neur)) ##### only plot 1000 random data points to have manageable image sizes
dens_2022_2023_neur <- kde2d(as.vector(geodist_matr_2022_2023_neur),as.vector(gendist_matr_2022_2023_neur), n=200)
plot(geodist_matr_2022_2023_neur[indices_neur], gendist_matr_2022_2023_neur[indices_neur], pch=20,cex=.5,xlab="Geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,3000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_neur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_neur)~as.vector(geodist_matr_2022_2023_neur)),col = "red", lty=1, lwd=2)

## Fig 2b - S_EUR2
samples_2022_2023_seur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_S_EUR+.args")
gendist_matr_2022_2023_seur <- as.dist(gendist_all[samples_2022_2023_seur,samples_2022_2023_seur])
geodist_matr_2022_2023_seur <- as.dist(geodist_all[samples_2022_2023_seur,samples_2022_2023_seur])
indices_seur <- runif(1000,1,length(gendist_matr_2022_2023_seur)) ##### only plot 1000 random data points to have manageable image sizes
dens_2022_2023_seur <- kde2d(as.vector(geodist_matr_2022_2023_seur),as.vector(gendist_matr_2022_2023_seur), n=200)
plot(geodist_matr_2022_2023_seur[indices_seur], gendist_matr_2022_2023_seur[indices_seur], pch=20,cex=.5,xlab="Geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,3000), ylim = c(0.0014,0.0025))
image(dens_2022_2023_seur, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023_seur)~as.vector(geodist_matr_2022_2023_seur)),col = "red", lty=1, lwd=2)

## Fig 2c - all-2022_2023
samples_2022_2023 <- readLines("../Datasets/tritici_2022_2023.args")
gendist_matr_2022_2023 <- as.dist(gendist_all[samples_2022_2023,samples_2022_2023])
geodist_matr_2022_2023 <- as.dist(geodist_all[samples_2022_2023,samples_2022_2023])
indices_all <- runif(1000,1,length(gendist_matr_2022_2023)) ##### only plot 1000 random data points to have manageable image sizes
dens_2022_2023 <- kde2d(as.vector(geodist_matr_2022_2023),as.vector(gendist_matr_2022_2023), n=200)
plot(geodist_matr_2022_2023[indices_all], gendist_matr_2022_2023[indices_all], pch=20,cex=.5,xlab="Geographic distance", 
     ylab="Genetic distance", cex.lab = 1.6, cex.axis=2.0, xlim = c(0,3000), ylim = c(0.0014,0.0025))
image(dens_2022_2023, col=transp(myPal(200),.7), add=TRUE)
abline(lm(as.vector(gendist_matr_2022_2023)~as.vector(geodist_matr_2022_2023)),col = "red", lty=1, lwd=2)


## Fig 2d - FEEMS plot is one of the outputs of the script ../FEEMS/run_feems.py. For details, see ../FEEMS/FEEMS.md



