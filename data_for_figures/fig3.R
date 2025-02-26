#### Fig 3 ####
library(tidyverse)
library(adegenet)
library(MASS)
library(fields)
library(vegan)

set.seed(123)

## Fig 3 a,b : Output of the script ../windscape/windscape_plot.R. For details, see ../windscape/windscape.md 

## Fig 3c,d - IBD along different axes, N_EUR, S_EUR2

gendist_all <- as.matrix(read.csv("../distance_matrix/gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))

# n_eur
geo_dist_ew <- as.matrix(read.csv("../Isolation_by_distance/geo_dist_X_n_eur_2022-2023.csv", header = TRUE, row.names = 1))
geo_dist_ns <- as.matrix(read.csv("./Isolation_by_distance/geo_dist_Y_n_eur_2022-2023.csv", header = TRUE, row.names = 1))
samples_2022_2023_neur <- readLines("./Isolation_by_distance/tritici_2022+2023_fs_level4_N_EUR.args")

# s_eur2
geo_dist_seur_ew <- as.matrix(read.csv("./Isolation_by_distance/geo_dist_X_s_eur2_2022-2023.csv", header = TRUE, row.names=1))
geo_dist_seur_ns <- as.matrix(read.csv("./Isolation_by_distance/geo_dist_Y_s_eur2_2022-2023.csv", header = TRUE, row.names=1))
samples_2022_2023_seur <- readLines("../Isolation_by_distance/tritici_2022+2023_fs_level4_S_EUR+.args")


# extract samples from matrices and get in correct order
# n_eur
gendist_matr_neur <- as.dist(gendist_all[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_ew_neur <- as.dist(geo_dist_ew[samples_2022_2023_neur,samples_2022_2023_neur])
geodist_matr_ns_neur <- as.dist(geo_dist_ns[samples_2022_2023_neur,samples_2022_2023_neur])

#s_eur2
gendist_matr_seur <- as.dist(gendist_all[samples_2022_2023_seur,samples_2022_2023_seur])
geodist_matr_ew_seur <- as.dist(geo_dist_seur_ew[samples_2022_2023_seur,samples_2022_2023_seur])
geodist_matr_ns_seur <- as.dist(geo_dist_seur_ns[samples_2022_2023_seur,samples_2022_2023_seur])

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
