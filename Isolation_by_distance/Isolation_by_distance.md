# Isolation by Distance
We tested for isolation by geographic distance and isolation by wind distance in our data. We performed these analyses at three levels:
1. the Europe+_recent dataset
2. the Europe+_2022_2023 datasets including only European samples, i.e. excluding samples from Turkey and Israel
3. within the N_EUR and S_EUR2 populations, as classified by fineSTRUCTURE level 4, and including only the samples collected in 2022-2023    

We used the pairwise genetic distance matrix as described [here](../distance_matrix/distance_matrix.md), and the symmetric pairwise wind-distance matrix as described [here](../windscape/windscape.md). The geographic distance matrix `geo_dist_matrix_tritici_europe_recent.csv` was constructed using the `rdist.earth` function in the R package `fields` as follows:
```
library(fields)
# Read in sample list and metadata file
tritici_eur_recent <- readLines("tritici_recent_extended_europe_2022+2023+ncsu.args")
all_meta <- read.csv("2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv")

rownames(all_meta) <- all_meta$Sample.Name
# select long, lat
meta_eur_recent <- all_meta[tritici_eur_recent,c(9,8)]

# compute distance matrix (in km)
dist_test <- rdist.earth(meta_eur_recent[,1:2], miles = FALSE, R = NULL)
rownames(dist_test)<- rownames(meta_eur_recent)
colnames(dist_test) <- rownames(meta_eur_recent)
geodist_test <- as.dist(dist_kms_eur)
write.csv(dist_test, "geo_dist_matrix_tritici_europe_recent.csv")
```

We performed the mantel test for correlation between matrices, as implemented in the R package `adegent` with the function `mantel.randtest` using 999 repititions. This gave us correlations between each of the three distance matrices for each dataset.

The script `isolation_by_distance_and_wind_new.R` was used for the above steps as well as for generating the correlation histograms found in this directory. Here are the scatter plots for [genetic distance-geographic distance](ibd_dens_compare_new.pdf), [genetic distance-wind distance](ibw_dens_compare_new.pdf) and [wind distance-geographic distance](corr_wind_geo_dens_compare_new.pdf) for each of the four datasets used.


### R packages
```
R 4.2.1
adegenet 2.1.0
fields 15.2
tidyverse 1.3.2
MASS 7.3.60
```
