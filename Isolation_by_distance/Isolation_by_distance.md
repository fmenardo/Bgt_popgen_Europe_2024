# Isolation by Distance
We tested for isolation by geographic distance and isolation by wind distance in our data. We performed these analyses at three levels:
1. the Europe+_recent dataset
2. the Europe+_2022_2023 datasets including only European samples, i.e. excluding samples from Turkey and Israel
3. within the N_EUR and S_EUR populations, as classified by fineSTRUCTURE level 4, and including only the samples collected in 2022-2023    

We used the pairwise genetic distance matrix as described [here](../distance_matrix/distance_matrix.md), and the pairwise wind-distance matrix as described [here](../windscape/windscape.md). We used the `rdist.earth` function in the R package `fields` to construct pairwise geographic distance matrices. 

We performed the mantel test for correlation between matrices, as implemented in the R package `adegent` with the function `mantel.randtest` using 999 repititions. This gave us correlations between each of the three distance matrices for each dataset.

The script `ibd_ibw_wind_geo_corr.R` was used for the above steps as well as for generating the correlation histograms found in this directory. Here are the scatter plots for [genetic distance-geographic distance](ibd_dens_compare.pdf), [genetic distance-wind distance](ibw_density_compare.pdf) and [wind distance-geographic distance](corr_wind_geo_dens_compare.pdf) for each of the four datasets used.


### R packages
```
R 4.2.1
adegenet 2.1.0
fields 15.2
tidyverse 1.3.2
MASS 7.3.60
```
