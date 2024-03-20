# Isolation by Distance
We tested for isolation by distance in the Europe+_recent and Europe+_2022_2023 [datasets](../Datasets/Datasets.md). We used the genetic distance matrix as described [here](../distance_matrix/distance_matrix.md), and used the `rdist.earth` function in the R package `fields` to construct geographic distance matrices. 

We performed the mantel test for correlation between matrices, as implemented in the R package `adegent` with the function `mantel.randtest` using 999 repititions.  

The script `isolation_by_distance_clean.R` was used for the above steps as well as for generating the p-value histograms (`ibd_hist_tritici_europe+_recent.pdf`, `ibd_hist_tritici_europe+_2022_2023.pdf`) and the scatter plots (`ibd_density_tritici_europe+_recent.png`, `ibd_density_tritici_europe+_2022_2023.png`).
