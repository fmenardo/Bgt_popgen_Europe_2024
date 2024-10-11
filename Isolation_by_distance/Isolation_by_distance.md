# Isolation by Distance
We tested for isolation by geographic, wind and climatic distances in our data. We performed these analyses for three datasets:
1. the _Europe+_2022_2023_ dataset
2. all samples in _Europe+_2022_2023_ that belonged to the N_EUR population (fineSTRUCTURE level 4 classification)
3. all samples in _Europe+_2022_2023_ that belonged to the S_EUR2 population (fineSTRUCTURE level 4 classification) 

The list of samples belonging to each of these datasets can be found in the *.args files in this directory.

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
write.csv(dist_test, "geo_dist_matrix_tritici_europe_recent.csv")
```
For climatic distances, we used the 12 bioclim variables described [here](../RDA/RDA.md). We performed a PCA of these 12 variables using the prcomp function in R and computed the euclidean distance between all pairs of samples based on the first 7 principal components using the ‘dist’ function in R to obtain a pairwise climatic distance matrix:

```R
md <- read.csv("~/projects/project_bgt_popgen/analysis/rda/clim/Variables_without_climcorrelation.csv", row.name = 1)  ## read in clim data values for all isolates
my_pca <- prcomp(md[,16:27], scale = TRUE)  ## perform PCA on 12 clim variables
pca_coords <- as.data.frame(my_pca$x[,1:7]) ## retain first 7 PCs
clim_dist_unw <- dist(pca_coords, method = "euclidean")  ## calculate pairwise euclidean distances 
clim_dist <- as.matrix(clim_dist_unw)
write.csv(clim_dist, "clim_dist_unweighted_7pc.csv")
```

We performed the mantel test for correlation between matrices, as implemented in the R package `adegent` with the function `mantel.randtest` using 999 repititions. This gave us correlations between each of the four distance matrices for each dataset. 

Additionally, we also tested for isolation by geographic distance along the east-west and nort-south axes separately for the _N_EUR_2022_2023_ and _S_EUR2_2022_2023_ datasets. The georaphic distance matrices were computed using the script `geo_dist_axes.R` and the genetic distance matrix was the same as described above. 

The script to perform all mantel tests and generate corresponding density plots is `all_mantel_tests.R`. The test results are reported in [this table](mantel_test_results.csv).

## Logistic Regression
We tested how well geographic, wind and climatic distances could predict population structure using logistic regression. Using samples belonging to the N_EUR_2022_2023 and S_EUR2_2022_2023 datasets, we modelled which factors could predict whether two individuals belonged to the same or different populations. For each pair of individuals, the response variable (Diff pop) was 0 if they belonged to the same population and 1 if different. The distance matrices were the same as described above. Logistic regression was performed using the rms package in R using the script `logsitic_regression.R`.


### R packages
```
R 4.2.1
adegenet 2.1.0
fields 15.2
tidyverse 1.3.2
MASS 7.3.60
vegan 2.6.4
permute 0.9.7
lattice 0.21.8
```
