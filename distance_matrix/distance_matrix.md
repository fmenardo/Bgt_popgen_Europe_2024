# Distance matrix

A distance matrix based on pairwise differences between individuals (number of loci for which they differ) was computed using the `dist.gene` function in the R package `ape` using the method `pairwise`. The matrix was computed parallely for each of the 11 chromosomes using the script `dist_matrix_parallel.R` and `input_dist_matrix` as input.  The distances were then summed over all chromosomes to obtain a genome-wide distance matrix `gw_dist_mat_2022+before2022+2023+ncsu.csv`.

Clonal isolates were identified based on the values of the genome-wide distance matrix using `get_clones_from_dist_matr.R`. Samples that had a genetic distance of less than 120 SNPs were classified as clonal and a list of clonal groups was outputted `2022+before2022+2023+ncsu_all_list_clones_merged.txt`. 
