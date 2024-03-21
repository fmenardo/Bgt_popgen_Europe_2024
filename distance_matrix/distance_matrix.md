# Distance matrix

A distance matrix based on pairwise differences between individuals (number of loci for which they differ) was computed using the `dist.gene` function in the R package `ape` using the method `pairwise`. The matrix was computed parallely for each of the 11 chromosomes using the script `dist_matrix_parallel.R` that took the per-chromosome SNP VCFs produced in [WGS_pipeline](../WGS_pipeline/WGS_pipeline.md) as input. The script was called using `input_dist_matrix`. An example from the input file is shown below:
```
Rscript dist_matrix_parallel.R -c 2 -i 2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_LR026984.1_chr1.vcf.gz -o chr1
```
The distances were then summed over all chromosomes to obtain a genome-wide distance matrix `gw_dist_mat_2022+before2022+2023+ncsu.csv`.

Clonal isolates were identified based on the values of the genome-wide distance matrix using `get_clones_from_dist_matr.R`. Samples that had a genetic distance of less than 120 SNPs were classified as clonal and a list of clonal groups was outputted `2022+before2022+2023+ncsu_all_list_clones_merged.txt`. 

![distr_gw_dist_vline120-1](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/c7a9f6dc-2162-46eb-98f4-68c4232ffb5e)

