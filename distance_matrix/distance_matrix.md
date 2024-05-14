# Distance matrix

A distance matrix based on pairwise differences between all 711 individuals (number of loci for which they differ) was computed using the `dist.gene` function in the R package `ape` using the method `pairwise`. Since there was great variation in the number of missing data and variants called among samples, we normalised this distance by the actual number of loci compared for each pair. This was done by setting `pairwise.deletion = TRUE` and recording the variance for each pair (`variance=TRUE`). The number of loci was then calculated using the formula D^2/(D-V), where D = distance and V = variance, following the [official manual](https://search.r-project.org/CRAN/refmans/ape/html/dist.gene.html) of the tool. 

To speed up computation, we divided the genome (chromosomes 1-11) into 143 windows of ~ 1 million bases. The per-chromosome 'all-sites' VCF files generated in the [WGS pipeline](../WGS_pipeline/WGS_pipeline.md) were divided into windows using `bcftoools`, for example:
```
bcftools view -r LR026984.1_chr1:1-1000000 2022+before2022+2023+ncsu_covg15_recoded_LR026984.1_chr1.vcf.gz -Oz -o 2022+before2022+2023+ncsu_covg15_recoded_LR026984.1_chr1_1-1000000.vcf.gz && tabix -p vcf 2022+before2022+2023+ncsu_covg15_recoded_LR026984.1_chr1_1-1000000.vcf.gz

```

The distance matrix was then computed parallely for each of these windows using the script `dist_matrix_parallel_loci.R` that took the corresponding VCF file as input. The script was called using `input_dist_matrix_windows`. An example from the input file is shown below:
```
Rscript dist_matrix_parallel_loci.R -c 2 -i 2022+before2022+2023+ncsu_covg15_recoded_LR026984.1_chr1_1-1000000.vcf.gz -o chr1_1-1000000
```  
The script produced a 'difference' matrix (number of loci that were different) and a 'loci' matrix (total number of loci compared) for each window. The differences and loci were summed over all windows and a genome-wide distance matrix was obtained by dividing the total differences by the total number of loci `gw_dist_mat_prop_2022+before2022+2023+ncsu.csv`. 
Clonal isolates were identified based on the values of the genome-wide distance matrix. Samples that had a genetic distance of less than 9e-05 were classified as clonal and a list of clonal groups was outputted `2022+before2022+2023+ncsu_all_list_clones_merged.txt`. 

The above steps were performed using the script `dist_matrix_WG_from_windows.R`.
![2022+before2022+2023+ncsu_dist_prop_distr](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/ab3b4cbe-5810-4835-96d6-dd575f847017)

Overall, we identified 35 clonal groups of _B.g. tritici_. 33 of these groups contained 2-3 isolates each that were collected from the same pot or field. We retained only one isolate from each such group, eliminating 39 isolates. The remaining two clonal groups contained isolates that were collected from distant locations but handled together in the laboratory. These groups, containing 9 isolates overall, were excluded completely as they were suspected to be contaminated.   

### R packages
```
R 4.2.3
vcfR 1.14.0
adegenet 2.1.10
ape 5.7.1
argparser 0.7.1
tidyverse 1.3.2
```
