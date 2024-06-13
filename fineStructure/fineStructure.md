# fineStructure
## Data preparation

For the fineStructure analyses we used the [Europe+](../Datasets/Datasets.md) datatset with 415 individuals sampled from Europe, the Middle East and Caucasus.
We kept only SNPS without any missing data, as fineStructure cannot handle them: 

```
gatk-4.4.0.0/gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe_large_tritici_no_clones_no_miss_new.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals MT880591.1 \                    #exclude mitochondrion
     --exclude-intervals LR026995.1_Un \                 #exclude conting with sequenced non assined to any chromosome
     --exclude-intervals Bgt_MAT_1_1_3 \                 #exclude contig with alternative mating type
     --sample-name tritici_extended_europe_2022+before2022+2023+ncsu.args      
```

This resulted in 1’201’198  SNPs.

To generate the input files for fineStructure we need to know the per base recombination rates, these where [obtained](../recombination_map/recombination_map.md) starting from a genetic map and are stored in the file `../recombination_map/THUN12x96224_bp_recombination_rates.txt`

We generate the id and phase files for fineStructure, we also generate a .pos file:

```
zcat Europe_large_tritici_no_clones_no_miss.vcf.gz > Europe_large_tritici_no_clones_no_miss.vcf
python make_input_files_4_fs.py -vcf Europe_large_tritici_no_clones_no_miss.vcf -o Europe_large
```

Starting from the .pos file and the recombination map we generate the recombination file for fineStructure

```
python make_input_rec_file_4_fs.py -rec THUN12x96224_bp_recombination_rates.txt -i Europe_large.pos_file -o Europe_large
```


As fineStructure cannot deal with sample names starting with digits we rename these two isolates:
```
sed -i 's/96224/CHE_96224/g' Europe_large.id_file
sed -i 's/94202/CHE_94202/g' Europe_large.id_file
```

## fineStructure

We run fineStructure with default parameters, except that we increase the number of iterations in the EM algorithm to 50 (default 10)

```
fs Europe_large -idfile Europe_large.id_file -phasefiles Europe_large.hap_file -recombfiles Europe_large_cp_rec_file.txt -ploidy 1 -v -n -hpc 1 -s1args:-in\ -iM\ -i\ 50\ --emfilesonly -go
```
With the command above fineStructure generates lists of commands to run at different stages, these list can be submitted as batch jobs to a computer cluster.

When all stages are completed the three most important outputs are the chromopainter chunkcounts file (`Europe_large_linked_hap.chunkcounts.out`), the fineStructure mcmc file (`Europe_large_linked_hap_mcmc.xml`), and the fineStructure tree file (`Europe_large_linked_hap_tree.xml`).

## Plot results and PCA

This code is based on the example provided by authors of fineStructure. You need FinestructureLibrary.R and FinestructureDendrogram.R provided [here](https://people.maths.bris.ac.uk/~madjl/finestructure/toolsummary.html).

```
Rscripts plot_fs_results.R
```

This generates the [coacestry matrix plot](Europe_large_Coancestry.pdf) and the [PCA plot](Europe_large_PCA.pdf)

## Software versions
```
gatk 4.4.0.0
fineStructure 4.4.4

- python and python modules

python 3.10.9
numpy 1.23.5    
argparser 1.4.0

- R and packages

sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RColorBrewer_1.1-3    Polychrome_1.5.1      patchwork_1.1.2       ggpubr_0.6.0          pals_1.8              ggplot2_3.4.2         psych_2.3.12         
 [8] GPArotation_2023.11-1 paran_1.5.2           MASS_7.3-60           XML_3.99-0.14         ape_5.7-1            

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0     purrr_1.0.1          lattice_0.21-8       carData_3.0-5        colorspace_2.1-0     vctrs_0.6.3          generics_0.1.3       utf8_1.2.3          
 [9] rlang_1.1.3          pillar_1.9.0         glue_1.6.2           withr_2.5.0          lifecycle_1.0.3      munsell_0.5.0        ggsignif_0.6.4       gtable_0.3.3        
[17] ragg_1.2.5           mapproj_1.2.11       labeling_0.4.2       parallel_4.2.1       fansi_1.0.4          broom_1.0.5          Rcpp_1.0.12          scales_1.3.0        
[25] backports_1.4.1      scatterplot3d_0.3-44 abind_1.4-5          systemfonts_1.0.4    farver_2.1.1         textshaping_0.3.6    mnormt_2.1.1         digest_0.6.33       
[33] rstatix_0.7.2        dplyr_1.1.2          grid_4.2.1           cowplot_1.1.1        cli_3.6.1            tools_4.2.1          magrittr_2.0.3       maps_3.4.1          
[41] tibble_3.2.1         dichromat_2.0-0.1    crayon_1.5.2         tidyr_1.3.0          car_3.1-2            pkgconfig_2.0.3      rstudioapi_0.15.0    R6_2.5.1            
[49] nlme_3.1-157         compiler_4.2.1    

```
