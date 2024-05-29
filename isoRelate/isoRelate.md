# isoRelate

We used [isoRelate](https://github.com/bahlolab/isoRelate) to perfome genome scans of selection. This method can identify loci that have been under positive selection in the ecent past (roughly within 20 sexual generations).

## Data preparation

For the isoRelate analysis we used the [Europe+_recent](../Datasets/Datasets.md) datatset with 368 individuals sampled from Europe, the Middle East and Caucasus not earlier than 2015.
We ran the analyses separately on the 5 populations (ME, N_EUR, S_EUR+, S_EUR1, TUR) defined by the fineStructure analysis (level 4). LINK to metainfo.

As an example, here we report the code to perform the analysis for one population (N_EUR).

First we identify intervals in which markers cannot be mapped unambiguosly on the genetic map. This needs to be done only once, and not for each population. 

```
python find_ambiguity_in_rec_map.py -rec ../recombination_map/THUN12x96224_genetic_map_in_cM_+_phy_distance
```

Then we selected kept only SNPs without any missing data. We also excluded positions that cannot be mapped unambiguosly on the genetic map, and SNPs with a minor allele frequency less than 5%. 

```
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe+_recent_tritici_no_clones_no_miss_no_amb_cM_NE.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals THUN12x96224_genetic_map_in_cM_+_phy_distance.ambiguous_intervals.list \
     --exclude-intervals MT880591.1 \
     --exclude-intervals LR026995.1_Un \
     --exclude-intervals Bgt_MAT_1_1_3 \
     --select "AF > 0.05" \
     --sample-name tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_N_EUR.args

```
We generated the ped and map file needed by isoRelate as follow:
```
plink --allow-extra-chr --vcf Europe+_recent_tritici_no_clones_no_miss_no_amb_cM_NE.vcf.gz --recode 12 --double-id --out BgtE+r_N_EUR --threads 1

awk '$5="1" && $6="0"' BgtE+r_N_EUR.ped >  BgtE+r_N_EUR_mod.ped
sed -E 's/\S+_chr//g' BgtE+r_N_EUR.map > BgtE+r_N_EUR_mod.map


python add_cM_to_map.py -map BgtE+r_N_EUR_mod.map -rec ../recombination_map/THUN12x96224_genetic_map_in_cM_+_phy_distance -o BgtE+r_N_EUR_mod

```
## isoRelate
We ran isoRelate in two steps for better compuational efficiency, the first step infers identity by descent segments, while the second step computes the proportion of IBD pairs for each SNP and its significance level.


```
Rscript run_ibd_step1.R -o BgtE+r_N_EUR_2cM -p BgtE+r_N_EUR_mod.ped -m BgtE+r_N_EUR_mod_cM.map -c 10 -C 2
Rscript run_ibd_step2.R -o BgtE+r_N_EUR_2cM

```
These scripts produce a set of output files, most importantly the file `BgtE+r_N_EUR_2cM_iR_table.txt` contains the p-value for each SNP.

## Plot results
We plotted the p-values in Manhattan plots, based on the table above. One can choose to zoom in a specific region using the second script. The plots for the 5 populations are available in this folder.

```
MH_isoRelate_noCap_plots_loop.R
MH_isoRelate_noCap_zoom_in_region.R
```

R packages used for ploting:

```
ggplot2 v3.4.4
cowplot v1.1.3
dplyr v1.1.2
stringr v1.5.0
```

## Software versions
```
gatk v4.4.0.0
plink v1.90b6.21

- python and python modules

python v3.11.4

- R and packages

sessionInfo()
R version 4.3.1 (2023-06-16)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8    LC_NUMERIC=C            LC_TIME=C              
 [4] LC_COLLATE=en_US.UTF-8  LC_MONETARY=C           LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=C              LC_NAME=C               LC_ADDRESS=C           
[10] LC_TELEPHONE=C          LC_MEASUREMENT=C        LC_IDENTIFICATION=C    

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] argparser_0.7.1   doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2    
[5] isoRelate_0.1.0  

loaded via a namespace (and not attached):
 [1] vctrs_0.6.3      cli_3.6.1        rlang_1.1.1      generics_0.1.3  
 [5] ggnetwork_0.5.12 glue_1.6.2       colorspace_2.1-0 scales_1.2.1    
 [9] fansi_1.0.4      grid_4.3.1       munsell_0.5.0    tibble_3.2.1    
[13] lifecycle_1.0.3  compiler_4.3.1   dplyr_1.1.2      codetools_0.2-19
[17] igraph_1.5.1     Rcpp_1.0.11      pkgconfig_2.0.3  R6_2.5.1        
[21] tidyselect_1.2.0 utf8_1.2.3       pillar_1.9.0     magrittr_2.0.3  
[25] gtable_0.3.4     ggplot2_3.4.3  
```
