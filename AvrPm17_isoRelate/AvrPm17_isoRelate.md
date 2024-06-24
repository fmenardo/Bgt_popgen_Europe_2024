# Identity by descent clusters analysis for AvrPm17

In the genome-wide analysis with isoRelate we found that in many populations there is an excess of IBD pairs at the locus containining AvrPm17. Here we explore which pairs of isolates are in IBD over this locus.

First we run again isoRelate, only this time we focus on the locus of AvrPm17 and we use the compete [Europe+ dataset](../Datasets/Datasets.md) (we do not run separate analysis for different populations as we did for the genome-wide analysis)

We prepared the input files:

```
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe+_recent_chr1_avrpm17.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals THUN12x96224_genetic_map_in_cM_+_phy_distance.ambiguous_intervals.list \
     --intervals LR026984.1_chr1:3000000-6000000 \
     --select "AF > 0.05" \
     --sample-name ../Datasets/tritici_recent_extended_europe_2022+2023+ncsu.args


plink --allow-extra-chr --vcf Europe+_recent_chr1_avrpm17.vcf.gz --recode 12 --double-id --out BgtE+r_chr1_avrpm17 --threads 1

awk '$5="1" && $6="0"' BgtE+r_chr1_avrpm17.ped >  BgtE+r_chr1_avrpm17_mod.ped
sed -E 's/\S+_chr//g' BgtE+r_chr1_avrpm17.map > BgtE+r_chr1_avrpm17_mod.map


python ../isoRelate/add_cM_to_map.py -map BgtE+r_chr1_avrpm17_mod.map -rec ../recombination_map/THUN12x96224_genetic_map_in_cM_+_phy_distance -o BgtE+r_chr1_avrpm17_mod
```
We ran isoRelate

```
Rscript ../isoRelate/run_ibd_step1.R -o BgtE+r_avrpm17 -p BgtE+r_chr1_avrpm17_mod.ped -m BgtE+r_chr1_avrpm17_mod_cM.map -c 5 -C 2
```

And finally we produced clusters and plots with `plot_cluster.R` in the following R environment:

```
sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C      

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] patchwork_1.2.0    RColorBrewer_1.1-3 igraph_1.5.0       dplyr_1.1.2        maps_3.4.1         ggplot2_3.4.2      isoRelate_0.1.0   

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.12       rstudioapi_0.15.0 magrittr_2.0.3    tidyselect_1.2.0  munsell_0.5.0     colorspace_2.1-0  R6_2.5.1          rlang_1.1.3      
 [9] foreach_1.5.2     fansi_1.0.4       tools_4.2.1       grid_4.2.1        gtable_0.3.3      utf8_1.2.3        cli_3.6.1         withr_2.5.0      
[17] iterators_1.0.14  tibble_3.2.1      lifecycle_1.0.3   vctrs_0.6.3       ggnetwork_0.5.12  codetools_0.2-19  glue_1.6.2        compiler_4.2.1   
[25] pillar_1.9.0      generics_0.1.3    scales_1.3.0      pkgconfig_2.0.3
```
