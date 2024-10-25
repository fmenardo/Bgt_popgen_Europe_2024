# Principal component analysis (PCA) 
We performed PCAs for two of our [datasets](../Datasets/Datasets.md): World and Europe+. Both datasets were filtered to include only biallelic SNPs, no singletons and maximum 10% missing data using GATK `SelectVariants`.
```
#### world dataset ####
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \   ### VCF file with all B.g.tritici isolates, all biallelic SNPs
    -XL MT880591.1 \     ### exclude mitochrondrion
    -XL Bgt_MAT_1_1_3 \  ### exclude alternate mating type locus
    --select "AC>1 && AC<567" \  ### remove singletons
    --max-nocall-fraction 0.1 \  ### max missing data
    -O 2022+before2022+2023+ncsu_11chr_no_sing_maxmiss0.1.vcf.gz

#### europe+ dataset ####
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \   ### VCF file with all B.g.tritici isolates, all biallelic SNPs
    --sample-name tritici_extended_europe_2022+before2022+2023+ncsu.args \  ### Europe+ dataset
    -XL MT880591.1 \     ### exclude mitochrondrion
    -XL Bgt_MAT_1_1_3 \  ### exclude alternate mating type locus
    --select "AC>1 && AC<415" \  ### remove singletons
    --max-nocall-fraction 0.1 \  ### max missing data
    -O tritici_extended_eur_2022+before2022+ncsu_nosing_maxmiss0.1.vcf.gz
```

The script to run the PCA is `pca_clean.R`. It takes as input the filtered VCF files, a [metadata file](../Datasets/2022+before2022+2023+ncsu_metadata+fs+admxK7_19032024.csv) and the number of cores to be used, as shown below:
```
Rscript pca_clean.R -s -c 6 -l 2022+before2022+2023+ncsu_metadata+fs+admxK7_19032024.csv -C 10 -i 2022+before2022+2023+ncsu_11chr_no_sing_maxmiss0.1.vcf.gz -o tritici_world
```

The PC scores and eigenvalues for all samples are written to .csv files, which can then be used for plotting using the script `plot_pca_clean.R`. The plots for the PCA are [here](fig_s2_world_pca.pdf) for the _World_ dataset and [here](fig_s5_pca_eur+.pdf) for the _Europe+_ dataset.

### R packages 
```
R 4.2.1
vcfR 1.14.0
adegenet 2.1.0
argparser 0.7.1
ggplot2 3.5.0
patchwork 1.2.0
pals 1.7
```
