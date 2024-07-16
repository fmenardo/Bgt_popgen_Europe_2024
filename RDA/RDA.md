# Redundancy Analysis

We used RDA explore what factors shape the distribution of genetic variation in powdery mildew in Europe. We tested for the effects of local climatic conditions and host ploidy.

### 1. Climate          
We used the _Europe+_recent_ dataset to see how climate variables affected genetic diversity. We used only biallelic SNPs, removed all missing data and filtered for minor allele frequency 0.05 using GATK.
```bash
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe+_recent_tritici_no_miss_maf_0.05.vcf.gz \
     --select "AF>0.05 && AF<0.95" \
     --max-nocall-fraction 0 \
     --exclude-intervals MT880591.1 \
     --exclude-intervals LR026995.1_Un \
     --exclude-intervals Bgt_MAT_1_1_3 \
     --sample-name tritici_recent_extended_europe_2022+2023+ncsu.args
```
This VCF file was then converted to a binary genotype matrix using the R package vcfR.      
```R
data <- read.vcfR("Europe+_recent_tritici_no_miss_maf_0.05.vcf.gz")
genotypes <- extract.gt(data)
genotypes <- as.data.frame(t(genotypes)) #transpose matrix
genotypes <- mutate_all(genotypes, as.factor)
write.csv(genotypes, file="genotypes.csv")
```                
We downloaded climate data from [CHELSA](https://chelsa-climate.org/). We used the average monthly climate data CHELSA V2.1 from climatologies 1981-2010. The 19 bioclim variables were complemented with 16 additional variables that were chosen based on biological relevance, resulting in a total of 35 climate variables, as listed [here](clim_variables_list.csv). For each variable we downloaded the corresponding tif-file. Each file consists of numeric data for every pixel. Information for every sample site was extracted with the coordinates and stacked using the extract and stack function respectively from the raster package in R. To avoid overfitting due to collinearity / multicollinearity amongst the climate variables, we excluded variables that had absolute pairwise correlation value >= 0.85 with another variable, finally retaining 12 variables. The script to perform these steps with the downloaded climate data is `data_prep.R`

For the RDA, we included the first three wind coordinates (as obtained from the [windscape](../windscape/windscape.md) analysis), sampling coordinates and country as covariates. The values of these variables as well as the 12 climate variables is given in [this table](Variables_without_climcorrelation.csv). Full and partial RDAs, including the step of forward variable selection, were performed using the script `rda.R` which was ran using 
```bash 
Rscript rda.R -v Variables_without_climcorrelation.csv -g genotypes.csv -o "output/"
```
This was followed by ANOVA and variance partitioning for each of the models using the scripts `anova_new.R` and `variance_partitioning.R` that were run as
```bash
Rscript anova_new.R -i "output/" -o "output/"
Rscript variance_partitioning.R -i "output/" -o "output/"
```

### 2. Host    
We investigated if the kind of host Bgt was sampled from (hexaploid bread wheat or tetraploid durum wheat) could explain some of the patterns of genetic variation. We used [131](2022_2023_field_ploidy_list.args) samples from the Europe+_2022_2023 dataset that were all sampled from fields and whose host information was available. We subset the VCF file to only include these samples and also filtered for zero missing data and minor allele frequency 0.05 using GATK.  
```bash
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O tritici_2022_2023_host_rda_no_miss_maf_0.05_2.vcf.gz \
     --select "AF>0.05 && AF<0.95" \
     --max-nocall-fraction 0 \
     --exclude-intervals MT880591.1 \
     --exclude-intervals LR026995.1_Un \
     --exclude-intervals Bgt_MAT_1_1_3 \
     --sample-name 2022_2023_list_host_rda.args
```
This VCF file was converted to a binary genotype matrix using vcfR as described in #1. For the RDA, we included the 12 climate variables mentioned above, as well as the 3 wind coordinates, smapling location and country as covariates. Full and partial RDA models for each of the variables were run using `rda_host.R`, followed by ANOVA and variance partitioning using `anova_new_host.R` and `variance_partitioning_host.R` as:
```bash
Rscript rda_host.R -v Variables_host_without_climcorr.csv -g genotypes.csv -o "output/"
Rscript anova_new_host.R -i "output/" -o "output/"
# can run variance_paritioning_host.R interactively.
```
