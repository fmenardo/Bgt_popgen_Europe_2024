# Redundancy Analysis

We used RDA explore what factors shape the distribution of genetic variation in powdery mildew in Europe. We tested for the effects of local climatic conditions, wind connectivity and host ploidy.

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
We downloaded climate data from [CHELSA](https://chelsa-climate.org/). We used the average monthly climate data CHELSA V2.1 from climatologies 1981-2010. The 19 bioclim variables were complemented with 16 additional variables that were chosen based on biological relevance, resulting in a total of 35 climate variables, as listed [here](). For each variable we downloaded the corresponding tif-file. Each file consists of numeric data for every pixel. Information for every sample site was extracted with the coordinates and stacked using the extract and stack function respectively from the raster package in R. To avoid overfitting due to collinearity / multicollinearity amongst the climate variables, we excluded variables that had absolute pairwise correlation value >= 0.85 with another variable, finally retaining 12 variables. The script to perform these steps with the downloaded climate data is `data_prep.R`
