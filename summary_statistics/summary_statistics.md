# Summary Statistics

We computed measures of intra and inter-population genetic diversity for the [Europe+_recent](../Datasets/Datasets.md) dataset. 
We used the 'level 4' population classification, as inferred from fineSTRUCTURE, which divided the dataset into five populations, namely ME, N_EUR, S_EUR1, S_EUR2 and TUR.

### Within-population statistics
We calculated pi, Watterson's theta and tajima's D for each population in windows of 10kb covering all 11 chromosomes. 

1. Pi  
We used [pixy](https://pixy.readthedocs.io/en/latest/) for calculating pi. First, we generated per-population VCF files containing both variant and invariant sites ("allsites VCF") and excluded sites that had been flagged by [previous filters](../WGS_pipeline/WGS_pipeline.md)  or those that had > 50% missing data using GATK `SelectVariants`. We handled the 11 chromosomes separately and parallely to increase computation speed.
```bash
gatk SelectVariants \
  -R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
  -V 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf.gz \
  --exclude-filtered \      # exclude sites that had failed previous site-level filters
  --sample-name tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_ME.args \   # include only samples in the desired population
  --max-nocall-fraction 0.5 \   # filter out sites with > 50% missing data
  -O tritici_fs4_ME_all_sites_filtered_maxmiss_0.5_$CHROMOSOME.vcf.gz
```
Next, we ran pixy to calculate pi in windows of 10 kb for each population separately. 
```bash
pixy --stats pi \
 --vcf tritici_fs4_ME_all_sites_filtered_maxmiss_0.5_$CHROMOSOME.vcf.gz \
 --populations tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_ME.pops \
 --window_size 10000 \
 --n_cores 2 \
 --output_folder output_fs_level4_10kb_per_pop \
 --output_prefix pixy_tritici_ext_eur_recent_fs_level4_maxmiss0.5_ME_$CHROMOSOME
```
This returned .txt files with average per-site pi values and the actual number of sites with valid genotypes in each window. These two values were multiplied to obtain a 'pi_in_window' for each window. 

2. Watterson's theta    
We first created per-population VCF files containing only biallelic SNPs and sites with less than 50% missing data.
```bash
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    --max-nocall-fraction 0.5 \
    --sample-name tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_ME.args \
    --select "AC>0 && AC<42" \
    -O tritici_fs4_ME_maxmiss0.5_biallelic_snps.vcf.gz
```
Next, we used VCFtools `SNPdensity` to get the number of SNPs in each 10 kb window across the genome.
```bash
vcftools --gzvcf tritici_fs4_ME_maxmiss0.5_biallelic_snps.vcf.gz --out tritici_fs4_ME_maxmiss0.5 --temp ~/scratch/ --SNPdensity 10000
```
The raw number of SNPs in each window was divided by the number of individuals and the number of sites in each window to obtain a normalised measure of Watterson's theta.

3. Tajimas's D     
We used the above-mentioned estimates of pi and Watterson's theta to calculate window-wise Tajima's D for all 5 populations using the script `calculate_tajimas_d_maxmiss.R`. The calculations were based on Tajima (1989), as also described [here](tajimasD_calculation.pdf).

The box-plots of the three within-population statistics, calculated in 10 kb windows over the entire genome, can be found [here](fs4_pi_theta_tajimasD_maxmiss0.5_boxplot_no_outliers.pdf). `fs4_pi_theta_tajimasD_maxmiss0.5.csv` is the table with compiled statistics for all populations.

Pi, theta and Tajima's D were also visualised as Manhattan (MH) plots. The input file (`fs4_pi_theta_tajimasD_maxmiss0.5.csv`) was split by population using `separate_pops_in_new_files.R` and the corresponding manhattan plots were generated using `MH_pi.R`, `MH_TajimasD.R` and `MH_wattersons_theta.R`. 

### Between-population statistics
Measures of genetic variation between populations namely, dxy and Fst were computed using pixy in 10 kb windows spanning 11 chromosomes for each population pair. 
Chromosomal all-site VCF files containing all samples and sites with less than 50% missing data were created using GATK `SelectVariants`. The 11 chromosomes were handled independently and parallely. 
```bash
gatk SelectVariants \
  -R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
  -V 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf.gz \
  --exclude-filtered \    # exclude sites that had failed previous site-level filters
  --sample-name tritici_recent_extended_europe_2022+2023+ncsu.args \   # Europe+_recent dataset
  --max-nocall-fraction 0.5 \  # filter out sites with > 50% missing data
  -O tritici_ext_eur_recent_all_sites_filtered_maxmiss_0.5_$CHROMOSOME.vcf.gz

```
These VCF files were used as input for pixy to calculate Weir and Cockerham’s Fst and dxy (average nucleotide difference between populations).
```bash
pixy --stats fst dxy \
 --vcf tritici_ext_eur_recent_all_sites_filtered_maxmiss_0.5_$CHROMOSOME.vcf.gz \
 --populations tritici_ext_eur_recent_fs_level4_populations \
 --window_size 10000 \
 --n_cores 2 \
 --output_folder output_fs_level4_10kb \
 --output_prefix pixy_tritici_ext_eur_recent_fs_level4_$CHROMOSOME
```
The outputs of both statistics and all chromosomes were combined and [plotted](fs4_dxy_wc-fst_boxplot_no_outliers_grey.pdf) using `plot_pairwise_comparisons.R`
