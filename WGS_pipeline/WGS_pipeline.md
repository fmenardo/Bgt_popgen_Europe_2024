# WGS pipeline

The pipeline starts with raw WGS (short) reads and returns a VCF file with all biallelic SNPs for all *B.g. tritici* isolates in the [World](../Datasets/Datasets.md) dataset (CHECK AT THE END, this will depend on the vcf that we share). 

### Software used
1. Python3
2. [fastp 0.23.2](https://github.com/OpenGene/fastp) 
3. [bwa 0.7.17-r1198-dirty](https://github.com/lh3/bwa)
4. [htslib 1.17](https://github.com/samtools/htslib/releases/tag/1.17)
5. [samtools 1.17](https://github.com/samtools/samtools/releases/tag/1.17)
6. [bcftools 1.17](https://github.com/samtools/bcftools/releases/tag/1.17)
7. [GATK 4.4.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.4.0.0)
8. R

### Workflow

1. `pipeline_with_gatk_statscsv.py` 
This takes as input the path to the raw fastq-files and reference genome along with some quality-based trimming and adapter trimming parameters. It returns a per-chromosome VCF file called by GATK `HaplotypeCaller` and some summary statistics about the mapping. It was run as an array job for 737 samples (including all *formae speciales*) (ADD FILE WITH LIST OF ALL ACCESSIONS) using the same input parameters for all samples. For example:
```
python3 ../scripts/pipeline_with_gatk_statscsv.py -ref GCA_900519115.1_2022_bgt_ref_mating_type.fa -minlen 50 -rw 5 -fw 1 -rq 20 -fq 20 -i /home/jjigis/projects/bgt_sequence_data/2023_collection/CHNY072301_R1.fastq.gz
```
2. Samples with coverage less than 15x were identified [(n = 26)](coverage_below_15) and excluded from all subsequent analyses. 
3. The VCF files for all remaining samples (n = 711) were combined (per-chromosome) using GATK `CombineGVCFs`. 
4. Variants were called on the combined VCF files for all chromosomes using GATK `GenotypeGVCFs`.
```
gatk --java-options "-Xmx20g" GenotypeGVCFs \
-R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
-V 2022+before2022+2023+ncsu_all_combinegvcfs_$CHROMOSOME.vcf.gz \
-O 2022+before2022+2023+ncsu_all_combinegvcfs_genotyped_$CHROMOSOME.vcf.gz \
--include-non-variant-sites \
-A StrandBiasBySample 
```
5. In order to decide threshold values for filtering variants, the distribution of annotation values for SNPs were visualised. For each chromosome, SNPs were first selected from the output of step #4  using GATK `SelectVariants` `--select-type-to-include SNP` and their annotation values were written to a table using GTAK `VariantsToTable`. Histograms were plotted for the genome-wide (chromosomes 1-11) values of the annotations `QD`, `FS`, `SOR`, `MQ`, `MQRankSum` and `ReadPosRankSum` using R `ggplot2`.

![2022+before2022+2023+ncsu_WG_gatk_info_distr-1](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/8e636ad7-1f92-4808-8250-f6d72ebaeb85)

6. Hard filtering for sites was done using GATK `VariantFiltration`. The threshold values (informed by the distributions plotted in step #5 and [GATK's recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)) were `QD < 10.0`, `FS > 50.0`, `MQ < 45.0`, `ReadPosRankSum < -5.0 || ReadPosRankSum > 5.0`.
```
gatk --java-options "-Xmx4g" VariantFiltration \
   -R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
   -V 2022+before2022+2023+ncsu_all_combinegvcfs_genotyped_$CHROMOSOME.vcf.gz \
   -O 2022+before2022+2023+ncsu_all_combinegvcfs_genotyped_varfilt_$CHROMOSOME.vcf.gz \
   --filter-name "QualByDepth10" \
   --filter-expression "QD < 10.0" \
   --filter-name "FisherStrand50" \
   --filter-expression "FS > 50.0" \
   --filter-name "RMSMapQual45" \
   --filter-expression "MQ < 45.0" \
   --filter-name "ReadPosRankSumTest5" \
   --filter-expression "ReadPosRankSum < -5.0 || ReadPosRankSum > 5.0" 
```
7. `recode_multivcf_after_gatk_2024_new_clean.py` This script takes the genotyped VCF file produced by GATK in step #3 as input and recodes the 'GT' field value as '.' for all sites at which (a) depth of high-quality informative reads < user-defined minimum depth or (b) variant call is supported by < 90% of the high-quality informative reads. It outputs the recoded VCF (which can then be gzipped and indexed using `bgzip` and `tabix -p vcf`) and a csv file with the number of positions failing each of the above filters as well as the total number or missing positions and the total number of variant positions for each sample. The genome-wide distribution (obtained after summing over chromosomes 1-11) was plotted using R `ggplot2`.

![2022+before2022+2023+ncsu_recoding_stats-1](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/59844197-a2c1-46e0-93e5-da85b9386ce9)

The samples with > 200,000 'heterozygous positions', i.e. positions at which the variant support was < 90%, were excluded from all further analyses [(n=13)](200k_het_pos_exclude_dact.args)  
8. SNPs were selected from the chromosomal VCF files using GATK `SelectVariants` with options `--select-type-to-include SNP`, `--restrict-alleles-to ALL` and sites failing the filters from step #6 were excluded using the option `--exclude-filtered`. The resulting VCFs contained some [spanning deletions](https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele) denoted by '\*' . As these would have caused problems in downstream analyses, all sites with an '*' were removed using the script `recode_asterisk_count_snp.py` which returned a modified VCF file that was gzipped and indexed using `bgzip` and `tabix -p vcf` respectively. The VCF files (with filtered variants and no asterisks) for chromosomes 1-11, the alternate mating type locus and the mitochondrion were merged using the `concat` option in bcftools. The Bgt_Un "chromosome" with contigs not assigned to any other chromosomes was excluded from all further analyses.  
9. The resulting VCF files were then subset to include only biallelic SNPs and only *B.g. tritici* isolates. This step was performed using GATK `SelectVariants` and samples in `2022+before2022+2023+ncsu_tritici_list.args` were included while samples with >200k het pos `2022+before2022+2023+ncsu_200k_hetpos_to_exclude_list.args` and clones (as decided based on the [dist matrix analysis](../distance_matrix/distance_matrix.md) ) `2022+before2022+2023+ncsu_tritici_clones_to_exclude_list.args` were excluded. The final list of samples in this resulting VCF file made up the [World](../Datasets/Datasets.md) dataset (CHEK at the end, this will change depending on what we share). 
```
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V ../project_data_prep/data/2022+before2022+2023+ncsu_recoded_snps_filtered_no_asterisk_11chr_mt_MAT.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    --sample-name ../project_data_prep/data/2022+before2022+2023+ncsu_tritici_list.args \
    --exclude-sample-name ../project_data_prep/data/2022+before2022+2023+ncsu_tritici_clones_to_exclude_list.args \
    --exclude-sample-name ../project_data_prep/data/2022+before2022+2023+ncsu_200k_hetpos_to_exclude_list.args \
    --select "AC>0 && AC<568" \
    -O tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz
```

