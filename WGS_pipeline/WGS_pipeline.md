# WGS pipeline

## Sampling and data

We sampled 276 isolates of *B.g. tritici* across Europe and the Mediterranean and whole-genome sequenced their DNA with short reads. The sequences can be found under the BioProject accession PRJEB75381. Additionally, we downloaded previously published, publicly available sequences of *B.g. tritici* (N = 375) using the `fasterq-dump` (v 3.0.5) method of [SRA tools](https://github.com/ncbi/sra-tools). For example,
```
fasterq-dump --split-files SRR11548116
```
The list of accessions of all downloaded sequences along with the metadata of all isolates (previously published and newly collected), including sampling location, year of collection and details on their hosts, can be found [here](../Datasets/S1_Data.csv). 

## Variant calling pipeline 

The pipeline starts with paired-end, raw WGS (short) reads and returns a VCF file with all biallelic SNPs for all *B.g. tritici* isolates in the [World](../Datasets/Datasets.md) dataset and 5 _B.g. secalis_ isolates used as outgroups for some analyses. The VCF file can be found at this [Zenodo repository](https://doi.org/10.5281/zenodo.13903934). 

### Software used
1. Python3
2. [fastp 0.23.2](https://github.com/OpenGene/fastp) 
3. [bwa 0.7.17-r1198-dirty](https://github.com/lh3/bwa)
4. [htslib 1.17](https://github.com/samtools/htslib/releases/tag/1.17)
5. [samtools 1.17](https://github.com/samtools/samtools/releases/tag/1.17)
6. [bcftools 1.17](https://github.com/samtools/bcftools/releases/tag/1.17)
7. [GATK 4.4.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.4.0.0)
8. R 4.2.1

### Workflow

1. `pipeline_with_gatk_statscsv.py` 
This takes as input the path to the raw paired-end fastq-files and reference genome along with some quality-based trimming and adapter trimming parameters. It returns a per-chromosome VCF file called by GATK `HaplotypeCaller` and some summary statistics about the mapping. It was run as an array job using the same input parameters for all samples. For example:
```bash
python3 pipeline_with_gatk_statscsv.py -ref GCA_900519115.1_2022_bgt_ref_mating_type.fa -minlen 50 -rw 5 -fw 1 -rq 20 -fq 20 -i CHNY072301_R1.fastq.gz
```
For further details on the steps executed in the script, refer to the Methods section of our paper. For more information on the input parameters of the script, run ```python3 pipeline_with_gatk_statscsv.py --help```  
2. Samples with coverage less than 15x were identified [(n = 26)](coverage_below_15) and excluded from all subsequent analyses. Mating types were assigned by comparing the coverage over the two alternate mating type genes:
```bash
while read p; do
        covga=$(samtools coverage -r LR026984.1_chr1:6734201-6735339 aln_${p}_compiled_marked_dup.bam | cut -f 7 | tail -1)  ## coverage over the 'reference' mating type
        covgb=$(samtools coverage -r Bgt_MAT_1_1_3 aln_${p}_compiled_marked_dup.bam | cut -f 7 | tail -1) ## coverage over the alternate mating type contig
        echo $p,$covga,$covgb >> mating_type_coverage
done < isolate_list
```
3. The VCF files for all remaining samples were combined (per-chromosome) using GATK `CombineGVCFs`. This step was performed in batches to avoid memory issues.
```bash
gatk --java-options "-Xmx44g -Xms44g" CombineGVCFs \
   -R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
   -V 2022_all_merge_$CHROMOSOME.list \
   -O 2022_all_combinegvcfs_$CHROMOSOME.vcf.gz
```
4. Variants were called on the combined VCF files for all chromosomes using GATK `GenotypeGVCFs`.
```bash
gatk --java-options "-Xmx20g" GenotypeGVCFs \
-R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
-V 2022+before2022+2023+ncsu_all_combinegvcfs_$CHROMOSOME.vcf.gz \
-O 2022+before2022+2023+ncsu_all_combinegvcfs_genotyped_$CHROMOSOME.vcf.gz \
--include-non-variant-sites \
-A StrandBiasBySample 
```
5. In order to decide threshold values for filtering variants, the distributions of site-based quality metrics were visualised. For each chromosome, SNPs were first selected from the output of step #4  using GATK `SelectVariants` `--select-type-to-include SNP` and their annotation values were written to a table using GATK `VariantsToTable`. Histograms were plotted for the genome-wide (chromosomes 1-11) values of the metrics `QD`, `FS`, `SOR`, `MQ`, `MQRankSum` and `ReadPosRankSum` using R `ggplot2`.
```bash
# select SNPs
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
     -V 2022+before2022+2023+ncsu_all_combinegvcfs_genotyped_$CHROMOSOME.vcf.gz \
     --select-type-to-include SNP \
     -O snps_ALL_$CHROMOSOME.vcf.gz

# output site-based quality metrics to a table
gatk VariantsToTable \
    -V snps_ALL_$CHROMOSOME.vcf.gz \
    -F CHROM -F POS -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum \
    -O snps_ALL_info_$CHROMOSOME.table

```

![2022+before2022+2023+ncsu_WG_gatk_info_distr-1](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/8e636ad7-1f92-4808-8250-f6d72ebaeb85)

6. Hard filtering for sites was done using GATK `VariantFiltration`. The threshold values (informed by the distributions plotted in step #5 and [GATK's recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)) were `QD < 10.0`, `FS > 50.0`, `MQ < 45.0`, `ReadPosRankSum < -5.0 || ReadPosRankSum > 5.0`.
```bash
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
7. `recode_multivcf_after_gatk_2024_new_clean.py` This script takes the VCF file produced by GATK in step #6 as input and recodes the 'GT' field value as '.' for all sites at which (a) depth of high-quality informative reads is less than 8 or (b) variant call is supported by < 90% of the high-quality informative reads. It outputs the recoded VCF (which can then be gzipped and indexed using `bgzip` and `tabix -p vcf`). 
```bash
python3 recode_multivcf_after_gatk_2024_new.py -i 2022+before2022+2023+ncsu_all_combinegvcfs_genotyped_varfilt_$CHROMOSOME.vcf.gz -ip $IN_PATH -o 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME -op $OUT_PATH -mc 8

bgzip 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf && tabix -p vcf 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf.gz
```
The script also outputs a csv file with the number of positions failing each of the above filters as well as the total number or missing positions and the total number of variant positions for each sample. The genome-wide distributions (obtained after summing over chromosomes 1-11) for these statistics were plotted using R `ggplot2`.

![2022+before2022+2023+ncsu_recoding_stats_distr_18042024](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/9bb03543-db20-477f-985c-0a04a16bb833)

The samples with > 200,000 'heterozygous positions', i.e. positions at which the variant support was < 90%, or those with a number of variants : number of heterozygous positions ratio of < 1  were excluded from all SNP-based analyses [(n=12)](2022+before2022+2023+ncsu_200k_hetpos_to_exclude_list.args). The exception was 96224, which had a low variant : heterozygous positions ratio because it is the reference isolate, and was thus retained in all analyses.

  
8. SNPs were selected from the chromosomal VCF files using GATK `SelectVariants` with options `--select-type-to-include SNP`, `--restrict-alleles-to ALL` and sites failing the filters from step #6 were excluded using the option `--exclude-filtered`.
```bash
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type_$CHROMOSOME.fa \
     -V 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf.gz \
     --select-type-to-include SNP \
     --restrict-alleles-to ALL \
     --exclude-filtered \
     -O 2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_$CHROMOSOME.vcf.gz
```
9. The resulting VCFs contained some [spanning deletions](https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele) denoted by '\*' . As these would have caused problems in downstream analyses, all sites with an '*' were removed using the script `recode_asterisk_count_snp.py` which returned a modified VCF file that was gzipped and indexed using `bgzip` and `tabix -p vcf` respectively.
```bash
python3 recode_asterisk_count_snp.py -i 2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_$CHROMOSOME.vcf.gz \
 -o 2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_$CHROMOSOME

bgzip 2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_$CHROMOSOME.vcf
tabix -p vcf 2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_$CHROMOSOME.vcf.gz
```
10. The VCF files (with filtered variants and no asterisks) for chromosomes 1-11, the alternate mating type locus and the mitochondrion were merged using the `concat` option in bcftools. The Bgt_Un "chromosome" with contigs not assigned to any other chromosomes was excluded from all further analyses.
```bash
bcftools concat -f 2022+before2022+2023+ncsu_snp_no_asterisk_11_chr_mt_MAT_list \ # list with the names of the VCF files
 -Oz -o 2022+before2022+2023+ncsu_recoded_snps_filtered_no_asterisk_11chr.vcf.gz
```
11. The resulting VCF files were then subset to include only biallelic SNPs. This step was performed using GATK `SelectVariants`. The tritici samples failing the "heterozygous" filter, as described in step #7, were excluded (n=12) `2022+before2022+2023+ncsu_200k_hetpos_to_exclude_list.args`. Further, tritici clones (as decided based on the [dist matrix analysis](../distance_matrix/distance_matrix.md) ) `2022+before2022+2023+ncsu_tritici_clones_to_exclude_list_18042024.args` were also excluded. The VCF file finally contained 568 _B.g. tritici_ (that made up the [World](../Datasets/Datasets.md) dataset) and 5 _B.g. secalis_ isolates that would be used as outgroups. The list of samples is given in `tritici_2022+before2022+2023+ncsu_no_clones_+_rye_old.args` (n=573).
```bash
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V ../project_data_prep/data/2022+before2022+2023+ncsu_recoded_snps_filtered_no_asterisk_11chr_mt_MAT.vcf.gz \
    --restrict-alleles-to BIALLELIC \
    --sample-name  tritici_2022+before2022+2023+ncsu_no_clones_+_rye_old.args\
    --select "AC>0 && AC<573" \
    -O tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz
```

