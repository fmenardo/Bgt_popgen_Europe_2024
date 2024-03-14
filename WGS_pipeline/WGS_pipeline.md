# WGS pipeline

The pipeline starts with raw WGS (short) reads and returns a VCF file with all biallelic SNPs for all *B.g. tritici* isolates in the [World](../Datasets/Datasets.md) dataset. 

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
This takes as input the path to the raw fastq-files and reference genome along with some quality-based trimming and adapter trimming parameters. It returns a per-chromosome VCF file called by GATK HaplotypeCaller and some summary statistics about the mapping. It was run as an array job for 737 samples (ADD FILE WITH LIST OF ALL ACCESSIONS) using the same input parameters for all samples. For example:
```
python3 ../scripts/pipeline_with_gatk_statscsv.py -ref GCA_900519115.1_2022_bgt_ref_mating_type.fa -minlen 50 -rw 5 -fw 1 -rq 20 -fq 20 -i /home/jjigis/projects/bgt_sequence_data/2023_collection/CHNY072301_R1.fastq.gz
```
2. Samples with coverage less than 15x were identified (n = 26) and excluded from all subsequent analyses.
3. The VCF files for all remaining samples were combined (per-chromosome) using GATK `CombineGVCFs`.
4. Variants were called on the combined VCF files for all chromosomes using GATK `GenotypeGVCFs` with options `--include-non-variant-sites` and `-A StrandBiasBySample`.
5. In order to decide threshold values for filtering variants, the distribution of annotation values for SNPs were visualised. For each chromosome, SNPs were first selected from the output of step #3  using GATK `SelectVariants` `--select-type-to-include SNP` and their annotation values were written to a table using GTAK `VariantsToTable`. Histograms were plotted for the genome-wide values of the annotations `QD`, `FS`, `SOR`, `MQ`, `MQRankSum` and `ReadPosRankSum` using R `ggplot2`.

![2022+before2022+2023+ncsu_WG_gatk_info_distr-1](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/8e636ad7-1f92-4808-8250-f6d72ebaeb85)

6. Hard filtering for sites was done using GATK `VariantFiltration`. The threshold values (informed by the distributions plotted in step #5 and [GATK's recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)) were `QD < 10.0`, `FS > 50.0`, `MQ < 45.0`, `ReadPosRankSum < -5.0 || ReadPosRankSum > 5.0`. 
7. `recode_multivcf_after_gatk_2024_new_clean.py` This script takes the genotyped VCF file produced by GATK in step #3 as input and recodes the 'GT' field value as '.' for all sites at which (a) depth of high-quality informative reads < user-defined minimum depth or (b) variant call is supported by < 90% of the high-quality informative reads. It also outputs a file with the final number of missing sites for each sample as well as other statistics.
8. 

