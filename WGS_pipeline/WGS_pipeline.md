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
8. 

### Workflow

1. `pipeline_with_gatk_statscsv.py` 
This takes as input the path to the raw fastq-files and reference genome along with some quality-based trimming and adapter trimming parameters. It returns a per-chromosome VCF file called by GATK HaplotypeCaller and some summary statistics about the mapping. It was run as an array job using the same input parameters for all samples. For example:
```
python3 ../scripts/pipeline_with_gatk_statscsv.py -ref GCA_900519115.1_2022_bgt_ref_mating_type.fa -minlen 50 -rw 5 -fw 1 -rq 20 -fq 20 -i /home/jjigis/projects/bgt_sequence_data/2023_collection/CHNY072301_R1.fastq.gz
```

2. The VCF files for all samples were combined (per-chromosome) using GATK CombineGVCFs
