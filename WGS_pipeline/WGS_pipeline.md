# WGS pipeline

The pipeline starts with raw WGS (short) reads and returns a VCF file with all biallelic SNPs for all *B.g. tritici* isolates in the [World](../Datasets/Datasets.md) dataset. 

### Software used 
1. GATK v 4.4 (link)
2. xyz

### Workflow

1. `pipeline_with_gatk_statscsv.py` 
This takes as input the path to the raw fastq-files and reference genome along with some quality-based trimming and adapter trimming parameters. It returns a per-chromosome VCF file called by GATK HaplotypeCaller and some summary statistics about the mapping. It was run as an array job using the same input parameters for all samples. For example:
```
python3 ../scripts/pipeline_with_gatk_statscsv.py -ref GCA_900519115.1_2022_bgt_ref_mating_type.fa -minlen 50 -rw 5 -fw 1 -rq 20 -fq 20 -i /home/jjigis/projects/bgt_sequence_data/2023_collection/CHNY072301_R1.fastq.gz
```

2. The VCF files for all samples were combined (per-chromosome) using GATK CombineGVCFs
