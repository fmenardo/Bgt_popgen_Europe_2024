# Haplotype analysis for AvrPm17

The genomic location and structure of the gene encoding the powdery mildew effector AvrPm17 was characterized in [Müller, Kunz et al 2022](https://doi.org/10.1073/pnas.2108808119). The effector was found to be encoded by a paralogous gene pair on chromosome 1, and variation in gene copy number was also reported. 

In order to characterize the number of gene copies and the corresponding haplotypes in all isolates of the [Europe+ dataset](../Datasets/Datasets.md), we performed the following steps:
1. Mapped raw sequencing reads to a reference genome containing just the first copy of the gene (Chr1: 4365017-4365402) (+-2Kb) and then called variants using steps similar to [WGS_pipeline](../WGS_pipeline/WGS_pipeline.md). The reference was prepared using the script `extract_avrpm17_region_from_ref.py`. The required index and dict files of the reference were generated using:
```bash
bwa index avrpm17_locus+-2Kb.fa
samtools faidx bwa avrpm17_locus+-2Kb.fa
gatk CreateSequenceDictionary -R avrpm17_locus+-2Kb.fa 
```
The script to perform the subsequent steps is `pipeline_with_gatk_statscsv_clean_avrpm17.py` which can be called for each isolate separately, for example:
```bash
python3 pipeline_with_gatk_statscsv_clean_avrpm17.py -ref avrpm17_locus+-2Kb.fa -minlen 50 -rw 5 -fw 1 -rq 20 -fq 20 -i 94202_R1.fastq.gz
```
2. The bam file produced for each isolate in step #1 was used to obtain the read depth for the genic region using `samtools coverage`.
```bash
#!/bin/bash
while read p; do
        covg=$(samtools coverage -r avrpm17_locus+-2Kb:2001-2386 ${p} | cut -f 7 | tail -1)  ### field 7 is mean depth of coverage
        iso=$(echo ${p} | sed 's/mapped_aln_//g' | sed 's/_compiled_marked_dup.bam//g')  ### get name of isolate from bam file name
        echo $iso,$covg >> avrpm17_covg  ### file with gene coverage
done < avrpm17_bam_list
```
The depth in the genic region was divided by the genome-wide average coverage (as obatained [here](../WGS_pipeline/WGS_pipeline.md)) to get the number of gene copies in each isolate. 
```R
#R
## read in genome-wide stats
gw_stats <- read.csv("2022+before2022+2023+ncsu_737_gatkpl_stats_new.csv")
gw_covg <- gw_stats[,c(1,5)]

## gene-coverage
gene_covg <- read.csv("avrpm17_covg.csv", header = FALSE)

covg <- merge(gw_covg, gene_covg, by.x = "Isolate", by.y = "V1")
covg$ratio <- covg$V2 / covg$average.genome.wide.coverage
hist(covg$ratio)
```
![image](https://github.com/fmenardo/Bgt_popgen_Europe_2024/assets/90404355/407a9a54-d49c-4faa-8f1f-eba248c161b8)


3. The per-sample VCF files generated in step #1 were combined, genotyped and filtered using GATK.
```bash
## combine vcfs
gatk --java-options "-Xmx8g -Xms8g" CombineGVCFs \
   -R avrpm17_locus+-2Kb.fa \
   -V tritici_ext_eur_avrpm17_vcf.list \
   -O tritici_ext_eur_avrpm17_combined.vcf.gz

## joint genotyping
gatk --java-options "-Xmx8g" GenotypeGVCFs \
-R avrpm17_locus+-2Kb.fa \
-V tritici_ext_eur_avrpm17_combined.vcf.gz \
-O tritici_ext_eur_avrpm17_combined_genotyped.vcf.gz \
--include-non-variant-sites \
-A StrandBiasBySample 
 
## site filtration
gatk --java-options "-Xmx8g" VariantFiltration \
   -R avrpm17_locus+-2Kb.fa \
   -V tritici_ext_eur_avrpm17_combined_genotyped.vcf.gz \
   -O tritici_ext_eur_avrpm17_combined_genotyped_filtered.vcf.gz \
   --filter-name "QualByDepth10" \
   --filter-expression "QD < 10.0" \
   --filter-name "FisherStrand50" \
   --filter-expression "FS > 50.0" \
   --filter-name "RMSMapQual45" \
   --filter-expression "MQ < 45.0" \
   --filter-name "ReadPosRankSumTest5" \
   --filter-expression "ReadPosRankSum < -5.0 || ReadPosRankSum > 5.0"
```

4. To distinguish between haplotypes of the different gene copies, we used sites that were heterozygous, i.e. those that had less than 90% support for a variant call. We recoded the final mutli-sample VCF file obtained from step #3 to have diploid genotype calls and reflect our definition of "heterozygous" using the script `recode_to_diploid.py` called as
```bash
python3 recode_to_diploid.py -i tritici_ext_eur_avrpm17_only_gene_combined_genotyped_filtered.vcf.gz -ip map_call_pl/ -o tritici_ext_eur_avrpm17_only_gene_diploid -op map_call_pl/ -mc 8
bgzip tritici_ext_eur_avrpm17_only_gene_diploid.vcf && tabix -p vcf tritici_ext_eur_avrpm17_only_gene_diploid.vcf.gz
```
5. This diploid VCF file was phased using [WhatsHap](https://whatshap.readthedocs.io/en/latest/). This required one bam file containing the alignments of all the samples. This was generated using:
```bash
samtools merge -b avrpm17_bam_list -o avrpm17_ext_eur_merged.bam && samtools index avrpm17_ext_eur_merged.bam
```
And the phasing was performed using:
```bash
whatshap phase -o avrpm17_only_gene_phased.vcf.gz --reference=avrpm17_locus+-2Kb.fa tritici_ext_eur_avrpm17_only_gene_diploid.vcf.gz avrpm17_ext_eur_merged.bam
```
6. WhatsHap was able to phase all samples that had at least two biallelic heterozygous sites. For the samples it could not phase (as they had just one usable heterozygous site), we created two haplotypes - one with the reference and one with the alternate allele using `bcftools consensus`.
```bash
while read p; do
        bcftools consensus -a '-' --missing 'N' --haplotype R -s ${p} -f avrpm17_locus+-2Kb.fa avrpm17_only_gene_phased.vcf.gz | sed "s/avrpm17_locus+-2Kb/${p}_1_avrpm17/g" > ${p}_avrpm17_hap1.fa
done < unphased_het.list

while read p; do
        bcftools consensus -a '-' --missing 'N' --haplotype A -s ${p} -f avrpm17_locus+-2Kb.fa avrpm17_only_gene_phased.vcf.gz | sed "s/avrpm17_locus+-2Kb/${p}_2_avrpm17/g" > ${p}_avrpm17_hap2.fa
done < unphased_het.list
```
For samples that could be phased, we generated consensus calls using:
```bash
while read p; do
        bcftools consensus -a '-' --missing 'N' --haplotype 1pIu -s ${p} -f avrpm17_locus+-2Kb.fa avrpm17_only_gene_phased.vcf.gz | sed "s/avrpm17_locus+-2Kb/${p}_1_avrpm17/g" > ${p}_avrpm17_hap1.fa
done < all_sample_except_one_hetpos_list.args

while read p; do
        bcftools consensus -a '-' --missing 'N' --haplotype 2pIu -s ${p} -f avrpm17_locus+-2Kb.fa avrpm17_only_gene_phased.vcf.gz | sed "s/avrpm17_locus+-2Kb/${p}_2_avrpm17/g" > ${p}_avrpm17_hap2.fa
done < more_than_one_het_pos.list
```
Some samples (n=7) had a missing genotype call in the genic region. Since this would cause ambiguity in haplotype assignment, we excluded these 7 samples (14 haplotypes). The remaining 516 haplotypes were combined to generate a multi-fasta file:
```bash
xargs cat < hap_list_no_missing > avrpm17_msa_ext_eur_no_missing.fa
```
7. Haplotypes for just the genic region (excluding flanking 2Kb on either side) were generated using `extract_and_rev_comp_fasta.py`. Haplotypes for the coding region (excluding the intron) were generated `extract_cds_translate.py`. Since this copy of the gene is on the reverse strand, the amino acid haplotype was generated by translating the reverse complement. The haplotype files are `avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa` (DNA) and `avrpm17_ext_eur_only_cds_trns_unique_no_miss.fa` (amino acid). To retain only unique haplotypes from each sample, we used `get_unique_hap_in_cds_no_miss.py`.
8. Identical haplotypes across samples were grouped using `classify_haplotypes.py`. We classified haplotypes at the nucleotide level, using the complete coding sequence (`hap_class_cds_RC_no_miss.csv`), and at the aminoacid level, using only the mature protein (i.e., after removal of signal petide; (`hap_class_mp_trns_no_miss.csv`)) .
9. Finally we inferred and plotted the haplotype network with `plot_hap_networks.R`.
   The R environment for this last steop was:
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
[1] igraph_1.5.0    dplyr_1.1.2     maps_3.4.1      ggplot2_3.4.2   isoRelate_0.1.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.12       rstudioapi_0.15.0 magrittr_2.0.3    tidyselect_1.2.0  munsell_0.5.0     colorspace_2.1-0  R6_2.5.1          rlang_1.1.3      
 [9] foreach_1.5.2     fansi_1.0.4       tools_4.2.1       grid_4.2.1        gtable_0.3.3      utf8_1.2.3        cli_3.6.1         withr_2.5.0      
[17] iterators_1.0.14  tibble_3.2.1      lifecycle_1.0.3   vctrs_0.6.3       ggnetwork_0.5.12  codetools_0.2-19  glue_1.6.2        compiler_4.2.1   
[25] pillar_1.9.0      generics_0.1.3    scales_1.3.0      pkgconfig_2.0.3 
```
