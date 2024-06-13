# Relate
We used Relate ([Speidel et al. 2019](https://doi.org/10.1038/s41588-019-0484-x)) to reconstruct genome-wide genealogies for the dataset [Europe+_recent](../Datasets/Datasets.md).

For each chromosome we selected sites without missing data and excluded all sites axcept biallelic SNPs and invariant sites.

```
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf.gz \
     --exclude-filtered \
     --sample-name tritici_2022-2023_Bgs.args\
     --remove-unused-alternates\
     -O 2022-2023+out_$CHROMOSOME.vcf.gz

gatk SelectVariants \
     -R ~/projects/vcf_project_tritici/GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V 2022-2023+out_$CHROMOSOME.vcf.gz \
     -O 2022-2023+out_$CHROMOSOME.BISNP.vcf.gz\
     --restrict-alleles-to BIALLELIC\
     --select-type-to-include SNP

gatk SelectVariants \
     -R ~/projects/vcf_project_tritici/GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V 2022-2023+out_$CHROMOSOME.vcf.gz \
     -O 2022-2023+out_$CHROMOSOME.NOVAR.vcf.gz\
     --select-type-to-exclude MNP\
     --select-type-to-exclude INDEL\
     --select-type-to-exclude SYMBOLIC\
     --select-type-to-exclude MIXED\
     --select-type-to-exclude SNP \
     --select-type-to-include NO_VARIATION


bcftools concat -a 2022-2023+out_$CHROMOSOME.NOVAR.vcf.gz 2022-2023+out_$CHROMOSOME.BISNP.vcf.gz -D -O z > 2022-2023+out_$CHROMOSOME.BISNP+NOVAR.vcf.gz

```
Then we prepared the input files (*.sample, *.poplabels, *.haps, and *.dist) for Relate, in this step we also polarized SNPs using 5 B.g. secalis strains as outgroup.

```
python prep_Relate_input.py -vcf 2022-2023+out_$CHROMOSOME.BISNP+NOVAR.vcf.gz -anc tritici_outgroups.args -meta 2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv -o 2022_2023_$CHROMOSOME -onc 2022_2023
```



