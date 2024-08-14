# Demographic inference

## Data preparation
### Analysis not considering the probability of misidentification of derived alleles

We focus on two populations, the list of samples can be found in `N_EUR2_Bgs.args` and `E_EUR2_Bgs.args`. Both lists contain 5 Bgs isolates to be used as outgroup to polarize SNPs.

For each chromosome we selected sites without missing data and excluded all sites axcept biallelic SNPs and invariant sites.
For example for population N_EUR2:

```
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V 2022+before2022+2023+ncsu_covg15_recoded_$CHROMOSOME.vcf.gz \
     --exclude-filtered \
     --sample-name N_EUR2_Bgs.args\
     --remove-unused-alternates\
     -O N_EUR2_$CHROMOSOME.vcf.gz

gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V N_EUR2_$CHROMOSOME.vcf.gz \
     -O N_EUR2_$CHROMOSOME.BISNP.vcf.gz\
     --restrict-alleles-to BIALLELIC\
     --select-type-to-include SNP

gatk SelectVariants \
     -R ~/projects/vcf_project_tritici/GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V N_EUR2_$CHROMOSOME.vcf.gz \
     -O N_EUR2_$CHROMOSOME.NOVAR.vcf.gz\
     --select-type-to-exclude MNP\
     --select-type-to-exclude INDEL\
     --select-type-to-exclude SYMBOLIC\
     --select-type-to-exclude MIXED\
     --select-type-to-exclude SNP \
     --select-type-to-include NO_VARIATION


bcftools concat -a 2022-2023+out_$CHROMOSOME.NOVAR.vcf.gz 2022-2023+out_$CHROMOSOME.BISNP.vcf.gz -D -O z > 2022-2023+out_$CHROMOSOME.BISNP+NOVAR.vcf.gz

python parse_vcf_genomics.py -o $CHROMOSOME.N_EUR2.genomic.fa -vcf N_EUR2_$CHROMOSOME.BISNP+NOVAR.vcf.gz

```
This code generate a fasta file containing all sites without missing data in Bgt which are invariant or biallelic SNPs. Additionally only sites that can be polarized are included (i.e., at least one outgroup isolates has a high quality call at the site, and the site is invariant in Bgs). The Bgs isolates are collapsed in one consensus sequences named ANC (anyway only sites that are monomorphic in Bgs are included).

### Analysis including the probability of misidentification of derived alleles
