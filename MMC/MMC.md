# Demographic inference

## Data preparation
### Data preparation for analysis not considering the probability of misidentification of derived alleles

We focus on two populations, the list of samples can be found in `N_EUR2_Bgs.args` and `E_EUR2_Bgs.args`. Both lists contain 5 Bgs isolates to be used as outgroup to polarize SNPs `list_Bgs.args`.

For example for population N_EUR2:

```bash
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
```
We then merge all chromosomes in one file and process it. This code generates a fasta file containing all biallelic SNPs without missing data in Bgt. Additionally only sites that can be polarized are included (i.e., at least one outgroup isolates has a high quality call at the site, and the site is invariant in Bgs). The Bgs isolates are collapsed in one consensus sequences named ANC (anyway only sites that are monomorphic in Bgs are included).

```bash
ls | grep BISNP+NOVAR.vcf.gz > list_vcf

bcftools concat -f list_vcf -O z -o E_EUR2_all_chr_NOVAR+BISNP.vcf.gz

python ../parse_vcf_pol.py -vcf E_EUR2_all_chr_NOVAR+BISNP.vcf.gz -o E_EUR2_all_chr_BISNP.fasta -anc ../list_Bgs.args

```

### Data preparation for analysis including the probability of misidentification of derived alleles

For this analysis we have to keep sites that are polymorphic in Bgs. This code generate one fasta file per chromosome containing all SNPs and invariant sites. The Bgs isolates are included in the alignments.

```bash
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V E_EUR2_$CHROMOSOME.vcf.gz \
     -O E_EUR2_$CHROMOSOME.MULTISNP.vcf.gz\
     --select-type-to-include SNP

bcftools concat -a E_EUR2_$CHROMOSOME.NOVAR.vcf.gz E_EUR2_$CHROMOSOME.MULTISNP.vcf.gz -D -O z > E_EUR2_$CHROMOSOME.MULTISNP+NOVAR.vcf.gz

python ../parse_vcf_genomics.py -o $CHROMOSOME.E_EUR2_MULTIALLELIC.genomic.fa -vcf E_EUR2_$CHROMOSOME.MULTISNP+NOVAR.vcf.gz

```
## Analysis

The fasta files generated above are used for subsequent analysis where we test if the assumption of small variance in reproductive success is appropriate for Bgt. The script for the analysis not considering the probability of misidentification of derived alleles is `mmc_vs_km_noalleleconf.R`, and that where the misidentification is accounted for is `mmc_vs_km_alleleconf.R`. Additional scripts `funct_and_consts.R` and `model_select.R`, and the git repo https://github.com/fabfreund/usfs_mmc/ are required for this analysis. `tajD_from_sfs.R` is a R function to compute Tajima's D.
