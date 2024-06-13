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
Then we prepared the input files (*.sample, *.poplabels, *.haps, and *.dist) for Relate, in this step we also polarized SNPs using 5 B.g. secalis strains as outgroup. The recombination map was generated [here](../recombination_map/recombinaion_map.md).

```
python prep_Relate_input.py -vcf 2022-2023+out_$CHROMOSOME.BISNP+NOVAR.vcf.gz -anc tritici_outgroups.args -meta 2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv -o 2022_2023_$CHROMOSOME -onc 2022_2023
```
Overall we obained 1306307 SNPs over the 11 chromosomes.

We ran the first step of Relate for each chromosome:

```
Relate --mode All -N 3290 -m 0.0000005 --haps 2022-2023+out_$CHROMOSOME.haps --map ../recombination_map/THUN12x96224_4Relate_chr$CHROMOSOME\_Mb_recombination_rates.txt --sample 2022_2023.sample --seed 54455 -o 2022_2023_theta_chr$CHROMOSOME --dist 2022-2023+out_$CHROMOSOME.dist
```

The mutation rate (-u) was obained from [Sotiropulos et al. 2022](https://doi.org/10.1038/s41467-022-31975-0), and the starting guess for Ne (-N) was calculated as theta/2*u with the script `calc_theta.R`.

In the first step Relate estimates genealogies under a constant population size model, in a second step we used the Relate script `EstimatePopulationSize.sh` to rerestimate population sizes and brach lengths under a stepwise population size model. 

```
EstimatePopulationSize.sh  -i 2022_2023_theta -o 2022_2023_theta_popsize --mu 0.0000005 --years_per_gen 1 --first_chr 1 --last_chr 11 --seed 56464 --poplabels 2022_2023.poplabels --threads 10``
```
And finall we convert the output of Relate in a tree sequence format: 
```
for i in {1..11}
do
 gzip -d 2022_2023_theta_popsize_chr$i.anc.gz
 gzip -d 2022_2023_theta_popsize_chr$i.mut.gz
 RelateFileFormats --mode ConvertToTreeSequence -i 2022_2023_theta_popsize_chr$i -o 2022_2023_theta_popsize_chr$i
done

```


