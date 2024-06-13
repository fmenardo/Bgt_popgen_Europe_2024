# Identity by descent clusters analysis for AvrPm17

In the genome-wide analysis with isoRelate we found that in many populations there is an excess of IBD pairs at the locus containining AvrPm17. Here we explore which pairs of isolates are in IBD over this locus.

First we run again isoRelate, only this time we focus on the locus of AvrPm17 and we use the compete Europe+_recent dataset (we do not run separate analysis for different populations as we did for the genome-wide analysis)

We prepare the input files:

```
gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe+_recent_chr1_avrpm17.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals THUN12x96224_genetic_map_in_cM_+_phy_distance.ambiguous_intervals.list \
     --intervals LR026984.1_chr1:3000000-6000000 \
     --select "AF > 0.05" \
     --sample-name ../Datasets/tritici_recent_extended_europe_2022+2023+ncsu.args


plink --allow-extra-chr --vcf Europe+_recent_chr1_avrpm17.vcf.gz --recode 12 --double-id --out BgtE+r_chr1_avrpm17 --threads 1

awk '$5="1" && $6="0"' BgtE+r_chr1_avrpm17.ped >  BgtE+r_chr1_avrpm17_mod.ped
sed -E 's/\S+_chr//g' BgtE+r_chr1_avrpm17.map > BgtE+r_chr1_avrpm17_mod.map


python ../isoRelate/add_cM_to_map.py -map BgtE+r_chr1_avrpm17_mod.map -rec ../recombination_map/THUN12x96224_genetic_map_in_cM_+_phy_distance -o BgtE+r_chr1_avrpm17_mod
```
We run isoRelate

```
Rscript run_ibd_step1.R -o BgtE+r_avrpm17 -p BgtE+r_chr1_avrpm17_mod.ped -m BgtE+r_chr1_avrpm17_mod_cM.map -c 5 -C 2
```

