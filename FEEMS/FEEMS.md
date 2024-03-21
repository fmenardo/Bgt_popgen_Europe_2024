# Fast Estimation of Effective Migration Surfaces (FEEMS)
We use the method developed by [Marcus *et al* (2021)](https://elifesciences.org/articles/61927) to estimate effective migration surfaces for our [Europe+_2022_2023 dataset](../Datasets/Datasets.md).  
We used the conda installation instructions and adapt the code provided by the authors of the tool [here](https://github.com/NovembreLab/feems).

#### Preparation of input files 
1. We constructed a discrete global grid of a suitable resolution over our sampling using the R package `dggridR`. The grid files produced are `dgg_res7.shp`, `dgg_res7.shx` and `dgg_res7.dbf`. 
```
%% R
library(dggridR)
### constructing grid for region of interest

dggs <- dgconstruct(metric=TRUE, resround='nearest', aperture = 4, topology = "TRIANGLE",spacing = 50)
mygrid <- dgrectgrid(dggs, minlat = 27.260104, minlon = -27.406657,
                     maxlat = 68.992133, maxlon = 47.651937, savegrid = "dgg_res7.shp")
```
2. We subset the all-sample, biallelic SNP VCF file to include only samples from the Europe+_2022_2023 dataset and filtered out all singletons and missing data using GATK `SelectVariants`.
```
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
    -XL MT880591.1 \   # exclude mitochondrion
    -XL Bgt_MAT_1_1_3 \   # exclude the mating type locus
    -sample-name tritici_2022_2023.args \   # list of samples to include
    --select "AC>1 && AC<255" \   # remove singletons
    --max-nocall-fraction 0 \    # remove all missing data
    -O tritici_2023_no_sing_no_miss.vcf.gz
```
3. Next, we performed LD-based pruning using PLINK. At each step (1 SNP), the program notes the pairs of variants in the window (25kb) that exceed the correlation value (0.1) and removes variants from the window until no such pairs remain. 
```
vcf=tritici_2022+2023_no_sing_no_miss.vcf.gz
out=tritici_2022+2023_no_sing_no_miss_LDp_25kb_0.1

# LD pruning 

# make bed file with variant id
plink --allow-extra-chr --vcf $vcf --set-missing-var-ids @:#\$1,\$2 --make-bed --double-id --out prefilter --threads 2
# identify variants within 25Kb and r2 > 0.1
plink --indep-pairwise 25kb 1 0.1 --bed prefilter.bed --bim prefilter.bim --fam prefilter.fam --double-id --allow-extra-chr --threads 2  --out $out
# select those variants from vcf
plink --bed prefilter.bed --bim prefilter.bim --fam prefilter.fam --extract "${out}.prune.in" --recode vcf-iid bgz --out $out --double-id --allow-extra-chr --threads 2
# make bed
plink --vcf "${out}.vcf.gz" --double-id --allow-extra-chr --biallelic-only --make-bed --out $out
```
4. The sampling coordinates (in long lat format) for all samples were extracted from the [metadata file](../Datasets/2022+before2022+2023+ncsu_metadata+fs+admxK7_19032024.csv) and ordered according to the order of samples in the .fam file produced by PLINK to produce `tritici_2022+2023_long_lat.coord`. An outer polygon surrounding the sampling region was given by `extended_europe.outer`.

#### Running FEEMS
We ran FEEMS using the script `run_feems.py`.
```
python3 run_feems.py -i tritici_2022+2023_no_sing_no_miss_LDp_25kb_0.1 -sc tritici_2022+2023_long_lat.coord -oc extended_europe.outer -dgg dgg_res7.shp
```
This script produced the [cross validation error plot](tritici_2022+2023_no_sing_no_miss_LDp_25kb_0.1_CV_error_dgg_res7.pdf), the [FEEMS plot for the best lambda](tritici_2022+2023_no_sing_no_miss_LDp_25kb_0.1_feems_plot_dgg_res7.pdf) and the [plot comparing migration surfaces over four lambda values](tritici_2022+2023_no_sing_no_miss_LDp_25kb_0.1_feems_plot_dgg_res7_lambda_compare.pdf).
