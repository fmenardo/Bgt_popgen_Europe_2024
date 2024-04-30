# ADMIXTURE
We performed [ADMIXTURE](https://dalexander.github.io/admixture/download.html) analysis on the [*World*](../Datasets/Datasets.md) dataset. Since the model assumes linkage equilibrium among markers, we pruned our SNP dataset using PLINK, as shown below:
```
vcf=2022+before2022+2023+ncsu_11chr_no_sing_maxmiss0.1.vcf.gz
out=2022+before2022+2023+ncsu_no_sing_25kb_0.1_LDp

# LD pruning

# make bed file with variant id
srun ~/data/plink --allow-extra-chr --vcf "${vcf}" --set-missing-var-ids @:#\$1,\$2 --make-bed --double-id --out prefilter --threads 2
# identify variants within 25 Kb and r2 > 0.1
srun ~/data/plink --indep-pairwise 25kb 1 0.1 --bed prefilter.bed --bim prefilter.bim --fam prefilter.fam --double-id --allow-extra-chr --threads 2  --out $out
# select those variants from vcf
srun ~/data/plink --bed prefilter.bed --bim prefilter.bim --fam prefilter.fam --extract "${out}.prune.in" --recode vcf-iid bgz --out $out --double-id --allow-extra-chr --threads 2
```
This resulted in a VCF file with 156,047 SNPs over 568 samples. ADMIXTURE was run for 10 replicates parallely using the script `admixture_parallel.py`, where each replicate had a different random seed (`random_seeds_10_runs`). An example of how the script was called is shown here-
```
python3 admixture_parallel.py -B tritici_ALL_25kb_0.1_LDp -c 4 -o tritici_ALL_25kb_0.1_LDp_r1 -r 1 -K 10
```
We ran ADMIXTURE for 10 K values (1-10) and evaluated the cross-validation error from all replicates at each K. The CV error values for all ADMIXTURE runs is reported in `CV_errors_10reps_rep_col` and shown in this [boxplot](CV_error_boxplot_10_reps.pdf). 

For each value of K, we selected the ADMIXTURE run with the least cross validation error. These minimum CV errors are plotted [here](min_CV_error_10_reps.pdf). We plot the ancestry proportions of the corresponding runs for all individuals for K=4-9 as barplots using the script `plot_admixture_best_runs.R` interactively. The plots can be found [here](k4-9_admixture_barplot_w18.pdf).
