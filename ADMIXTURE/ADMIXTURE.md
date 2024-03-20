# ADMIXTURE
We performed ADMIXTURE analysis on the [*World*](../Datasets/Datasets.md) dataset. Since the model assumes linkage equilibrium among markers, we pruned our SNP dataset using PLINK, as shown below:
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
This resulted in a VCF file with 155,979 SNPs over 568 samples. ADMIXTURE was run using the script `admixture.py` 
```
python3 admixture.py -vcf 2022+before2022+2023+ncsu_no_sing_25kb_0.1_LDp.vcf.gz -c 4 -o tritici_ALL_25kb_0.1_LDp -r 1 -K 15
```
We ran ADMIXTURE for 15 K values (1-15) and evaluated the cross-validation error at each K. `tritici_ALL_25kb_0.1_LDp_admixture_K_r1.log` shows that the lowest CV error value was obtained for K=7, as also seen in this [plot](CV_error_admixture_r1.pdf).

We plot the ancestry proportions for all individuals for K=4-8 as barplots using the script `plot_admixture_revised_colours_r1.R` interactively. The plots can be found [here](admx_bar_plots_r1_k4-8_revised_colours.pdf).
