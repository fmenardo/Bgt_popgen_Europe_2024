# TreeMix
We used TreeMix to construct a population tree while allowing for admixture. We used 5 _B.g. secalis_ isolates as outgroup along with the [Europe+ dataset](../Datasets/Datasets.md) for analysis. We chose the fineSTRUCTURE level-8 of population classification which divides the _B.g. tritici_ Europe+ dataset into 19 populations. The population assignment can be seen [here](fs8_tritici_ext_eur+secalis_old.clust) The VCF file containing the 420 samples, filtered for singletons and missing data was generated using
```bash
gatk SelectVariants \
    -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
    -V tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz \
    -XL MT880591.1 \  #exclude mitochondria
    -XL Bgt_MAT_1_1_3 \  #exclude mating type contig
    -sample-name tritici_extended_europe_2022+before2022+2023+ncsu_+_old_rye.args \  #list of samples to include
    --select "AC>1 && AC<419" \  #remove singletons
    --max-nocall-fraction 0 \  #remove missing data
    -O tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss.vcf.gz
```
This was followed by LD based pruning as the program requires independent markers.
```bash
vcf=tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss.vcf.gz
out=tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss_25kb_0.1_LDp

# make bed file with variant id
plink --allow-extra-chr --vcf "${vcf}" --set-missing-var-ids @:#\$1,\$2 --make-bed --double-id --out prefilter --threads 2
# identify variants within 25>Kb and r2 > 0.1
plink --indep-pairwise 25kb 1 0.1 --bed prefilter.bed --bim prefilter.bim --fam prefilter.fam --double-id --allow-extra-chr --threads 2  --out $out
# select those variants from vcf
plink --bed prefilter.bed --bim prefilter.bim --fam prefilter.fam --extract "${out}.prune.in" --recode vcf-iid bgz --out $out --double-id --allow-extra-chr --threads 2
```
The chromosomes were renamed for the subsequent steps using
``` bash
zcat tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss_25kb_0.1_LDp.vcf.gz | sed 's/LR0269[0-9][0-9].1_chr//g' > tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss_25kb_0.1_LDp_renamed_chr.vcf && bgzip tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss_25kb_0.1_LDp_renamed_chr.vcf && tabix -p vcf tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss_25kb_0.1_LDp_renamed_chr.vcf.gz
```
This VCF file and the clust file `fs8_tritici_ext_eur+secalis_old.clust` are used as input for the script `vcf2treemix_mod.sh` that prepares the input files needed for treemix. This script requires the python script `plink2treemix.py` available for download [here](https://bitbucket.org/nygcresearch/treemix/downloads/).

Treemix was then run allowing for 1-8 migration edges with 5 replicates each
```bash
#!/bin/bash  
infile=tritici_ext_eur+secalis_biallelic_snps_old_no_sing_no_miss_25kb_0.1_LDp_renamed_chr.treemix.frq.gz
outstem=output_dir/fs8_no_sing_no_miss_LDp

for m in {1..8}
        do
        for i in {1..5}
                do
                treemix -i "${infile}" -o "${outstem}.r${i}.m${m}" -global -m ${m} -root secalis -bootstrap -k 500 -seed ${RANDOM}
                done
        done
```
The [trees](fs8_m1-6_maxLL_tree.pdf) and [residuals](fs8_m1-6_maxLL_residuals.pdf) for the maximum-likelihood replicates of each migration edge(m=1-6) were plotted using the script `plotting_funcs.R` (included in the treemix package).
```R
source("plotting_funcs.R")
maxll_reps <- c(1,4,1,2,3,1)  #max-likelihood replicate

par(mfrow=c(2,3))
for (m in 1:6){
  plot_tree(paste0("fs8_no_sing_no_miss_LDp.r",maxll_reps[m],".m",m))
  title(paste(m, "edges"))
}

par(mfrow=c(2,3))
for (m in 1:6){
  plot_resid(paste0("fs8_no_sing_no_miss_LDp.r",maxll_reps[m],".m",m), "pops_fs8", cex = 0.8)
  title(paste(m, "edges"))
}
```

The R package [OptM](https://cran.r-project.org/web/packages/OptM/index.html) was used to assess the optimum number of migration edges to include. 
```R
library(OptM)
test.optM <- optM("output_dir/")
plot_optM(test.optM, method = "Evanno")
```
The [optM plots](optM_plots.pdf) showed the highest increase in likelihood for m=2 which was thus chosen as the 'best m'. 

We then performed 500 bootstrap replicates of treemix with m=2. We used [PHYLIP](https://phylipweb.github.io/phylip/) to construct a consensus tree based on the bootstrap replicates. Treemix was then run again by loading this consensus tree and setting m=2. These steps were performed using the wrapper script `Treemix_bootstrap.sh` as provided in the R package [BITE](https://github.com/marcomilanesi/BITE). This R package was also used to visualise the [consensus tree](fs8_m2_bs100_cons_tree.pdf).
```R
library(BITEV2)
BITEV2::treemix.bootstrap(in.file = "fs8_m2_100bs",out.file = "bite_out_m2",phylip.file = "fs8_m2_100bs_outtree.newick",
                  nboot=100, plotmig = T, boot.legend.location='bottomright')

```
## Software
```
GATK 4.4.0.0
PLINK 1.90b7
TreeMix 1.13
htslib 1.17
PHYLIP 3.697

R packages:
- OptM 0.1.6
- BITE 2.1.0
- BITE dependencies
    SNPRelate >= 1.32.2
    ggplot2 >= 3.4.2
    gridExtra >= 2.3
    stringr >= 1.5.0
    shiny >= 1.7.4.1
    dplyr >= 1.1.2
    plotly >= 4.10.2
    ggrepel >= 0.9.3
    data.table >= 1.14.8
    RCircos >= 1.2.2
    poolfstat >= 2.1.2
    RColorBrewer >= 1.1-3
```
