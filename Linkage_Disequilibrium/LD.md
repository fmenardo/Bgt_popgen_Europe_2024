# Linkage Disequilibrium (LD) decay

We calculated the Linkage Disequilibrium decay in the _Europe+_recent_ dataset separately for the five separate populations of fine structure level 4 subdivision as well as for all samples together.
- We used the *tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz* as input file and filtered for maximum 10% missing data while creating the separate vcf files for each population.
- We then transformed the files to plink (bed) format 
- and run LD decay, using plink2 in 10 kb windows.
- We finally plotted the results

The script used is this:
```
bash LD_decay_level4_10kb_100kb.sh
```

The scripts used for the plots are:
```
LD_decay_plink_10kb_onePlot_legend_fixed.R
```

**R packages**
```
R 4.3.1
ggplot2 3.4.4
dplyr 1.1.2
stringr 1.5.0
```
