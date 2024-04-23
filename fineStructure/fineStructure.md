# fineStructure
## Data preparation

For the fineStructure analyses we used the [Europe+](../Datasets/Datasets.md) datatset with 415 individuals sampled from Europe, the Middle East and Caucasus.
We kept only SNPS without any missing data, as fineStructure cannot handle them: 

```
gatk-4.4.0.0/gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe_large_tritici_no_clones_no_miss_new.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals MT880591.1 \                    #exclude mitochondrion
     --exclude-intervals LR026995.1_Un \                 #exclude conting with sequenced non assined to any chromosome
     --exclude-intervals Bgt_MAT_1_1_3 \                 #exclude contig with alternative mating type
     --sample-name list_Europe_large_no_clones.arg      
```

This resulted in 1’403’790  SNPs.

To generate the input files for fineStructure we need to know the per base recombination rates, these where [obtained](../recombination_map/recombination_map.md) starting from a genetic map and are stored in the file `../recombination_map/THUN12x96224_bp_recombination_rates.txt`

We generate the id and phase files for fineStructure, we also generate a .pos file:

```
zcat Europe_large_tritici_no_clones_no_miss.vcf.gz > Europe_large_tritici_no_clones_no_miss.vcf
python make_input_files_4_fs.py -vcf Europe_large_tritici_no_clones_no_miss.vcf -o Europe_large
```

Starting from the .pos file and the recombination map we generate the recombination file for fineStructure

```
python make_input_rec_file_4_fs.py -rec THUN12x96224_bp_recombination_rates.txt -i Europe_large.pos_file -o Europe_large
```


As fineStructure cannot deal with sample names starting with digits we rename these two isolates:
```
sed -i 's/96224/CHE_96224/g' Europe_large.id_file
sed -i 's/94202/CHE_94202/g' Europe_large.id_file
```

## fineStructure

We run fineStructure with default parameters, except that we increase the number of iterations in the EM algorithm to 50 (default 10)

```
fs Europe_large -idfile Europe_large.id_file -phasefiles Europe_large.hap_file -recombfiles Europe_large_cp_rec_file.txt -ploidy 1 -v -n -hpc 1 -s1args:-in\ -iM\ -i\ 50\ --emfilesonly -go
```
With the command above fineStructure generates lists of commands to run at different stages, these list can be submitted as batch jobs to a computer cluster.

When all stages are completed the three most important outputs are the chromopainter chunkcounts file (`Europe_large_linked_hap.chunkcounts.out`), the fineStructure mcmc file (`Europe_large_linked_hap_mcmc.xml`), and the fineStructure tree file (`Europe_large_linked_hap_tree.xml`).

## Plot results and PCA

This code is based on the example provided by authors of fineStructure. You need FinestructureLibrary.R and FinestructureDendrogram.R provided [here](https://people.maths.bris.ac.uk/~madjl/finestructure/toolsummary.html).

```
Rscripts plot_fs_results.R
```

This generates the [coacestry matrix plot](Europe_large_Coancestry.pdf) and the [PCA plot](Europe_large_PCA.pdf)

## Software versions:
```
gatk 4.4.0.0
fineStructure 4.4.4

- python and python modules

python 3.10.9
numpy 1.23.5    
argparser 1.4.0

- R and R packages
R 4.2.1
psych 2.3.12
GPArotation 2023.11-1
paran 1.5.2

```
