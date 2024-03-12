# fineStructure

For the fineStructure analyses we used the extended Europe dataset (link with XX individuals etc.).
We filtered all SNPS with any missing data, as fs cannot deal with them, and 

```
gatk-4.4.0.0/gatk SelectVariants \
     -R GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe_large_tritici_no_clones_no_miss_new.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals MT880591.1 \                    #exclude mitochondrion
     --exclude-intervals LR026995.1_Un \                 #exclude conting with sequenced non assined to any chromosome
     --exclude-intervals Bgt_MAT_1_1_3 \                 #exclude contig with alternative mating type
     --sample-name list_Europe_large_no_clones.args      
```



To generate the input files for fineStructure we need to know the per base recombination rates, these where (obtained)[../recombination_map/recombination_map.md] starting from a genetic map and are stored in the file `../recombination_map/THUN12x96224_bp_recombination_rates.txt`


```
zcat Europe_large_tritici_no_clones_no_miss.vcf.gz > Europe_large_tritici_no_clones_no_miss.vcf
python make_input_files_4_fs.py -vcf Europe_large_tritici_no_clones_no_miss.vcf -o Europe_large
python make_input_rec_file_4_fs.py -rec THUN12x96224_bp_recombination_rates.txt -i Europe_large.pos_file -o Europe_large
```


As fineStructure cannot deal with sample names starting with digit we rename these two isolates.
```
sed -i 's/96224/CHE_96224/g' Europe_large.id_file
sed -i 's/94202/CHE_94202/g' Europe_large.id_file
```

software versions:

gatk 4.4.0.0

python 3.10.9

numpy 1.23.5
