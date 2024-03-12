# fineStructure

For the fineStructure analyses we used the extended Europe dataset (link with XX individuals etc.).
fineStructure uses dense 


```
gatk-4.4.0.0/gatk SelectVariants \
    -R ~/projects/project_tritici_fabrizio/data/GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V ~/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe_large_tritici_no_clones_no_miss_new.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals MT880591.1 \
     --exclude-intervals LR026995.1_Un \
     --exclude-intervals Bgt_MAT_1_1_3 \
     --sample-name list_Europe_large_no_clones.args
```



cp ~/projects/project_recombination_rate_map/analysis/THUN12x96224_bp_recombination_rates.txt ~/projects/project_tritici_fabrizio/data/.

source activate biopython_env

srun zcat Europe_large_tritici_no_clones_no_miss.vcf.gz > Europe_large_tritici_no_clones_no_miss.vcf
srun python ../2022_2023/make_input_files_4_fs.py -vcf Europe_large_tritici_no_clones_no_miss.vcf -o Europe_large
srun python ../2022_2023/make_input_rec_file_4_fs.py -rec ~/projects/project_tritici_fabrizio/data/THUN12x96224_bp_recombination_rates.txt -i Europe_large.pos_file -o Europe_large

sed -i 's/96224/CHE_96224/g' Europe_large.id_file
sed -i 's/94202/CHE_94202/g' Europe_large.id_file

