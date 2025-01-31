#!/bin/bash

#SBATCH --time=0-24:00:00
#SBATCH --mem 10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --job-name=LD_decay_plink2_10kb_100kb

# gatk filter and subset
module load anaconda3
source activate conda_env_gatk

#change the args and names

# whole dataset
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
-O /home/nminad/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1.vcf.gz \
-sn ~/projects/vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1

# N_EUR
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz \
-O /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_N_EUR.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_N_EUR.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 

# S_EUR+
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz \
-O /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR+.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_S_EUR+.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1

# TUR
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz \
-O /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_TUR.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_TUR.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1

# ME
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz \
-O /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_ME.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_ME.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1

# S_EUR1
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps.vcf.gz \
-O /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR1.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_S_EUR1.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1

# vcf to plink format
# whole dataset
/home/nminad/data/plink2 --vcf /home/nminad/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1.vcf.gz\
--make-bed --out /home/nminad/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1

# N_EUR
/home/nminad/data/plink2 --vcf /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_N_EUR.vcf.gz \
--make-bed --out /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_N_EUR

# S_EUR+
/home/nminad/data/plink2 --vcf /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR+.vcf.gz \
--make-bed --out /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR+

# TUR
/home/nminad/data/plink2 --vcf /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_TUR.vcf.gz \
--make-bed --out /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_TUR

# ME
/home/nminad/data/plink2 --vcf /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_ME.vcf.gz \
--make-bed --out /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_ME

# S_EUR1
/home/nminad/data/plink2 --vcf /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR1.vcf.gz \
--make-bed --out /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR1

# run LD decay ld-window 10 kb

# whole dataset
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1 \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_ext_eur_recent_10kb
# N_EUR
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_N_EUR \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_N_EUR_10kb

# S_EUR+
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR+ \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_S_EUR+_10kb

# TUR
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_TUR \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_TUR_10kb

# ME
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_ME \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_ME_10kb

# S_EUR1
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR1 \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_S_EUR1_10kb

# run LD decay ld-window 100 kb

# whole dataset
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps_filtered_miss0.1 \
--r2-unphased --ld-window-r2 0 --thin 0.2 --ld-window-kb 100 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_ALL_100kb

# N_EUR
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_N_EUR \
--r2-unphased --ld-window-r2 0 --thin 0.2 --ld-window-kb 100 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_N_EUR_100kb

# S_EUR+
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR+ \
--r2-unphased --ld-window-r2 0 --thin 0.2 --ld-window-kb 100 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_S_EUR+_100kb

# TUR
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_TUR \
--r2-unphased --ld-window-r2 0 --thin 0.2 --ld-window-kb 100 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_TUR_100kb

# ME
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_ME \
--r2-unphased --ld-window-r2 0 --thin 0.2 --ld-window-kb 100 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_ME_100kb

# S_EUR1
/home/nminad/data/plink2 --bfile /home/nminad/projects/nikos/LD_decay/0_data/tritici_recent_extended_europe_fs_level4_S_EUR1 \
--r2-unphased --ld-window-r2 0 --thin 0.2 --ld-window-kb 100 --out /home/nminad/projects/nikos/LD_decay/2_output/LD_plink2_S_EUR1_100kb







