#!/bin/bash

#SBATCH --time=0-5:00:00
#SBATCH --mem 300M
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load anaconda3
source activate conda_env_gatk

# Filter and subset for iHS

# whole dataset
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/nikos/selection_scans/polarization/2_output/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps_polarized.vcf.gz \
-O /home/nminad/projects/nikos/selection_scans/iHS/0_data/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps_filtered.vcf.gz \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AF > 0.05"

# N_EUR
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/nikos/selection_scans/polarization/2_output/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps_polarized.vcf.gz \
-O /home/nminad/projects/nikos/selection_scans/iHS/0_data/tritici_recent_extended_europe_fs_level4_N_EUR.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_N_EUR.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AF > 0.05"

# S_EUR2
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/nikos/selection_scans/polarization/2_output/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps_polarized.vcf.gz \
-O /home/nminad/projects/nikos/selection_scans/iHS/0_data/tritici_recent_extended_europe_fs_level4_S_EUR2.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_S_EUR+.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AF > 0.05"

# TUR
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/nikos/selection_scans/polarization/2_output/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps_polarized.vcf.gz \
-O /home/nminad/projects/nikos/selection_scans/iHS/0_data/tritici_recent_extended_europe_fs_level4_TUR.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_TUR.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AF > 0.05"

# ME
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/nikos/selection_scans/polarization/2_output/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps_polarized.vcf.gz \
-O /home/nminad/projects/nikos/selection_scans/iHS/0_data/tritici_recent_extended_europe_fs_level4_ME.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_ME.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AF > 0.05"

# S_EUR1
/home/nminad/data/gatk-4.4.0.0/gatk SelectVariants -V /home/nminad/projects/nikos/selection_scans/polarization/2_output/tritici_2022+before2022+2023+ncsu_ALL_+outgroup_biallelic_snps_polarized.vcf.gz \
-O /home/nminad/projects/nikos/selection_scans/iHS/0_data/tritici_recent_extended_europe_fs_level4_S_EUR1.vcf.gz \
-sn /home/nminad/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_recent_ext_eur_fs_level4_S_EUR1.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AF > 0.05"
