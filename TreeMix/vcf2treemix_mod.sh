#!/bin/bash

# Script to convert vcf file to Treemix format
# by Joana Meier (https://github.com/speciationgenomics/scripts/blob/master/vcf2treemix.sh)
# requires vcftools, plink and plink2treemix.py from Treemix package to be installed
# takes as input the vcf file and the clust file
# modified 10052024 jigisha

module load anaconda3
source activate pixy_env

if [ $# -ne 3 ]
 then
 echo "Please provide the following arguments: <vcf file> <clust file> <fs_level(int)>"
 echo "The .clust file contains three columns: samplename\tsamplename\tgroup"
 exit 1
fi

clust=$2
file=${1%.gz}
file=${file%.vcf}
fs=$3

# Use VCFtools to make a map and a ped file, using only bi-allelic SNPs with mac 2 (also creates a log file)
if [ -s $file.vcf.gz ]
then

 # Get a .map and .ped file
 vcftools --gzvcf $file".vcf.gz" \
         --plink --mac 2 --remove-indels --max-alleles 2 \
         --out "fs${fs}/fs_${fs}_${file}"

else
 file=${file%.vcf}
 vcftools --vcf $file".vcf" \
         --plink --mac 2 --remove-indels --max-alleles 2  \
         --out "fs${fs}/fs_${fs}_${file}"

fi

# convert it to a stratified frq file, also creates .bed, .bim, .fam, .log, .nosex
~/data/plink --file "fs${fs}/fs_${fs}_${file}" --make-bed --out "fs${fs}/fs_${fs}_${file}" --allow-no-sex --allow-extra-chr 0
~/data/plink --bfile "fs${fs}/fs_${fs}_${file}" --freq --missing --within $2 --out "fs${fs}/fs_${fs}_${file}" --allow-no-sex --allow-extra-chr 0

# zip it
gzip "fs${fs}/fs_${fs}_${file}.frq.strat"

# create input file for Treemix
python2 plink2treemix.py "fs${fs}/fs_${fs}_${file}.frq.strat.gz" "fs${fs}/fs_${fs}_${file}.treemix.frq.gz"

