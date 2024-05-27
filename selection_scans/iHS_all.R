library(tidyverse)
library(rehh)
library(vcfR)


setwd("/home/nminad/projects/nikos/selection_scans/iHS/0_data/")


# Open and read the file list including all the variable names
f_files_1 <- list.files(full.names = T, pattern = ".vcf.gz")
f_files <- f_files_1[c(TRUE, FALSE)]
f_names <- list("ALL", "ME", "N_EUR", "S_EUR+", "S_EUR1", "TUR")

# Start the loop
for (i in 1:length(f_files)) {

  # calculate iHS
  
  # chromosomes
  fs_level2_ALL_hh_chr1 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026984.1_chr1", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr2 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026985.1_chr2", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr3 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026986.1_chr3", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr4 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026987.1_chr4", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr5 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026988.1_chr5", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr6 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026989.1_chr6", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr7 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026990.1_chr7", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr8 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026991.1_chr8", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr9 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026992.1_chr9", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr10 <- data2haplohh(hap_file = f_files[i],
                                         chr.name = "LR026993.1_chr10", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level2_ALL_hh_chr11 <- data2haplohh(hap_file = f_files[i],
                                         chr.name = "LR026994.1_chr11", polarize_vcf = TRUE, vcf_reader = "vcfR")
  
  
  ## perform scans
  fs_level2_ALL_hh_chr1_scan <- scan_hh(fs_level2_ALL_hh_chr1, polarized = TRUE)
  fs_level2_ALL_hh_chr2_scan <- scan_hh(fs_level2_ALL_hh_chr2, polarized = TRUE)
  fs_level2_ALL_hh_chr3_scan <- scan_hh(fs_level2_ALL_hh_chr3, polarized = TRUE)
  fs_level2_ALL_hh_chr4_scan <- scan_hh(fs_level2_ALL_hh_chr4, polarized = TRUE)
  fs_level2_ALL_hh_chr5_scan <- scan_hh(fs_level2_ALL_hh_chr5, polarized = TRUE)
  fs_level2_ALL_hh_chr6_scan <- scan_hh(fs_level2_ALL_hh_chr6, polarized = TRUE)
  fs_level2_ALL_hh_chr7_scan <- scan_hh(fs_level2_ALL_hh_chr7, polarized = TRUE)
  fs_level2_ALL_hh_chr8_scan <- scan_hh(fs_level2_ALL_hh_chr8, polarized = TRUE)
  fs_level2_ALL_hh_chr9_scan <- scan_hh(fs_level2_ALL_hh_chr9, polarized = TRUE)
  fs_level2_ALL_hh_chr10_scan <- scan_hh(fs_level2_ALL_hh_chr10, polarized = TRUE)
  fs_level2_ALL_hh_chr11_scan <- scan_hh(fs_level2_ALL_hh_chr11, polarized = TRUE)
  
  # iHS (withing population scans)
  
  ## Sabeti et al. 2007:
  ## integrated Haplotype Score (iHS) tests, rely on the principle that,
  ## under positive selection, an allele may rise to high frequency rapidly enough
  ## that long-range association with nearby polymorphisms the long-range haplotype will not have time
  ## to be eliminated by recombination. These tests control for local variation in recombination rates
  ## by comparing long haplotypes to other alleles at the same locus.
  ## As a result, they lose power as selected alleles approach fixation (100% frequency),
  ## because there are then few alternative alleles in the population
  fs_level2_ALL_hh_chr1_ihs <- ihh2ihs(fs_level2_ALL_hh_chr1_scan)
  fs_level2_ALL_hh_chr2_ihs <- ihh2ihs(fs_level2_ALL_hh_chr2_scan)
  fs_level2_ALL_hh_chr3_ihs <- ihh2ihs(fs_level2_ALL_hh_chr3_scan)
  fs_level2_ALL_hh_chr4_ihs <- ihh2ihs(fs_level2_ALL_hh_chr4_scan)
  fs_level2_ALL_hh_chr5_ihs <- ihh2ihs(fs_level2_ALL_hh_chr5_scan)
  fs_level2_ALL_hh_chr6_ihs <- ihh2ihs(fs_level2_ALL_hh_chr6_scan)
  fs_level2_ALL_hh_chr7_ihs <- ihh2ihs(fs_level2_ALL_hh_chr7_scan)
  fs_level2_ALL_hh_chr8_ihs <- ihh2ihs(fs_level2_ALL_hh_chr8_scan)
  fs_level2_ALL_hh_chr9_ihs <- ihh2ihs(fs_level2_ALL_hh_chr9_scan)
  fs_level2_ALL_hh_chr10_ihs <- ihh2ihs(fs_level2_ALL_hh_chr10_scan)
  fs_level2_ALL_hh_chr11_ihs <- ihh2ihs(fs_level2_ALL_hh_chr11_scan)
  
  fs_level2_ALL_ihs <- rbind(fs_level2_ALL_hh_chr1_ihs$ihs,
                             fs_level2_ALL_hh_chr2_ihs$ihs,
                             fs_level2_ALL_hh_chr3_ihs$ihs,
                             fs_level2_ALL_hh_chr4_ihs$ihs,
                             fs_level2_ALL_hh_chr5_ihs$ihs,
                             fs_level2_ALL_hh_chr6_ihs$ihs,
                             fs_level2_ALL_hh_chr7_ihs$ihs,
                             fs_level2_ALL_hh_chr8_ihs$ihs,
                             fs_level2_ALL_hh_chr9_ihs$ihs,
                             fs_level2_ALL_hh_chr10_ihs$ihs,
                             fs_level2_ALL_hh_chr11_ihs$ihs)
  
  # Write the results in .csv
  write.csv(fs_level2_ALL_ihs, paste("../2_output/iHS_fs_level2_", f_names[i], ".csv", sep = ""))
}
