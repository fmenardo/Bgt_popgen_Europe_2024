library(tidyverse)
library(rehh)
library(vcfR)

setwd("/home/nminad/projects/nikos/selection_scans/xp-EHH/0_data/")


### calculate xp-EHH

## read in data for each clade
# create a list for the 5 populations
f_files_1 <- list.files(full.names = T, pattern = ".vcf.gz")
f_files <- f_files_1[c(TRUE, FALSE)]
f_names_1 <- list.files(full.names = F, pattern = ".vcf.gz")
f_names <- f_names_1[c(TRUE, FALSE)]
#f_names <- list("ALL", "ME", "N_EUR", "S_EUR2", "S_EUR1", "TUR")

#population_listoflists <- list()

# Start the loop
for (i in 1:length(f_files)) {
  
  # calculate iHS
  
  # chromosomes
  fs_level4_ALL_hh_chr1 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026984.1_chr1", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr2 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026985.1_chr2", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr3 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026986.1_chr3", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr4 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026987.1_chr4", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr5 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026988.1_chr5", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr6 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026989.1_chr6", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr7 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026990.1_chr7", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr8 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026991.1_chr8", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr9 <- data2haplohh(hap_file = f_files[i],
                                        chr.name = "LR026992.1_chr9", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr10 <- data2haplohh(hap_file = f_files[i],
                                         chr.name = "LR026993.1_chr10", polarize_vcf = TRUE, vcf_reader = "vcfR")
  fs_level4_ALL_hh_chr11 <- data2haplohh(hap_file = f_files[i],
                                         chr.name = "LR026994.1_chr11", polarize_vcf = TRUE, vcf_reader = "vcfR")
  
  
  ## perform scans
  fs_level4_ALL_hh_chr1_scan <- scan_hh(fs_level4_ALL_hh_chr1, polarized = TRUE)
  fs_level4_ALL_hh_chr2_scan <- scan_hh(fs_level4_ALL_hh_chr2, polarized = TRUE)
  fs_level4_ALL_hh_chr3_scan <- scan_hh(fs_level4_ALL_hh_chr3, polarized = TRUE)
  fs_level4_ALL_hh_chr4_scan <- scan_hh(fs_level4_ALL_hh_chr4, polarized = TRUE)
  fs_level4_ALL_hh_chr5_scan <- scan_hh(fs_level4_ALL_hh_chr5, polarized = TRUE)
  fs_level4_ALL_hh_chr6_scan <- scan_hh(fs_level4_ALL_hh_chr6, polarized = TRUE)
  fs_level4_ALL_hh_chr7_scan <- scan_hh(fs_level4_ALL_hh_chr7, polarized = TRUE)
  fs_level4_ALL_hh_chr8_scan <- scan_hh(fs_level4_ALL_hh_chr8, polarized = TRUE)
  fs_level4_ALL_hh_chr9_scan <- scan_hh(fs_level4_ALL_hh_chr9, polarized = TRUE)
  fs_level4_ALL_hh_chr10_scan <- scan_hh(fs_level4_ALL_hh_chr10, polarized = TRUE)
  fs_level4_ALL_hh_chr11_scan <- scan_hh(fs_level4_ALL_hh_chr11, polarized = TRUE)
  
  ## create a list and append to the population list of lists
  chromosome_scans_list <- list(fs_level4_ALL_hh_chr1_scan,
                          fs_level4_ALL_hh_chr2_scan,
                          fs_level4_ALL_hh_chr3_scan,
                          fs_level4_ALL_hh_chr4_scan,
                          fs_level4_ALL_hh_chr5_scan,
                          fs_level4_ALL_hh_chr6_scan,
                          fs_level4_ALL_hh_chr7_scan,
                          fs_level4_ALL_hh_chr8_scan,
                          fs_level4_ALL_hh_chr9_scan,
                          fs_level4_ALL_hh_chr10_scan,
                          fs_level4_ALL_hh_chr11_scan,
                          f_files[i])
  saveRDS(chromosome_scans_list, paste0("../../xp-EHH/2_output/population_",f_names[i],"_scans.RDS"))
  if (i==1){
    population_listoflists <- list(chromosome_scans_list)
  } else {
    population_listoflists <- append(population_listoflists, list(chromosome_scans_list))
  }
}

saveRDS(population_listoflists, "../../xp-EHH/2_output/population_list_scans.RDS")

