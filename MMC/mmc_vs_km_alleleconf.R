#' Script to perform model selection and parameter estimation 
#' for population E_EUR2 for Beta vs Kingman+exp. growth
#' including accounting for potential switching of 
#' ancestral and derived alleles (2nd version from Freund et al,
#' where one accounts for difference in occurence of transition and
#' transversion mutations)

#' Assumption: Fasta files including outgroup and 
#' positions that are multiallelic between outgroup and ingroup 
#' is in subfolder data/

#' This script uses the expected coalescent branch lengths tables
#' (phi tables) that are generated when running 
#' mmc_vs_km_noalleleconf.R with option phi <- TRUE in l. 14
#' for the E_EUR2 data set w/o outgroup sequences
#' This other script uses the code and computation techniques
#' from the git repository https://github.com/fabfreund/usfs_mmc/
 
#' This script sources two additional R scripts 
#' funct_and_consts.R and model_select.R,which were in a 
#' subdirectory aux_tstv. Adjust the path if necessary. 

#' Libraries
library(adegenet)
library(pegas)
library(future.apply)

#' Source filter scripts - to filter fasta files for SNP counts
#' You may need to adjust the path
source("aux_tstv/funct_and_consts.R")

nsamp <- 17 #sample size is 17
sfs_name <- "E_EUR2_allchrom"

#' function: Make SFS from data table
sfs_as_sfs <- function(tab1){out1 <- rep(-1,nsamp-1)
for (i in 1:length(out1)){
  out1[i] <- unname(tab1[as.character(i)])}
if (any(is.na(out1))){out1[is.na(out1)] <- 0}
return(out1)}

SFS_ts <- rep(0,16)
SFS_tv <- rep(0,16)
Sneq_ts <- 0
Sneq_tv <- 0
for (c1 in c(1,2,4,5,6,7,8,9,10)){
#' We read in data chromosome per chromosome  
  set_fasta <- paste0("data/LRX_chr",c1,
                      ".E_EUR2_MULTIALLELIC.genomic.fa")
  #' Whether to run also Phi estimation
  dophi <- FALSE
  
#' In case data is not fully formatted but still has * in it  
  if (TRUE){
    com1 <- paste("cat",set_fasta)
    com1 <- paste(com1,"| tr '*' 'n'  >")
    com1 <- paste0(com1," data/powmil_",sfs_name,".fasta")
    system(com1)
  }
  
  #' Read in data
  data1 <- read.dna(file=paste0("data/powmil_",sfs_name,
                                ".fasta"),format="fasta")
  
  #' Extract outgroup (information on outgroup lines manually provided)
  out_rows <- 16:20 #which rows/lines in the fasta file are from outgroup?
  print(rownames(data1)[out_rows])
  print(rownames(data1)[-out_rows])
  
  #' first filter step
  keep_sit1 <- apply(as.character(data1),2,ccheck_awful)
  data2 <- data1[,keep_sit1]
  rm(data1)
  gc()
  
  #' Count bases - which ones can be oriented etc
  res1 <- apply(as.character(data2),2,ccount_bases)
  
  
  
  
    
#' Extract SFS from transitions and SFS from transversions
#' and add to the ones from the previous chromosomes  
SFS_ts <- SFS_ts + sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1" & res1["ts",]=="1"]))
SFS_tv <- SFS_tv + sfs_as_sfs(table(res1["nmut",res1["iclean",]=="1" & res1["ts",]=="0"]))


#' Extract non-orientable sites (3rd base in outgroup) - both for in-sample ts and tv
#' add to the ones from previous chromosomes
Sneq_ts <- Sneq_ts + sum(res1["ts",]=="1" & res1["pol",]=="0" & res1["iclean",]==1,
                    na.rm = TRUE)
Sneq_tv <- Sneq_tv + sum(res1["ts",]=="0" & res1["pol",]=="0" & res1["iclean",]==1,
                    na.rm = TRUE)
}
#' After collating the SFS across chromosome, run the MMC vs KM inference
SFS <- list(ts=SFS_ts,tv=SFS_tv)
Sneq <- list(ts=Sneq_ts,tv=Sneq_tv)
#' Run model selection and ML parameter optimisation
source("aux_tstv/model_select.R")
#' Save output  
save(res_list,file=paste0("tvts_E_EUR_allchrom_ma.RData"))
