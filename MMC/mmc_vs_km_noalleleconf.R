#' Script to perform model selection and parameter estimation 
#' for populations N_EUR2 and E_EUR2 for 
#' Beta vs Kingman+exp. growth without
#' accounting effects of confusing the ancestral and derived 
#' alleles
#' To run the script, you need the git repository
#' https://github.com/fabfreund/usfs_mmc/
#' Data should be provided as a Fasta file, with the ancestral sequence 
#' in row 1
index <- 2

#' Needs to be adjusted when rerun
data_path <- switch(index,
             "data/N_EUR2_all_chr_BISNP.fasta",
             "data/E_EUR2_all_chr_BISNP.fasta")


short_name  <- switch(index,
                      "N_EUR2","E_EUR2")

#' Switches for calculations
#' 1) should the data be read in and the SFS computed?
#' 2) Should the table of expected branch lengths (Phi table)
#' be computed?
#' 3) Should the maximum likelihood inference be performed?
#' 4) Should the expected (scaled) SFS for the 
#' best fitting KM and Beta model and the observed SFS be plotted 
#' (as in Fig. S12)
#' To run any step, all previous steps have to be run.
readin <- FALSE
phi <- FALSE
ml_comp <- FALSE
plotit <- TRUE

if (readin){
#' Observed data
#' Should contain the ancestral sequence in row 1
library(pegas)
#' If necessary:
#' Import data and get site frequency spectrum
#' pegas does not recognize * in Fasta - change to -
com1 <- paste("cat",data_path)
com1 <- paste(com1,"| tr '*' 'n'  >")
com1 <- paste0(com1," powmil_",short_name,".fasta")
system(com1)


data_powmild <- read.dna(paste0("powmil_",short_name,".fasta"),
                         format = "fasta")
data_powmild_char <- as.character(data_powmild)
#' check overall base frequencies
base.freq(data_powmild,all = TRUE)
#' get columns with any non-ACTG data 
cols_wo_fulldata <- apply(data_powmild_char,2,
                          function(v){bool_anc <- (v[1]%in% c("a","t","g","c"))
                          uv <- unique(v[-1])
                          bool1 <- all(uv %in% c("a","t","g","c"))
                          bool1 <- bool1 & length(uv)==2
                          bool1 <- bool1 & bool_anc
                          return(!bool1)})

#Remove all columns with multiallelic SNPs or missing data (including gaps)
data_powmild_char <- data_powmild_char[,!cols_wo_fulldata]

#get count per site
count_per_site <- function(v){sum(v[-1]!=v[1])}

sfs_obs <- apply(data_powmild_char,2,count_per_site)

#' Get SFS from counts
sfs_as_sfs <- function(countv,n=nrow(data_powmild_char)-1){
  out1 <- rep(-1,n-1)
  tab1 <- table(countv)
  for (i in 1:length(out1)){
    out1[i] <- unname(tab1[as.character(i)])}
  if (any(is.na(out1))){out1[is.na(out1)] <- 0}
  return(out1)}

sfs_obs <- sfs_as_sfs(sfs_obs)
#' Write out sfs
write.table(t(sfs_obs),
           file = paste0("powmil_",short_name,".sfs"),
            row.names = FALSE,col.names = FALSE,
            quote=FALSE,sep = "\t",eol = "\n")
}
if (phi){
sfs <- read.delim(paste0("powmil_",
                         short_name,".sfs"),
                  header = FALSE)
#' Get sample size from sfs
n <- length(sfs) + 1 

#' For the following two sections, we need code from
#' the git repository https://github.com/fabfreund/usfs_mmc/
#'
#' Compute Phi table - NEEDS path adjustment when rerun
com_phi <- paste0("./../../../git_repos/",
                  "usfs_mmc/main_sim_inf_tool/Beta/",
                  "MMC-Phi-LookUpTable/",
                  "MMC-Phi-LookUpTable.out", 
                  " -sampleSize ",n,
                  " -MaxRho 10 -NoStepsRho 40",
                  " -OUT powmil_",short_name,".phi")
system(com_phi)
system(paste0("mv powmil_",short_name,".phi* powmil_",short_name,".phi"))
}
#' Run Max-Likelihood comp
if (ml_comp){

#' Path adjustment necessary - path_to_input is the file where the sfs
#' and the phi-table are situated (which for us were situated in the folder
#' of this script). Also necessary to adjust the path to the executable
#' from the U-SFS git repository
path_to_input <- getwd()  
com_ml <- paste0("./../../../git_repos/usfs_mmc/",
                 "main_sim_inf_tool/Beta/",
                 "MMC-MaxLikelihoodInference-",
                 "GridSearch/MMC-MaxLikelihood",
                 "Inference-GridSearch.out", 
                 " -minRho 0 -maxRho 10",
                 " -noStepsRho 40 -minAlpha 1", 
                 " -maxAlpha 2 -noStepsAlpha 100",
                 " -printGrid -SFS ", 
                 path_to_input,
                 "/powmil_",
                 short_name,".sfs",
                 " -Phi ",path_to_input,"/",
                 "powmil_",short_name,".phi")   
system(com_ml)
}


#' Plot results
if (plotit){

library(tidyverse)
library(hrbrthemes)
library(viridis)
#' The ML grid computation produces several output files, 
#' we want the one whose name includes "Grid.txt"
#' In the case we have rerun at different timepoints,
#' we may even have several ones, so this is why we force
#' a specific entry with [1] 
grid_name <- list.files(pattern=paste0("powmil_",
                      short_name,"-L*"))

grid_name <- grid_name[grep(x = grid_name,
                   pattern = "Grid.txt")]  
grid1 <- read_table(grid_name[1],
                      skip = 19,col_names = TRUE)  
sfs <- read_table(paste0("powmil_",short_name,
                         ".sfs"),
                  col_names =FALSE)
sfs <- t(sfs)
cat("No of SNPS",sum(sfs[,1]))
sfs <- unname(sfs[,1])/sum(sfs[,1])

phi_tab_alpha <- read_table(file=paste0("powmil_",
                                        short_name,
                                        ".phi"),
                            col_names = FALSE)

#' Extract likelihood ratio between best fitting KM and best fitting Beta
maxpos_km <- which(grid1$LogL==max(grid1$LogL[abs(grid1$Alpha-2)<0.00001]))
maxpos_mmc <- which.max(grid1$LogL)
logratio_test <- grid1$LogL[maxpos_mmc] - grid1$LogL[maxpos_km]

#maxlik_km <- max(grid1$LogL[abs(grid1$Alpha-2)<0.00001])
#maxlik_mmc <- max(grid1$LogL)
#logratio <- maxlik_mmc - maxlik_km

sfs_beta <- unname(unlist(phi_tab_alpha[maxpos_mmc,-1]))
sfs_km <- unname(unlist(phi_tab_alpha[maxpos_km,-1]))

sfs <- tibble(observed=sfs,beta=sfs_beta,km=sfs_km)


plot_data <- tibble("SFS classes"=rep(seq(along=sfs$observed),3),
                    "Fraction of mutations"=c(sfs$observed,
                                              sfs$beta,sfs$km),
                    "type"=factor(c(rep("observed",nrow(sfs)),
                                    rep("best Beta",nrow(sfs)),
                                    rep("best Kingman",nrow(sfs)))))

#' To save the data for the plot, optional
saveRDS(sfs,file = paste0("plotdata_",short_name,".RDS"))
#' Plots
pdf(paste0("plots_",short_name,".pdf"))

  
  types <- c("Observed"="black","Beta"="red","KM"="blue")
  
plot1 <-  ggplot(data=sfs) + geom_point(aes(y=observed,
                                    x=seq(along = observed),
                                          colour="Observed")) +
    geom_line(aes(y=beta,x=seq(along = observed),colour = "Beta")) +
    geom_line(aes(y=km,x=seq(along = observed),colour="KM")) +
    labs(x = "SFS class", y= "Fraction") +
    labs(title=paste0("SFS ",short_name)) +
    scale_colour_manual(values=types) +
    theme_classic() +
    theme(legend.text = element_text(size = 11),
          legend.background = element_rect(linewidth = .1,
                                           colour = "black" 
                                            ),
          legend.title=element_blank(),
          legend.position="inside",
          legend.justification = c(1,1))

print(plot1)  

dev.off()
#' Return the log likelihood ratio between the 
#' best fitting Kingman+exp. growth model and the best 
#' fitting Beta coalescent
print(logratio)
}

