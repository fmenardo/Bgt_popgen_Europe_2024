library(vcfR)
library(adegenet)
library(ape)
library(ggplot2)
library(pals)
library(dplyr)
library(parallel)
library(argparser)


parser <- arg_parser("run dist matrix")

parser <- add_argument(parser,"-c", type="integer", default=1,help="number of cores")
parser <- add_argument(parser,"-i", nargs= 1,help="input vcf file")
parser <- add_argument(parser,"-o", nargs= 1,help="output stem")                                                                               
                                          
args <- parse_args(parser)       

vcf <- read.vcfR(args$i)

genl <- vcfR2genlight(vcf, n.cores = args$c)

dist_matr<-dist.gene(as.matrix(genl),method = "pairwise",pairwise.deletion=TRUE,variance=FALSE)
write.csv(as.matrix(dist_matr), paste0(args$o,"_2022+before2022+2023+ncsu_dist_mat.csv"))

jpeg(file=paste0(args$o,"_2022+before2022+2023+ncsu_dist_mat.jpeg"))
hist(dist_matr,col="steelblue")
dev.off()
