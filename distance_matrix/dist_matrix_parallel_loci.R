library(vcfR)
library(adegenet)
library(ape)
library(argparser)

parser <- arg_parser("run dist matrix")

parser <- add_argument(parser,"-c", type="integer", default=1,help="number of cores")
parser <- add_argument(parser,"-i", nargs= 1,help="input vcf file")
parser <- add_argument(parser,"-o", nargs= 1,help="output stem")                                                                               
                                          
args <- parse_args(parser)       

vcf <- read.vcfR(args$i)

genl <- vcfR2genlight(vcf, n.cores = args$c)

dist_matr<-dist.gene(as.matrix(genl),method = "pairwise",pairwise.deletion=TRUE,variance=TRUE)
dm <- as.matrix(dist_matr)
var <- as.matrix(attr(dist_matr, "variance"))
loci <- (dm^2)/(dm-var)

write.csv(dm, paste0(args$o,"_dist_mat_num_diffs.csv"))
write.csv(loci,paste0(args$o,"_dist_mat_num_loci.csv"))
