#install.packages("devtools",lib = "~/data/R_lib")
#library(devtools,lib = "~/data/R_lib")

#devtools::install_github("bahlolab/isoRelate",lib = "~/data/R_lib")

library(isoRelate,lib = "~/data/R_lib")
library(doParallel,lib = "~/data/R_lib")
library(argparser, lib = "~/data/R_lib")

parser <- arg_parser("run isoRelate")

parser <- add_argument(parser,"-p", nargs= 1,help="ped file")
parser <- add_argument(parser,"-m", nargs= 1,help="map file")
parser <- add_argument(parser,"-C", nargs= 1,help="min length in centimorgan for ibd segments to be included")
#parser <- add_argument(parser,"-M", nargs= 1,help="min allele frequenyE for a SNP to be included")
parser <- add_argument(parser,"-o", nargs= 1,help="outüut stem")
parser <- add_argument(parser,"-c", nargs= 1,help=" n cores")
#parser <- add_argument(parser,"-i", nargs= 1, default=0, help="if = 0 (default) performs complete analysis, if = 1 start from matrix calculation. i=2 analysis start after calculation of pairwise matrix ")



args <- parse_args(parser) 


ped <- read.table(file=args$p, header = FALSE, sep = "")


map <- read.table(file=args$m, header = FALSE, sep = "")


ped_map<- list(ped,map)


my_genotypes <- getGenotypes(ped.map = ped_map,
                             reference.ped.map = NULL,
#                            maf = as.numeric(args$M),
                             isolate.max.missing = 0,
                             snp.max.missing = 0,
                             chromosomes = NULL,
                             input.map.distance = "cM")
                             
                             
my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, number.cores = as.integer(args$c))


ibd <- getIBDsegments(ped.genotypes = my_genotypes,
                         parameters = my_parameters, 
                         number.cores = as.integer(args$c), 
                         minimum.snps = 50, 
                         minimum.length.bp = 50000,
                         error = 0.001)

M_threshold= as.integer(args$C)/100
	
my_ibd <- ibd[ibd$length_M > M_threshold,]

save(my_genotypes,my_ibd,my_parameters, file = paste0(args$o,"_geno.RData"))
	write.table(my_ibd, file = paste0(args$o,"_ibd_table.txt"),row.names = TRUE,col.names = TRUE)


pdf(paste0(args$o,"plot_ibd.pdf"))
plotIBDsegments(ped.genotypes = my_genotypes, 
                ibd.segments = my_ibd, 
                interval = NULL,
                annotation.genes = NULL,
                annotation.genes.color = NULL,
                highlight.genes = NULL, 
                highlight.genes.labels = FALSE,
                highlight.genes.color = NULL,
                highlight.genes.alpha = 0.1,
                segment.height = 0.6,
                number.per.page = 100, 
                fid.label = FALSE, 
                iid.label = FALSE, 
                ylabel.size = 9, 
                add.rug = FALSE,
                plot.title = "Distribution of IBD segments", 
                add.legend = TRUE, 
                segment.color = NULL)

dev.off()


