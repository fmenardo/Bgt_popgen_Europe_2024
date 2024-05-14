#install.packages("devtools",lib = "~/data/R_lib")
#library(devtools,lib = "~/data/R_lib")

#devtools::install_github("bahlolab/isoRelate",lib = "~/data/R_lib")

library(isoRelate,lib = "~/data/R_lib")
#library(doParallel,lib = "~/data/R_lib")
library(argparser, lib = "~/data/R_lib")

parser <- arg_parser("run isoRelate")

#parser <- add_argument(parser,"-i", nargs= 1,help="input file stem: output stem of step 1")
parser <- add_argument(parser,"-o", nargs= 1,help="output stem of step 1")

args <- parse_args(parser)


load(paste0(args$o,"_geno.RData"))

# generate a binary IBD matrix
my_matrix <- getIBDmatrix(ped.genotypes = my_genotypes, 
                          ibd.segments = my_ibd)


# calculate the proportion of pairs IBD at each SNP
my_proportion <- getIBDproportion(ped.genotypes = my_genotypes, 
                                  ibd.matrix = my_matrix, 
                                  groups = NULL)

write.table(my_proportion, file = paste0(args$o,"_proportion_table.txt"),row.names = TRUE,col.names = TRUE)




# plot the proportion of pairs IBD
pdf(paste0(args$o,"plot_ibd_prop.pdf"))


plotIBDproportions(ibd.proportions = my_proportion, 
                   interval = NULL, 
                   annotation.genes = NULL,
                   annotation.genes.color = NULL,
                   highlight.genes = NULL,
                   highlight.genes.labels = TRUE,
                   highlight.genes.color = NULL,
                   highlight.genes.alpha = 0.1,
                   add.rug = FALSE, 
                   plot.title = "Proportion of pairs IBD", 
                   add.legend = FALSE,
                   line.color = NULL, 
                   facet.label = TRUE, 
                   facet.scales = "fixed", 
                   subpop.facet = FALSE)

	dev.off()
save(my_matrix,my_proportion, file = paste0(args$o,"_matrix.RData"))



# calculate the significance of IBD sharing
my_iR <- getIBDiR(ped.genotypes = my_genotypes, 
                  ibd.matrix = my_matrix, 
                  groups = NULL)

iR <- my_iR

for(i in 1:nrow(my_iR)) {   
  iR$new_p <- -log10(pchisq(my_iR$iR, df = 1, log.p=FALSE, lower.tail = FALSE))
}


write.table(iR, file = paste0(args$o,"_iR_table.txt"),row.names = TRUE,col.names = TRUE)


# plot the iR statistics
pdf(paste0(args$o,"plot_ibd_IR.pdf"))


plotIBDiR(ibd.iR = my_iR, 
          interval = NULL, 
          annotation.genes = NULL,
          annotation.genes.color = NULL,
          highlight.genes = NULL,
          highlight.genes.labels = FALSE,
          highlight.genes.color = NULL,
          highlight.genes.alpha = 0.1,
          point.size = 1,
          point.color = NULL,
          add.rug = FALSE, 
          plot.title = "Significance of IBD sharing", 
          add.legend = FALSE,
          facet.label = TRUE, 
          facet.scales = "fixed")
dev.off()



save(my_iR, file = paste0(args$o,"_iR.RData"))
