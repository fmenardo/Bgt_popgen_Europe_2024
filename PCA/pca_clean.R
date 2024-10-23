library(vcfR)
library(adegenet)
library(argparser)
library(ggplot2)



parser <- arg_parser("run pca")

parser <- add_argument(parser,"-c", type="integer", default=1,help="number of cores")
parser <- add_argument(parser,"-s",flag = TRUE, help="save pca R object (needed for dapc)")
parser <- add_argument(parser,"-l",nargs=1,help="comma delimited table with samples names and metadata")
parser <- add_argument(parser,"-C",nargs= 1,type="integer",help="number of PC to be retained (as large as the n of samples)")
parser <- add_argument(parser,"-i", nargs= 1,help="input vcf file")
parser <- add_argument(parser,"-o", nargs= 1,help="output stem")
                                                                                  
                                          
args <- parse_args(parser)       

vcf <- read.vcfR(args$i)

genl <- vcfR2genlight(vcf, n.cores = args$c)


pca<-glPca(genl, center = TRUE, scale = FALSE, loadings = FALSE,nf=args$C, 
      alleleAsUnit = FALSE, useC = FALSE, parallel = TRUE,
      n.cores = args$c, returnDotProd=FALSE, matDotProd=NULL )

iso_det <- read.table(args$l,header=TRUE,sep=",")

var_frac <- pca$eig/sum(pca$eig)
pca_df <- as.data.frame(pca$scores)

pca_df$isolate <- rownames(pca_df)

pca_dataset <- merge(pca_df,iso_det,by.x="isolate",by.y="Sample.ID")
pca_dataset$Country <- as.factor(pca_dataset$Country)
#pca_dataset$Collection <- as.factor(pca_dataset$Collection)
pca_dataset$Strain <- as.factor(pca_dataset$Strain)


write.csv(pca_dataset,paste0(args$o,".pca.csv"),row.names = FALSE)
write.csv(as.data.frame(pca$eig),paste0(args$o,".eig.csv"),row.names=FALSE)

if (args$s==TRUE){
	save(pca,file=paste0(args$o,".glpca"))
	}

