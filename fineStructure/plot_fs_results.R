source("FinestructureLibrary.R") # read in the R functions, which also calls the needed packages
source("FinestructureDendrogram.R")

## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

### Define our input files
chunkfile<-"Europe_large_linked_hap.chunkcounts.out" ## chromopainter chunkcounts file
mcmcfile<-"Europe_large_linked_hap_mcmc.xml" ## finestructure mcmc file
treefile<-"Europe_large_linked_hap_tree.xml" ## finestructure tree file

###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 

###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame

###### READ IN THE TREE FILES

treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
## now is a good time to remove them via:
#     ttree$node.label<-NULL
## Will will instead remove "perfect" node labels
ttree$node.label[ttree$node.label=="1"] <-""
## And reduce the amount of significant digits printed:
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)

tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format

####################################

## PLOT 1: RAW DENDROGRAM PLOT




pdf(file="Europe_large_dendro.pdf",height=6,width=14)
par(mar=c(6,0,2,0),mfrow=c(1,1))
plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.3),edgePar=list(p.lwd=0,t.srt=90,t.off=-0.5),axes=F)
dev.off()

## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations

popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only


popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend

popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend

########################
## COANCESTRY MATRIX
fullorder<-labels(tdend)
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmatmax<-500 # cap the heatmap
tmpmat<-datamatrix 
tmpmat[tmpmat>tmatmax]<-tmatmax #
par(mar = c(0.1, 0.1, 0.1, 0.1))
pdf(file="Europe_large_Coancestry.pdf",height=12,width=12)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=0.6,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8))
dev.off()





########################

library(paran)
library(GPArotation)
library(psych)
library(patchwork)
library(Polychrome)
library(RColorBrewer)
## PCA Principal Components Analysis
pcares<-mypca(dataraw)
# For figuring out how many PCs are important; see Lawson & Falush 2012
# You need packages GPArotation and paran
#tmap<-optimalMap(dataraw)
#thorn<-optimalHorn(dataraw)
#c(tmap,thorn) # 11 and 5. Horn typically underestimates, Map is usually better
pcapops<-getPopIndices(rownames(dataraw),mapstatelist)
pcanames<-rownames(dataraw)


pca_df <- as.data.frame(pcares$vectors)
pca_colnames<- c()
for (x in 1:415) { 
  pca_colnames <- c(pca_colnames,paste0("PC",x))
}

names(pca_df) <- pca_colnames
pca_df$isolate <- pcanames


iso_det <- read.table("~/projects/vcf_project_tritici_old/2022+before2022+2023+ncsu_metadata+fs+admxK7_11032024.csv",header=TRUE,sep=",")
pca_dataset <- merge(pca_df,iso_det,by.x="isolate",by.y="Sample.Name")



library(ggplot2)
library(pals)
library(ggpubr)

pca_dataset$Collection <- as.factor(pca_dataset$Collection)
pca_dataset$Country <- as.factor(pca_dataset$Country)

var_frac <- pcares$values/sum(pcares$values)

## modified as.vector(polychrome(30)) to better distinguish between russia and turkey
my_pal <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B",
            "#C4451C", "#B10DA1", "#85660D", "#1C8356", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5")

# PC1 vs PC2
eur_pc12 <- ggplot(pca_dataset, aes(PC1, PC2, colour=Country,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(var_frac[1]*100, 3),"%", ")",  sep = ""), 
        y = paste("PC2", " ", "(", signif(var_frac[2]*100, 3), "%",")",  sep = ""),
        colour = "Country", shape = "Collection") + scale_color_manual(values = my_pal)

# PC1 vc PC3
eur_pc13 <- ggplot(pca_dataset, aes(PC1, PC3, colour=Country,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC1", " ", "(", signif(var_frac[1]*100, 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(var_frac[3]*100, 3), "%",")",  sep = ""),
        colour = "Country", shape = "Collection") + scale_color_manual(values = my_pal) 

#PC2 vs PC3
eur_pc23 <- ggplot(pca_dataset, aes(PC2, PC3, colour=Country,shape=Collection)) + geom_point(size=3) + theme_classic(base_size = 14)+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  labs (x = paste("PC2", " ", "(", signif(var_frac[2]*100, 3),"%", ")",  sep = ""), 
        y = paste("PC3", " ", "(", signif(var_frac[3]*100, 3), "%",")",  sep = ""),
        colour = "Country", shape = "Collection") + scale_color_manual(values = my_pal) 

eur_all <- eur_pc12 + eur_pc13 + eur_pc23 + guide_area() + plot_layout(guides = "collect") & guides(shape=guide_legend(order = 1, ncol = 2 ), colour = guide_legend(order = 2))

ggsave(eur_all, filename = "tritici_europe+_pca_fs_plots.pdf", height = 30, width = 45, unit = "cm")

ggsave(eur_all, filename = "tritici_europe+_pca_fs_plots.png", height = 30, width = 45, unit = "cm")










