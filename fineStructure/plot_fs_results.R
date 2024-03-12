source("FinestructureLibrary.R") # read in the R functions, which also calls the needed packages
source("FinestructureDendrogram.R")

## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

### Define our input files
chunkfile<-"Europe_large_linked_hap.chunkcounts.out" ## chromopainter chunkcounts file
mcmcfile<-"Europe_large_linked_hap_mcmc.xml" ## finestructure mcmc file
treefile<-"Europe_large_linked_hap_tree.xml" ## finestructure tree file

## Additional files that you can extract from finestructure
#mappopchunkfile<-"EastAsiaSimple.EMlinked.mapstate.csv" # population-by-population chunkcount file for the populations used in the MAP (i.e tree)
#system( paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile) )
#meancoincidencefile<-"EastAsiaSimple.EMlinked.meancoincidence.csv" # pairwise coincidence, .i.e. proportion of MCMC files where individuals are found in the same 
#system( paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile) )
## there are ways of generating these within R but are either slower or more annoying - its your call how you do it

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
## PCA Principal Components Analysis
pcares<-mypca(dataraw)
# For figuring out how many PCs are important; see Lawson & Falush 2012
# You need packages GPArotation and paran
tmap<-optimalMap(dataraw)
thorn<-optimalHorn(dataraw)
c(tmap,thorn) # 11 and 5. Horn typically underestimates, Map is usually better
pcapops<-getPopIndices(rownames(dataraw),mapstatelist)
pcanames<-rownames(dataraw)
rcols<-rainbow(max(pcapops))

pdf("Europe_large_PCA.pdf",height=16,width=12)
par(mfrow=c(3,2))
for(i in 1:3) for(j in (i+1):4) {
  plot(pcares$vectors[,i],pcares$vectors[,j],col=rcols[pcapops],xlab=paste("PC",i),ylab=paste("PC",j),main=paste("PC",i,"vs",j),pch=rcols)
  text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames,col=rcols[pcapops],cex=0.5,pos=1)
}
dev.off()
