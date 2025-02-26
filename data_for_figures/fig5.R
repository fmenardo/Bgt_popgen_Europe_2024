## Fig 5 ##
#### Fig 5a,b ####
# haplotype networks
library(adegenet)
library(pegas)
library(seqinr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

load("../AvrPm17_haplotypes/hap_net_xy.Rdata")  # pre cooked coordinates for hap net plotting, can be commented out, plot will be less pretty
ali<-fasta2DNAbin("../AvrPm17_haplotypes/avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa")

### fetch list of isolates in ali
system("grep \">\" ../AvrPm17_haplotypes/avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa | sed \'s/>//g\'  > list_samples_for_h")
names<-read.table("list_samples_for_h",header=FALSE)

temp <- gsub("_1_avrpm17_cds","",names[,1])
names[,2] <- gsub("_2_avrpm17_cds","",temp)

## read metainfo
meta<-read.csv("../Datasets/S1_Data.csv")
meta1<-data.frame(meta$Sample_ID,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$Year.of.Collection)
colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Year")

## assign pop to h
colnames(names)<-c("Hap.Name","Sample.Name")

for (i in 1:nrow(names)){
  for (z in 1:nrow(meta1)){
    if (names$Sample.Name[i] == meta1$Sample.Name[z]){names$pop[i] = meta1$Population[z]}
  }
}


##assign varuiant nomenclature of Kunz et al, to haplotype


#I = B
#II = C
#III = A
#IV = A
#V = A
#VI = C
#VII = C
#VIII = H
#IX = B
#X = A
#XI = C
#XII = I
#XIII = J
#XIV = K
#XV = F
#XVI = A
#XVII = L
#XVIII = C
#XIX = C
#XX = B
#XXI = B


dna_variants<-read.csv("../AvrPm17_haplotypes/hap_class_cds_RC_no_miss.csv")
colnames(dna_variants) <- c("dna_variants","Sample.Name")

mapping <- c('0' = 'B', '1' = 'C', '2' = 'A','3'= 'A','4'= 'A','5'= 'C','6'= 'C','7'= 'H','8'= 'B','9'= 'A',
             '10'= 'C', '11'= 'I','12'= 'J','13'= 'K','14'= 'F','15'= 'A','16'= 'L','17'= 'C','18'= 'C',
             '19'= 'B','20'= 'B')

dna_variants$Variants <- mapping[as.character(dna_variants$dna_variants)]
names_var<- merge(names,dna_variants,by.x="Hap.Name",by.y="Sample.Name")

### generate matrix of variants for plotting
summary_table <- names_var %>%
  group_by(dna_variants, Variants) %>%
  summarise(count = n()) %>%
  arrange(dna_variants, Variants)

wide_table <- summary_table %>%
  pivot_wider(names_from = Variants, values_from = count, values_fill = list(count = 0))

matrix_var <- as.matrix(wide_table[,-1])  # Exclude the first column (dna_variants) for the matrix
rownames(matrix_var) <- wide_table$dna_variants

rownames(matrix_var) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII",
                          "XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI")


### generate matrix of pop for plotting
summary_table <- names_var %>%
  group_by(dna_variants, pop ) %>%
  summarise(count = n()) %>%
  arrange(dna_variants, pop)

wide_table <- summary_table %>%
  pivot_wider(names_from = pop, values_from = count, values_fill = list(count = 0))


matrix_pop <- as.matrix(wide_table[,-1])  # Exclude the first column (dna_variants) for the matrix
rownames(matrix_pop) <- wide_table$dna_variants

rownames(matrix_pop) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII",
                          "XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI")


## generate basic haplonetwork
h<-haplotype(ali)

d <- dist.dna(h, "N")
nt <- haploNet(h)

plot(nt)
(sz <- summary(h))
(nt.labs <- attr(nt, "labels"))
sz <- sz[nt.labs]
sz_prop <- sz/sum(sz)

#plot(nt,xy=xy,size=sz_prop*10,labels=TRUE,pie = V,scale.ratio=1.5,bg=gbg,threshold = c(0,1))

#plot by pop
P <- matrix_pop[nt.labs, ]
V <- matrix_var[nt.labs, ]

setHaploNetOptions(link.width=1,link.width.alt=0.0,mutations.sequence.length=0.01,mutations.sequence.width=0.5)# this does not wrk properly

pdf("AvrPm17_hap_networks.pdf")
par(mfrow=c(1,2))
gbg <- c("#984EA3","#377EB8","#EA9999","#E41A1C","#E5B110")
plot(nt,xy=xy,size=sz_prop*10,labels=FALSE,pie = P,scale.ratio=1.5,bg=gbg,threshold = c(0,1))
#B,C,A,H
gbg <- c("#556B2F", "#8B4500","#76EEC6", "#EE9A49", "gray","gray","gray","gray","gray")
plot(nt,xy=xy,size=sz_prop*10,labels=FALSE,pie = V,scale.ratio=1.5,bg=gbg,threshold = c(0,1))
dev.off()


#### Fig 5c,d ####
# relatedness networks
library(isoRelate)
library(ggplot2)
library(maps)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(patchwork)
library(ggpubr)

load("../AvrPm17_isoRelate/BgtE+_avrpm17_geno.RData")

## make clusters with isoRelate
my_i_clusters <- getIBDiclusters(ped.genotypes = my_genotypes, 
                                 ibd.segments = my_ibd, 
                                 interval = c("1", 4365017, 4370466), # From first nucleotide of CDS of first copy, to last nucleotide of CDS of second copy
                                 prop=1, 
                                 hi.clust = FALSE)

g<-my_i_clusters$i.network

save(my_i_clusters,file="clusters_avrpm17.RData")
Samples <- my_genotypes$pedigree$fid
samples <- paste(Samples, Samples, sep = "/")

my_groups <- my_genotypes[[1]][,1:3]

# read variants (mature protein aa)
variants<-read.csv("../AvrPm17_haplotypes/hap_class_mp_trns_no_miss.csv")
colnames(variants) <- c("variants","Sample.Name")
variants$Sample.Name <- gsub("_1_avrpm17_cds", "",variants$Sample.Name)
variants$Sample.Name <- gsub("_2_avrpm17_cds", "",variants$Sample.Name)

# read dna variants
dna_variants<-read.csv("../AvrPm17_haplotypes/hap_class_cds_RC_no_miss.csv")
colnames(dna_variants) <- c("dna_variants","Sample.Name")
dna_variants$Sample.Name <- gsub("_1_avrpm17_cds", "",dna_variants$Sample.Name)
dna_variants$Sample.Name <- gsub("_2_avrpm17_cds", "",dna_variants$Sample.Name)

meta<-read.csv("../Datasets/S1_Data.csv")
meta1<-data.frame(meta$Sample_ID,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$year_of_collection)
colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Year")

groups<- merge(my_groups,meta1,by.x="fid",by.y="Sample.Name")
groups$merged <- paste(groups$fid, groups$iid, sep = "/")
colnames(groups) <- c("fid","iid","pid","Country","Population","Longitude","Latitude","Year","ID")

vertices<-V(g)$name

## identify samples not in the graph (not clustered) and add them
unique_elements <- setdiff(samples, vertices)
length(unique_elements)
g <- add_vertices(g, length(unique_elements), attr = list(name = unique_elements))
vertices<-V(g)$name
groups_s <- groups[match(vertices, groups$ID), ]

#### parse through variants and recode based on nomenclature Mueller et al.
var <- data.frame(Sample.Name=groups_s$fid)
var$var <-""
for (i in 1:nrow(var)){
  flag=0
  for (z in 1:nrow(variants)){
    if (var$Sample.Name[i] == variants$Sample.Name[z]){
      flag = flag+1
      if (flag > 1){
        var$var[i] = paste0(var$var[i],"/",variants$variants[z])
      }
      if (flag ==1){
        var$var[i] = variants$variants[z]
      }
    }
  }
}

var

for (i in 1:nrow(var)){
  if (var$var[i] == ""){var$var[i]  <- "NA"}
}
var$var <- gsub("0/0", "B",var$var)
var$var <- gsub("1/1", "C",var$var)
var$var <- gsub("2/2", "A",var$var)
var$var <- gsub("0/1", "B/C",var$var)
var$var <- gsub("0/2", "A/B",var$var)
var$var <- gsub("1/2", "A/C",var$var)
var$var <- gsub("1/4", "C/I",var$var)
var$var <- gsub("2/6", "A/K",var$var)
var$var <- gsub("2/7", "A/F",var$var)
var$var <- gsub("1/8", "C/L",var$var)
var$var <- gsub("1/5", "C/J",var$var)
var$var <- gsub("0", "B",var$var)
var$var <- gsub("1", "C",var$var)
var$var <- gsub("2", "A",var$var)
var$var <- gsub("3", "H",var$var)

# var
# 0 = varB
# 1 = varC
# 2 = varA
# 3 = varH  Y31H + A53V compared to varA
# 4 = varI  Y31H + A53V + R80S compared to varA
# 5 = varJ  R80S
# 6 = varK E55R + G61A
# 7 = varF
# 8 = varL E55R + G61A + R80S

#### parse through dna_variants
dna_var <- data.frame(Sample.Name=groups_s$fid)
dna_var$dna_var <-""
for (i in 1:nrow(dna_var)){
  flag=0
  for (z in 1:nrow(dna_variants)){
    if (dna_var$Sample.Name[i] == dna_variants$Sample.Name[z]){
      flag = flag+1
      if (flag > 1){
        dna_var$dna_var[i] = paste0(dna_var$dna_var[i],"/",dna_variants$dna_variants[z])
      }
      if (flag ==1){
        dna_var$dna_var[i] = dna_variants$dna_variants[z]
      }
    }
  }
}

for (i in 1:nrow(dna_var)){
  if (dna_var$dna_var[i] == ""){dna_var$dna_var[i]  <- "NA"}
}

unique(dna_var$dna_var)

# merge everything together and reorder
groups_sor<- merge(groups_s,var,by.x="fid",by.y="Sample.Name")
groups_sort<- merge(groups_sor,dna_var,by.x="fid",by.y="Sample.Name")

vertices<-V(g)$name
groups_sorted <- groups_sort[match(vertices, groups_sort$ID), ]

######## plot all clusters
## attach attributes to graph
vertex_attr(g, "Population", index = V(g)) <- groups_sorted$Population
vertex_attr(g, "Variant", index = V(g)) <- groups_sorted$var
vertex_attr(g, "dna_Variant", index = V(g)) <- groups_sorted$dna_var

#Define a color mapping
color_mapping_var <- c("A" = "#76EEC6", "B" = "#556B2F","C" = "#8B4500", "H" = "#EE9A49", "A/B" = "#EE3A8C", "A/C" = "#008B8B","B/C" = "#458B00","A/K" = "purple","A/F" = "violet","C/L" = "gold", "NA" = "gray", "C/I" = "#7FFF00","C/J" = "antiquewhite")
color_mapping_pop <- c("N_EUR" = "#377EB8", "S_EUR1" = "#EA9999","S_EUR2" = "#E41A1C", "TUR" = "#E5B110", "ME" = "#984EA3")
vertex_colors_pop <- color_mapping_pop[V(g)$Population]
vertex_colors_var <- color_mapping_var[V(g)$Variant]
load("../AvrPm17_isoRelate/layout.Rdata")
pdf("Avrpm17_ibd_all_clusters_main.pdf")
par(oma = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

par(mfrow = c(1, 2))
par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=3,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_pop
)

par(mar = c(0, 0, 0, 0) + 0.1) # Adjust as needed (bottom, left, top, right)

plot(g,
     layout=layout,
     #vertex.shapes=
     vertex.size=3,
     vertex.label.cex=0.5,
     vertex.label.dist=0.4,
     vertex.label.color="black",
     vertex.label=NA,
     vertex.frame.width=0.5,
     edge.width=0.3,
     vertex.color=vertex_colors_var
)

dev.off()
