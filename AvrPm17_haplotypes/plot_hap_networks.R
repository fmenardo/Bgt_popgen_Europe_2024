setwd("~/projects/project_avrpm17/haplotypes/")
library(adegenet)
library(pegas)
library(seqinr)
library(dplyr)
library(tidyr)


load("hap_net_xy.Rdata")  # pre cooked coordinates for hap net plotting, can be commented out, plot will be less pretty

ali<-fasta2DNAbin("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa")

fasta<-read.fasta("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa",as.string=TRUE,set.attributes=FALSE)

### fetch list of isolates in ali
system("grep \">\" avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa | sed \'s/>//g\'  > list_samples_for_h")
names<-read.table("list_samples_for_h",header=FALSE)


temp <- gsub("_1_avrpm17_cds","",names[,1])
names[,2] <- gsub("_2_avrpm17_cds","",temp)


## read metainfo

meta<-read.csv("~/projects/vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv")
meta1<-data.frame(meta$Sample.Name,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$Year.of.Collection)
colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Year")

## assign pop to h
colnames(names)<-c("Hap.Name","Sample.Name")

for (i in 1:nrow(names)){
  for (z in 1:nrow(meta1)){
    if (names$Sample.Name[i] == meta1$Sample.Name[z]){names$pop[i] = meta1$Population[z]}
  }
}


##assign varuiant nomenclature of Mueller et al, to haplotype


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


dna_variants<-read.csv("~/projects/project_avrpm17/haplotypes/hap_class_cds_RC_no_miss.csv")
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



#xy<- replot() # save coordinates and change manually, the order is the same of h

#xy$x[15] <- -4.5

#save(xy,file="hap_net_xy.Rdata")

