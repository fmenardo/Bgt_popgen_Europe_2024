setwd("~/projects/project_avrpm17/haplotypes/")
library(adegenet)
library(pegas)
library(seqinr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)



load("hap_net_xy.Rdata")  # pre cooked coordinates for hap net plotting, can be commented out, plot will be less pretty

ali<-fasta2DNAbin("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa")

#fasta<-read.fasta("avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa",as.string=TRUE,set.attributes=FALSE)

### fetch list of isolates in ali
system("grep \">\" avrpm17_ext_eur_only_cds_RC_unique_no_miss.fa | sed \'s/>//g\'  > list_samples_for_h")
names<-read.table("list_samples_for_h",header=FALSE)


temp <- gsub("_1_avrpm17_cds","",names[,1])
names[,2] <- gsub("_2_avrpm17_cds","",temp)


## read metainfo

meta<-read.csv("~/projects/vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK9_27062024.csv")
meta1<-data.frame(meta$Sample.Name,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$Year.of.Collection)
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



#xy<- replot() # save coordinates and change manually, the order is the same of h

#xy$x[15] <- -4.5

#save(xy,file="hap_net_xy.Rdata")


###############
## haplotypes distributions

names_var1<- merge(names_var,meta1,by="Sample.Name")

haplotypes <- subset(names_var1, select=-c(Hap.Name,dna_variants))
nrow(haplotypes)
unique(haplotypes$Variants)

I<- subset(haplotypes,haplotypes$Variants =="I")
length(unique(I$Sample.Name))

J<- subset(haplotypes,haplotypes$Variants =="J")
length(unique(J$Sample.Name))

K<- subset(haplotypes,haplotypes$Variants =="K")
length(unique(K$Sample.Name))

F<- subset(haplotypes,haplotypes$Variants =="F")
length(unique(F$Sample.Name))

L<- subset(haplotypes,haplotypes$Variants =="L")
length(unique(L$Sample.Name))

H<- subset(haplotypes,haplotypes$Variants =="H")
length(unique(H$Sample.Name))

C<- subset(haplotypes,haplotypes$Variants =="C")
length(unique(C$Sample.Name))

B<- subset(haplotypes,haplotypes$Variants =="B")
length(unique(B$Sample.Name))

A<- subset(haplotypes,haplotypes$Variants =="A")
length(unique(A$Sample.Name))


recent_h<- subset(haplotypes,haplotypes$Year>2014)
recent_h <- recent_h %>% filter(!is.na(Variants))

recent_h <- recent_h[!is.na(recent_h$Variants),]
nrow(recent_h)

# remove duplicates, these are isolates with two copies of same 
recent_h <-recent_h[!duplicated(recent_h), ]
nrow(recent_h)
write.csv(recent_h,file= "variants_h_avrpm17_coord.csv")

length(unique(recent_h$Sample.Name))

countA <-subset(recent_h,recent_h$Variants =="A")
length(unique(countA$Sample.Name))
59/362


countB <-subset(recent_h,recent_h$Variants =="B")
length(unique(countB$Sample.Name))
112/362

countC <-subset(recent_h,recent_h$Variants =="C")
length(unique(countC$Sample.Name))
218/362

#### by pop


recent_h$Variants <- factor(recent_h$Variants, levels = c("J", "I", "H","C","B","A" ))
recent_h$pop <- factor(recent_h$pop, levels = c("N_EUR","S_EUR2", "TUR","S_EUR1" ,"ME"))
nrow(recent_h)



p <- ggplot(recent_h, aes(x=pop, fill=Variants)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Number of isolates") +
  scale_fill_manual(values = c("B"= "#556B2F", "C" = "#8B4500", "A" = "#76EEC6", "H" = "#EE9A49", "J" = "gray", "I" = "gray" ),
                    labels=c("others","H","C","B","A")) +
  theme_classic() +
  theme(axis.text=element_text(size=9),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(colour = "black", fill = NA)) +
  theme(legend.text = element_text(size = 10),
        legend.title=element_text(size=10))

p
#### temporal variation
temporal_ds <- subset(haplotypes,haplotypes$Country == "Switzerland" |haplotypes$Country == "France" |haplotypes$Country == "United Kingdom" )
temporal_ds <- temporal_ds[!is.na(temporal_ds$Variants),]
temporal_ds <-temporal_ds[!duplicated(temporal_ds), ]

unique(temporal_ds$Sample.Name)


for (i in 1:nrow(temporal_ds)){
  if (temporal_ds$Year[i] < 2002){temporal_ds$period[i] = "old"}
  if (temporal_ds$Year[i] > 2021){temporal_ds$period[i] = "new"}
  if (temporal_ds$Year[i] > 2001 & temporal_ds$Year[i] < 2022){temporal_ds$period[i] = "middle"}
}

temporal_ds$period <- factor(temporal_ds$period, levels = c("old", "middle", "new"))
temporal_ds$Variants <- factor(temporal_ds$Variants, levels = c("J", "I", "H","C","B","A" ))
y <- ggplot(temporal_ds, aes(x=period, fill=Variants)) +
  geom_bar() +
  scale_x_discrete("", labels= c("1980-2001","2007-2020","2022-2023")) +
  scale_y_continuous("Number of isolates") +
  scale_fill_manual(values = c("B"= "#556B2F", "C" = "#8B4500", "A" = "#76EEC6", "H" = "#EE9A49", "J" = "gray", "I" = "gray" ),
                    labels=c("others","H","C","B","A")) +  theme_classic() +
  theme(axis.text=element_text(size=9),
        axis.title.y=element_text(size=12),
        panel.border = element_rect(colour = "black", fill = NA))+
  theme(legend.position = "none")
ggarrange(p,y)

ggsave("AvrPm17_Variants_barplot.pdf")




I_early<- subset(temporal_ds,temporal_ds$Variants =="I" & temporal_ds$period =="old")
length(unique(I_early$Sample.Name))
J_early<- subset(temporal_ds,temporal_ds$Variants =="J" & temporal_ds$period =="old")
length(unique(J_early$Sample.Name))
I_new<- subset(temporal_ds,temporal_ds$Variants =="I" & temporal_ds$period =="new")
length(unique(I_new$Sample.Name))
J_new<- subset(temporal_ds,temporal_ds$Variants =="J" & temporal_ds$period =="new")
length(unique(J_new$Sample.Name))



#### test variant A frequency


#### this is not the number of isolates, but the sum of haplotypes found in one isolate
tor_early=20
tot_new = 64

tot_sample_early<- subset(temporal_ds, temporal_ds$period =="old")
length(unique(tot__sample_early$Sample.Name))

tot_sample_new<- subset(temporal_ds, temporal_ds$period =="new")
length(unique(tot_sample_new$Sample.Name))



A_early<- subset(temporal_ds,temporal_ds$Variants =="A" & temporal_ds$period =="old")
length(unique(A_early$Sample.Name))

6/20

A_new<- subset(temporal_ds,temporal_ds$Variants =="A" & temporal_ds$period =="new")
length(unique(A_new$Sample.Name))



6/63


varA_freq<-matrix(c(6, 14, 6, 57),
                  nrow = 2,
                  dimnames = list(Var = c("varA", "not_varA"),
                                  period = c("old", "new")))





fisher.test(varA_freq, alternative = "greater")


#### test variant C frequency change

C_early<- subset(temporal_ds,temporal_ds$Variants =="C" & temporal_ds$period =="old")
length(unique(C_early$Sample.Name))

11/20

C_new<- subset(temporal_ds,temporal_ds$Variants =="C" & temporal_ds$period =="new")
length(unique(C_new$Sample.Name))

45/63


varC_freq<-matrix(c(11, 9, 45, 18),
                  nrow = 2,
                  dimnames = list(Var = c("varC", "not_varC"),
                                  period = c("old", "new")))

fisher.test(varC_freq, alternative = "less")

#### test variant B frequency change

B_early<- subset(temporal_ds,temporal_ds$Variants =="B" & temporal_ds$period =="old")
length(unique(B_early$Sample.Name))

2/20

B_new<- subset(temporal_ds,temporal_ds$Variants =="B" & temporal_ds$period =="new")
length(unique(B_new$Sample.Name))

9/63


varB_freq<-matrix(c(2, 18, 9, 54),
                  nrow = 2,
                  dimnames = list(Var = c("varB", "not_varB"),
                                  period = c("old", "new")))

fisher.test(varB_freq, alternative = "less")


#### test variant H frequency change

H_early<- subset(temporal_ds,temporal_ds$Variants =="H" & temporal_ds$period =="old")
length(unique(H_early$Sample.Name))

0/20

H_new<- subset(temporal_ds,temporal_ds$Variants =="H" & temporal_ds$period =="new")
length(unique(H_new$Sample.Name))

2/63


varH_freq<-matrix(c(0, 20, 2, 61),
                  nrow = 2,
                  dimnames = list(Var = c("varH", "not_varH"),
                                  period = c("old", "new")))

fisher.test(varH_freq, alternative = "less")



