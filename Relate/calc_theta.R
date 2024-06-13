
setwd("~/projects/project_tritici_fabrizio/analysis/Relate/")
table1 <- read.table("2022-2023+out_LR026984.1_chr1.dist", sep = " " )
table2 <- read.table("2022-2023+out_LR026985.1_chr2.dist", sep = " " )
table3 <- read.table("2022-2023+out_LR026986.1_chr3.dist", sep = " " )
table4 <- read.table("2022-2023+out_LR026987.1_chr4.dist", sep = " " )
table5 <- read.table("2022-2023+out_LR026988.1_chr5.dist", sep = " " )
table6 <- read.table("2022-2023+out_LR026989.1_chr6.dist", sep = " " )
table7 <- read.table("2022-2023+out_LR026990.1_chr7.dist", sep = " " )
table8 <- read.table("2022-2023+out_LR026991.1_chr8.dist", sep = " " )
table9 <- read.table("2022-2023+out_LR026992.1_chr9.dist", sep = " " )
table10 <- read.table("2022-2023+out_LR026993.1_chr10.dist", sep = " " )
table11 <- read.table("2022-2023+out_LR026994.1_chr11.dist", sep = " " )




genome_size =sum(table1[,2])+sum(table2[,2])+sum(table3[,2])+sum(table4[,2])+sum(table5[,2])+sum(table6[,2])+sum(table7[,2])+sum(table8[,2])+sum(table9[,2])+sum(table10[,2])+sum(table11[,2])


genome_size

#snp count
SNP<-(system("wc -l 2022-2023+out_LR0269*.1_chr*.dist| tr -s ' '|cut -f 2 -d \" \"",intern = TRUE))
#SNP1<-as.integer(system("wc -l 2022-2023+out_LR026985.1_chr1.dist| cut -f 1 -d \" \"",intern = TRUE))
SNP<-as.integer(SNP)

a1 <- function(n){i <- 1:(n-1) ; return(sum(1/i))}
SNP[12]

theta <-SNP[12]/(a1(255)*genome_size)

theta

u=0.0000005
Ne=theta/(2*u)
Ne

