library(ggplot2)
library(dplyr)
library(stringr)

setwd("~/projects/nikos/LD_decay/2_output/100kb/")

# ALL

# import the data
ld <- read.table("LD_plink2_ALL_100kb.vcor", sep="", header=F)

# check column names
colnames(ld) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort <- ld %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort$dists <- cut(ld_summary_sort$dist,
                             breaks=seq(from=min(ld_summary_sort$dist)-1,
                                        to=max(ld_summary_sort$dist)+1,
                                        by=1000)) # by=10 makes it more detailed

dfr1 <- ld_summary_sort %>% group_by(dists) %>% summarise(mean=mean(R2),median=median(R2))

dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


# N_Europe

# import the data
ld_N_Europe <- read.table("LD_plink2_N_EUR_100kb.vcor", sep="", header=F)

# check column names
colnames(ld_N_Europe) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort_N_Europe <- ld_N_Europe %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort_N_Europe$dists <- cut(ld_summary_sort_N_Europe$dist,
                             breaks=seq(from=min(ld_summary_sort_N_Europe$dist)-1,
                                        to=max(ld_summary_sort_N_Europe$dist)+1,
                                        by=1000)) # by=10 makes it more detailed

dfr1_N_Europe <- ld_summary_sort_N_Europe %>% group_by(dists) %>% summarise(mean=mean(R2),median=median(R2))

dfr1_N_Europe <- dfr1_N_Europe %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


# S_EUR+

# import the data
ld_S_EUR <- read.table("LD_plink2_S_EUR+_100kb.vcor", sep="", header=F)

# check column names
colnames(ld_S_EUR) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort_S_EUR <- ld_S_EUR %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort_S_EUR$dists <- cut(ld_summary_sort_S_EUR$dist,
                             breaks=seq(from=min(ld_summary_sort_S_EUR$dist)-1,
                                        to=max(ld_summary_sort_S_EUR$dist)+1,
                                        by=1000)) # by=10 makes it more detailed

dfr1_S_EUR <- ld_summary_sort_S_EUR %>% group_by(dists) %>% summarise(mean=mean(R2),median=median(R2))

dfr1_S_EUR <- dfr1_S_EUR %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

# S_EUR1

# import the data
ld_S_EUR1 <- read.table("LD_plink2_S_EUR1_100kb.vcor", sep="", header=F)

# check column names
colnames(ld_S_EUR1) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort_S_EUR1 <- ld_S_EUR1 %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort_S_EUR1$dists <- cut(ld_summary_sort_S_EUR1$dist,
                                   breaks=seq(from=min(ld_summary_sort_S_EUR1$dist)-1,
                                              to=max(ld_summary_sort_S_EUR1$dist)+1,
                                              by=1000)) # by=10 makes it more detailed

dfr1_S_EUR1 <- ld_summary_sort_S_EUR1 %>% group_by(dists) %>% summarise(mean=mean(R2),median=median(R2))

dfr1_S_EUR1 <- dfr1_S_EUR1 %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                    end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                    mid=start+((end-start)/2))


# TUR

# import the data
ld_TUR <- read.table("LD_plink2_TUR_100kb.vcor", sep="", header=F)

# check column names
colnames(ld_TUR) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort_TUR <- ld_TUR %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort_TUR$dists <- cut(ld_summary_sort_TUR$dist,
                                   breaks=seq(from=min(ld_summary_sort_TUR$dist)-1,
                                              to=max(ld_summary_sort_TUR$dist)+1,
                                              by=1000)) # by=10 makes it more detailed

dfr1_TUR <- ld_summary_sort_TUR %>% group_by(dists) %>% summarise(mean=mean(R2),median=median(R2))

dfr1_TUR <- dfr1_TUR %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                    end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                    mid=start+((end-start)/2))

# ME

# import the data
ld_ME <- read.table("LD_plink2_ME_100kb.vcor", sep="", header=F)

# check column names
colnames(ld_ME) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort_ME <- ld_ME %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort_ME$dists <- cut(ld_summary_sort_ME$dist,
                             breaks=seq(from=min(ld_summary_sort_ME$dist)-1,
                                        to=max(ld_summary_sort_ME$dist)+1,
                                        by=1000)) # by=10 makes it more detailed

dfr1_ME <- ld_summary_sort_ME %>% group_by(dists) %>% summarise(mean=mean(R2),median=median(R2))

dfr1_ME <- dfr1_ME %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

LD_decay <- ggplot() +
  geom_line(data=dfr1_ME, aes(x=start, y=mean, colour="#FF7F00"), linewidth=0.6, alpha=0.8)+
  geom_line(data=dfr1_N_Europe, aes(x=start, y=mean, colour="#377EB8"),linewidth=0.6, alpha=0.8)+
  geom_line(data=dfr1_S_EUR, aes(x=start, y=mean, colour="#E41A1C"), linewidth=0.6, alpha=0.8)+
  geom_line(data=dfr1_S_EUR1, aes(x=start, y=mean, colour="#EA9999"), linewidth=0.6, alpha=0.8)+
  geom_line(data=dfr1_TUR, aes(x=start, y=mean, colour="#E5B110"), linewidth=0.6, alpha=0.8)+
  geom_line(data=dfr1, aes(x=start, y=mean, colour="black"), linewidth=0.6, alpha=0.8)+
  labs(x="Distance (Kb)", y=expression(LD~(r^{2})))+
  scale_colour_manual(name = 'Population', 
                      values = c('#FF7F00'='#FF7F00', '#377EB8'='#377EB8', '#E41A1C'='#E41A1C',
                                 '#EA9999'='#EA9999', '#E5B110'='#E5B110', 'black'='black'),
                      labels = c('N_EUR', 'S_EUR2', 'TUR', 'S_EUR1', 'ME', 'all')) +
  scale_x_continuous(breaks=c(0,2*10^4,4*10^4,6*10^4,8*10^4,10*10^4), labels=c("0","20","40","60","80","100"))+
  theme_bw()


pdf(paste("plots/LD_decay_ALL_100kb_by1000_fs_level4_colors_legend.pdf"), width = 8, height = 6)
print(LD_decay)
dev.off()

png(paste("plots/LD_decay_ALL_100kb_by1000_fs_level4_colors_legend.png"), width = 8, height = 6, units = "in", res = 150)
print(LD_decay)
dev.off()