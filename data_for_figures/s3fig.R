#### S3 Fig ####
library(tidyverse)
library(data.table)

log<-read_delim("../ADMIXTURE/CV_errors_10reps_rep_col",delim=":",col_names = FALSE)[,c(2:3)]
log$X2<-gsub("\\(K=", "", log$X2)
log$X2<-gsub(")", "", log$X2)
#interpret K values as numerical
log$X2<-as.numeric(log$X2)
colnames(log)<-c("k","cv")
#find number of replicates
rep<-length(which(log$k==1))
k<-max(as.numeric(log$k))

mean_df <- log %>%  
  group_by(k) %>%  
  summarize(average = mean(cv)) %>% 
  ungroup()
  
p <- ggplot(log, aes(x=as.factor(k), y=cv)) + 
  geom_boxplot()+
  geom_point()+
  ylab("Cross-validation error")+
  xlab("K")+
  theme_classic(base_size=20)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

