library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(data.table)
setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/admixture/2022+before2022+2023+ncsu/tritici_ALL_25kb_0.1_LDp_admixture/")

#### cv error ####
# individual runs
cv <- read_delim("CV_errors_10reps_rep_col",delim=":",col_names = FALSE)
cv$X2<-gsub("\\(K=", "", cv$X2)
cv$X2<-gsub(")", "", cv$X2)
cv$X2<-as.numeric(cv$X2)
colnames(cv)<- c("repl","k","cv_error")
cv$repl<-gsub("r","",cv$repl)
cv$repl<- as.numeric(cv$repl)

DT <- data.table(cv)
min_runs <- DT[ , .SD[which.min(cv_error)], by = k]
min_runs$k <- as.numeric(min_runs$k)

ggplot(data=min_runs, aes(k,cv_error))+geom_line()+geom_point()+
  ylab("Minimum CV error across 10 replicates")+
  xlab("K")+
  scale_x_continuous(breaks = c(1:10))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


ggplot(data=cv, aes(as.factor(k),cv_error, colour = as.factor(repl)))+
  geom_line(aes(group=repl),linewidth=1)+
  ylab("Cross-validation error")+
  xlab("K")+
  labs(colour = "Replicate")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_brewer(palette = "Paired")

# average over runs
log<-read_delim("CV_errors_10reps_rep_col",delim=":",col_names = FALSE)[,c(2:3)]
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
  #geom_line(data = mean_df,mapping = aes(x = k, y = average, group=1),color="red", linetype = "dashed")+
  stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
  ylab("Cross-validation error")+
  xlab("K")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

#### bar plots ####

# read .fam file with sample names
samps<-read.table("tritici_ALL_25kb_0.1_LDp.fam")[,1]

# desired order of samples in plot
tritici_order <- read.delim("../tritici_sorted_25042024", header = FALSE)

admx_ord <- function(anc_prop){
  df <- read.delim(anc_prop, sep = " ", header = FALSE)
  rownames(df) <- samps
  df_ord <- df[match(tritici_order$V1,rownames(df)),]
  df_ord_t <- t(as.matrix(df_ord))
  return(df_ord_t)
}

list_files <- c("r10/tritici_ALL_25kb_0.1_LDp.k4.r10.Q", "r8/tritici_ALL_25kb_0.1_LDp.k5.r8.Q", "r3/tritici_ALL_25kb_0.1_LDp.k6.r3.Q", 
                "r8/tritici_ALL_25kb_0.1_LDp.k7.r8.Q", "r2/tritici_ALL_25kb_0.1_LDp.k8.r2.Q", "r10/tritici_ALL_25kb_0.1_LDp.k9.r10.Q")

tritici_admx <- lapply(list_files, admx_ord)

# dummy barplot to get label segment coordinates for subsequent plots
bp <- barplot(tritici_admx[[1]])
# regions : Northern/Central Europe, Southern Europe, Russia/Central Asia, N. Turkey, S. Turkey, Israel, Egypt, China, Japan, USA/Australia, Argentina
start_samples <- c(bp[1],bp[206],bp[277],bp[313],bp[370],bp[381],bp[409],bp[444],bp[506],bp[516],bp[561])
end_samples <- c(bp[205],bp[276],bp[312],bp[369],bp[380],bp[408],bp[443],bp[505],bp[515],bp[560],bp[568])
start_segments <- start_samples + 0.5
end_segments <- end_samples - 0.5
text_pos<- start_segments+((end_segments - start_segments)/2)

col_pals <- list(c("#999999","#377EB8","#E41A1C","#4DAF4A"),
                 c("#999999","#FFFF33","#E41A1C","#377EB8","#4DAF4A"),
                 c("#FF7F00","#377EB8","#E41A1C","#4DAF4A","#FFFF33","#999999"),
                 c("#FF7F00","#4DAF4A","#FFFF33","#F781BF","#E41A1C","#999999","#377EB8"),
                 c("#377EB8","#F781BF","#FF7F00","#FFFF33","#999999","#E41A1C","#4DAF4A","#984EA3"),
                 c("#E41A1C","#984EA3","#4DAF4A","#87CEEB","#FF7F00","#377EB8","#FFFF33","#F781BF","#999999"))


six_plots <- function(list_of_six, colpal){
  par(mfrow = c(6,1), mar = c(1,1,0.2,0), oma=c(1.5,0,0.3,0))
  for (i in 1:6){
    if (i<6){
      barplot(list_of_six[[i]],col=colpal[[i]], border=colpal[[i]], las = 2, cex.names = 0.3, ylab = "",
              ylim = c(-0.1,1), axes = FALSE,xaxt = "n")
      segments(x0=start_segments,
               x1=end_segments,
               y0=rep(-0.04,11),y1=rep(-0.04,11),lwd = 2, col = "black")
      title(ylab = paste0("K=",i+3), line = -4.5, cex.lab = 1.8)
      } else{
      barplot(list_of_six[[i]],col=colpal[[i]], border=colpal[[i]], las = 2, cex.names = 0.3, ylab = "",
              ylim = c(-0.1,1), axes = FALSE,xaxt = "n")
      segments(x0=start_segments,
               x1=end_segments,
               y0=rep(-0.04,11),y1=rep(-0.04,11),lwd = 2, col = "black")
      title(ylab = paste0("K=",i+3), line = -4.5, cex.lab = 1.8)
      text(x = text_pos, y = -0.17,
           labels = c(paste("North/Central", "Europe", sep = "\n"),paste("South", "Europe", sep = "\n"), 
                      "RUS+", paste("North", "Turkey", sep = "\n"), paste("S", "TUR",sep = "\n"), 
                      "ISR", "EGY", "China", "JPN", paste("USA/","AUS",sep = "\n"), "ARG"),
           xpd = NA, srt = 0, adj = 0.45, cex = 1.5 )
    }
  }
}


six_plots(tritici_admx,col_pals)


# only K = 9
par(mar = c(1,4,1.5,0), oma=c(1.5,1,1.5,0))
barplot(tritici_admx[[6]],col=col_pals[[6]], border=col_pals[[6]], las = 2, cex.names = 0.3, ylab = "Ancestry proportion",
        ylim = c(-0.1,1), xaxt = "n")
segments(x0=start_segments,
         x1=end_segments,
         y0=rep(-0.04,11),y1=rep(-0.04,11),lwd = 2, col = "black")
text(x = text_pos, y = -0.08,
     labels = c("N/C EUR","S EUR", "RUS+", "N TUR", paste("S", "TUR",sep = "\n"), 
                "ISR", "EGY", "CHN", "JPN", paste("USA/","AUS",sep = "\n"), "ARG"),
     xpd = NA, srt = 0, adj = 0.45, cex = 1.5 )

