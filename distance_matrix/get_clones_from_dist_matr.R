library(vcfR)
library(adegenet)
library(ape)
library(ggplot2)
library(plotly)
library(heatmaply)
library(Polychrome)
library(parallel)
library(pals)

setwd("~/projects/project_bgt_popgen/analysis/dist_mat")

##### get genome wide distance matrix from per-chromosome matrix 
chr_dist_mat_list <- list.files(pattern="ncsu_dist_mat.csv$",full.names = TRUE)

dist_mats <- lapply(chr_dist_mat_list, read.csv, header = TRUE, row.names = 1) 
dist_mats_matr <- lapply(dist_mats, as.matrix)

gw_dist_mat <- Reduce('+', dist_mats_matr)

head(gw_dist_mat)
dim(gw_dist_mat)

write.csv(gw_dist_mat, "gw_dist_mat_2022+before2022+2023+ncsu.csv")

##### plot histogram of distances 
# clones upto 119 SNPs. next values in matrix are 7439,14079,16539...
gw_dist_df <- as.data.frame(as.vector(gw_dist_mat))
ggplot()+aes(gw_dist_df$`as.vector(gw_dist_mat)`)+geom_histogram(bins = 500, fill = "slategray2", color="slategray2" ,alpha=0.4)+
  geom_vline(xintercept = 120, linetype = "dashed")+scale_x_continuous(breaks = seq(0,600000,50000))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Pairwise distance (no. of SNPs)")


#matr <- as.matrix(read.csv("dist_mat_2022+before2022_tritici.csv", header = TRUE, row.names = 1))
dist_matr <- as.dist(gw_dist_mat)

finv <- function (k, dist_obj) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (k >= 1) & (k <= n * (n - 1) / 2)
  k_valid <- k[valid]
  j <- rep.int(NA_real_, length(k))
  j[valid] <- floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k_valid
                                                               - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  #cbind(i, j)
  cbind(labels(dist_obj)[i], labels(dist_obj)[j])
}



# Function to merge overlapping lists (chatgpt)
merge_overlapping_lists <- function(lists) {
  result <- list()  # Initialize an empty list for the merged lists
  
  while (length(lists) > 0) {
    current_list <- lists[[1]]  # Take the first list
    lists <- lists[-1]  # Remove the first list from the original list
    
    # Check for overlap with other lists
    overlap_indices <- sapply(lists, function(lst) length(intersect(current_list, lst)) > 0)
    
    if (any(overlap_indices)) {
      # If there is an overlap, merge the current list with the overlapping lists
      overlapping_lists <- lists[overlap_indices]
      current_list <- unique(unlist(c(current_list, overlapping_lists)))
      lists <- lists[!overlap_indices]
    }
    
    # Add the current merged list to the result
    result <- append(result, list(current_list))
  }
  
  return(result)
}

#hist(dist_matr, xlim(0,5000))

### chose threshold as 120 SNPs (stays the same till 7438)
clone_pairs <-as.data.frame(finv(which(dist_matr <=120), dist_matr)) 


clone_pairs
list_clones<-unique(unlist(clone_pairs,use.names = FALSE))
list_clones
my_list <- list()

for(i in 1:length(list_clones)) {
  
  target_value <- list_clones[i]
  filtered_df <- clone_pairs %>% filter_all(any_vars(. == target_value))
  filtered_df
  
  ucl <- unique(unlist(filtered_df,use.names = FALSE))# Function to merge overlapping lists
  
  list(ucl)
  
  my_list <- c(my_list, list(ucl))
  
  
}

list_clones_merged<-merge_overlapping_lists(my_list)

sink("2022+before2022+2023+ncsu_all_list_clones_merged.txt")
print(list_clones_merged)
sink()


heatmaply(as.matrix(dist_matr),
row_dend_left = T,
plot_method = "plotly")#,
#          row_side_colors = clade_info$clade)


