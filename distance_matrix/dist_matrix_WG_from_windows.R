setwd("/shares/menardo.bgt.uzh/project_bgt_popgen/analysis/dist_mat/")
library(tidyverse)

modifiedSum <- function(x, y) {
  replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
}

loci_list <- list.files(path = "dist_mat_windows", pattern = "dist_mat_num_loci.csv", full.names = TRUE)
loci_mats <- lapply(loci_list, read.csv, header = TRUE, row.names = 1) 
loci_mats_matr <- lapply(loci_mats, as.matrix)
gw_loci_matr_2 <- as.matrix(Reduce(modifiedSum,loci_mats_matr))

num_diff_list  <- list.files(path = "dist_mat_windows", pattern="dist_mat_num_diffs.csv",full.names = TRUE)
diff_mats <- lapply(num_diff_list, read.csv, header = TRUE, row.names = 1) 
diff_mats_matr <- lapply(diff_mats, as.matrix)
gw_diff_mat <- Reduce('+', diff_mats_matr)

gw_dist_mat <- gw_diff_mat / gw_loci_matr_2

#write.csv(gw_dist_mat, "gw_dist_mat_prop_2022+before2022+2023+ncsu.csv")

#### get clones ####

gw_dist_mat <- as.matrix(read.csv("gw_dist_mat_prop_2022+before2022+2023+ncsu.csv", header = TRUE, row.names = 1))
gw_dist_df <- as.data.frame(as.vector(gw_dist_mat))
ggplot()+aes(gw_dist_df$`as.vector(gw_dist_mat)`)+geom_histogram(bins = 500, fill = "slategray2", color="slategray2" ,alpha=0.6)+
  geom_vline(xintercept = 9e-05, linetype = "dashed")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Pairwise distance (no. of SNPs normalised by total number of loci)", y = "count")+
  geom_text(aes(x=9e-05, label="9e-05\n", y=8000), colour="black", angle=90, size=4)
  #coord_cartesian(xlim = c(0,0.00015), ylim = c(0,250))#+
  


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

### chose threshold as 9e-05 

clone_pairs <-as.data.frame(finv(which(dist_matr <=9e-05), dist_matr)) 
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

