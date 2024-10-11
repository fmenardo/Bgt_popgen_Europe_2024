library(tidyverse)
library(rms)

## data preparation ##

# list of samples in the N_EUR and S_EUR2 datasets (only those collected in 2022-2023)
samps <- readLines("tritici_2022+2023_fs_level4_N_EUR_+_S_EUR2")

# file with sample name and corresponding population name
pops <- read.delim("tritici_ext_eur_recent_fs_level4_populations", header = FALSE, sep = "\t")
rownames(pops) <- pops$V1

# read in geographic distance matrix
geo_dist <- as.matrix(read.csv("geo_dist_matrix_tritici_europe_recent.csv", row.names = 1, header = TRUE))
geo_subs <- geo_dist[samps, samps]  # subset matrix 
geo_df <- broom::tidy(as.dist(geo_subs),diag=FALSE, upper = FALSE)   # convert to dataframe with each pair of individuals as a row

# response variable - 0 if isolates belong to same population, 1 if different
geo_df$diff_pop <- ifelse(pops[as.character(geo_df$item1),2]==pops[as.character(geo_df$item2),2], 0, 1)
geo_df$diff_pop <- as.factor(geo_df$diff_pop)

# read in wind dist matrix
wind_dist <- as.matrix(read.csv("windscape_2012-2021.wind_distance_sym.csv", row.names = 1, header = TRUE))
wind_subs <- wind_dist[samps, samps]  ## subset matrix
wind_df <- broom::tidy(as.dist(wind_subs),diag=FALSE, upper=FALSE)  # convert to dataframe

# read in clim dist matrix
clim_dist <- as.matrix(read.csv("clim_dist_unweighted_7pc.csv", row.names = 1, header  = TRUE))
clim_subs <- clim_dist[samps, samps] # subset matrix
clim_df <- broom::tidy(as.dist(clim_subs),diag=FALSE, upper=FALSE) # convert to dataframe

# combine data frames
all_df <- geo_df %>% inner_join(wind_df, by=c("item1", "item2")) %>% inner_join(clim_df, by = c("item1", "item2"))
colnames(all_df) <- c("Isolate1", "Isolate2", "GeoDist", "same_pop", "diff_pop", "WinDist","ClimDist")

## perform logistic regression ##
dd <- datadist(all_df)
options(datadist="dd")

geo_s <- lrm(diff_pop ~ GeoDist_sc, data = all_df)
wind_s <- lrm(diff_pop ~ WinDist_sc, data = all_df)
clim_s <- lrm(diff_pop ~ ClimDist_sc, data = all_df)
full_s <- lrm(diff_pop ~ ClimDist_sc+GeoDist_sc+WinDist_sc, data = all_df)  # multiple regression 
