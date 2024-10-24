## dist along north-south / east-west axes

library(fields)
all_meta <- read.csv("S1_Data.csv")  ## read in all metadata
s_eur <- readLines("tritici_2022+2023_fs_level4_S_EUR2.args")   ## list of samples in Europe+_2022_2023_S_EUR2
n_eur <- readLines("tritici_2022+2023_fs_level4_N_EUR.args")	## list of samples in Europe+_2022_2023_N_EUR

# subset metadata	(get long-lat)
meta_neur <- all_meta[n_eur,c(8,7)]
meta_seur <- all_meta[s_eur,c(8,7)]

## calculate xy distance components ###
## function chatgpt 
compute_components <- function(lat1, lon1, lat2, lon2) {
  R <- 6371  # Earth's radius in km
  
  d_lat <- (lat2 - lat1) * pi / 180  # Difference in latitude in radians
  d_lon <- (lon2 - lon1) * pi / 180  # Difference in longitude in radians
  
  mean_lat <- (lat1 + lat2) / 2 * pi / 180  # Mean latitude in radians
  
  x_component <- abs(R * d_lon * cos(mean_lat))  # East-West distance component  (have to correct for less distance at poles and more at equator)
  y_component <- abs(R * d_lat)  # North-South distance component   (simply the difference in latitude)
  
  return(c(x_component, y_component))
}

## N_EUR

# Create matrices to store x and y components
n <- nrow(meta_neur)
x_matrix <- matrix(0, n, n)
y_matrix <- matrix(0, n, n)

# Fill the matrices with computed components
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      components <- compute_components(meta_neur$Latitude[i], meta_neur$Longitude[i], meta_neur$Latitude[j], meta_neur$Longitude[j])
      x_matrix[i, j] <- components[1]
      y_matrix[i, j] <- components[2]
    }
  }
}

rownames(x_matrix) <- rownames(meta_neur)
colnames(x_matrix) <- rownames(meta_neur)
rownames(y_matrix) <- rownames(meta_neur)
colnames(y_matrix) <- rownames(meta_neur)

write.csv(x_matrix, "geo_dist_X_n_eur_2022-2023.csv")
write.csv(y_matrix, "geo_dist_Y_n_eur_2022-2023.csv")


## S_EUR2

ns <- nrow(meta_seur)
x_matrix <- matrix(0, ns, ns)
y_matrix <- matrix(0, ns, ns)

# Fill the matrices with computed components
for (i in 1:ns) {
  for (j in 1:ns) {
    if (i != j) {
      components <- compute_components(meta_seur$Latitude[i], meta_seur$Longitude[i], meta_seur$Latitude[j], meta_seur$Longitude[j])
      x_matrix[i, j] <- components[1]
      y_matrix[i, j] <- components[2]
    }
  }
}

rownames(x_matrix) <- rownames(meta_seur)
colnames(x_matrix) <- rownames(meta_seur)
rownames(y_matrix) <- rownames(meta_seur)
colnames(y_matrix) <- rownames(meta_seur)

write.csv(x_matrix, "geo_dist_X_s_eur2_2022-2023.csv")
write.csv(y_matrix, "geo_dist_Y_s_eur2_2022-2023.csv")

