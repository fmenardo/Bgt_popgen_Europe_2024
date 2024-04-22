#devtools::install_github("matthewkling/windscape")

.libPaths(c("/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib/R/site-library","/usr/local/lib/R/library" ))

library(argparser)#, lib = "~/data/R_lib")

parser <- arg_parser("run windscape analysis")


parser <- add_argument(parser,"-L", nargs= 1,help="left longitude in degrees (0-360)")
parser <- add_argument(parser,"-R", nargs= 1,help="right longitude in degrees (0-360)")
parser <- add_argument(parser,"-S", nargs= 1,help="min lat in degrees (0-90). will work out of the box only for N emisphere")             # will work out of the box only for N emisphere
parser <- add_argument(parser,"-N", nargs= 1,help="max lat in degrees (0-90). will work out of the box only for N emisphere")
parser <- add_argument(parser,"-o", nargs= 1,help="output stem ")
parser <- add_argument(parser,"-c", nargs= 1,help=" n cores", default=1)
parser <- add_argument(parser,"-s", nargs= 1,help=" csv file with list of samples and coordinates column names must be Sample.Name,\"Longitude\",\"Latitude\"")
parser <- add_argument(parser,"-P", nargs= 1,help=" comma delimited list of coordinates in Lon Lat Lon Lat.. format for each pair (point)  a plot of wind connectivity will be produced")
parser <- add_argument(parser,"-l", nargs= 1,help="file with list of .rds files with the wind data (output of download script). One per line ")
parser <- add_argument(parser,"-i", nargs= 1, default=0, help="if = 0 (default) performs complete analysis, if = 1 only produces plots and do not calculate distance matrix")



args <- parse_args(parser)

#print (args)





library(windscape)
library(raster)
library(tidyverse)
library(gdistance)
#library(abind)
library(dplyr)


# convert longitude from 0-360 to -180 to 180
if (as.integer(args$L[1]) > 180){min_x = as.integer(args$L[1])-360}else{min_x = as.integer(args$L[1])}
if (as.integer(args$R[1]) > 180){max_x = as.integer(args$R[1])-360}else{max_x = as.integer(args$R[1])}

min_y = as.integer(args$S[1])
max_y = as.integer(args$N[1])


# process list of focal sites for plotting

if (!is.null(args$P) & (!is.na(args$P))) {

sites<-as.data.frame(matrix(unlist(strsplit(args$P,",")),ncol=2,byrow=T))

colnames(sites)<-c("Longitude","Latitude")

sites_clean <- filter(sites, Latitude >= min_y & Latitude <= max_y, Longitude >= min_x & Longitude <= max_x)
}

## read wind data
list_files<-read.table(args$l, header=F)

list_obj <- list()
for (i in 1:nrow(list_files)){
  list_obj[[i]] <- readRDS(list_files[i,1])
}

#load(paste0(args$o,".wind.Rdata"))

wind <- stack(list_obj)

# summarize wind time series (n = 1152 hourly layers) into 
# a wind conductance or "windrose" raster (n = 8 directional layers)
conductance <- windrose_rasters(wind, order = "uuvv", p = 1)


# generate downwind and upwind dispersal surfaces for focal sites 


if (!is.null(args$P)&(!is.na(args$P))){
  library(maps)
  
  # Load map data
  world <- map_data("world")
  
  for (i in 1:nrow(sites_clean)) {
    site <- matrix(as.numeric(sites_clean[i,]), ncol = 2)
    downwind <- build_wind_graph(conductance, "downwind") %>%
      accCost(site)
    upwind <- build_wind_graph(conductance, "upwind") %>%
      accCost(site)

    # restructure data, and plot
    d1 <- stack(downwind, upwind) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      gather(direction, wind_hours, layer.1, layer.2) %>%
      mutate(direction = recode(direction, 
                            layer.1 = "downwind (outbound)", 
                            layer.2 = "upwind (inbound)"))

plot <- ggplot(d1) +
      facet_wrap(~direction) +
      geom_raster(aes(x, y, fill = wind_hours), alpha = 0.9) +
      geom_contour(aes(x, y, z = wind_hours), bins = 20, color = "white", linewidth = .25) +
      geom_point(data = as.data.frame(site), aes(V1, V2)) +
      scale_x_continuous() +
      scale_y_continuous() +
      coord_fixed(ratio = 1.2) +
      scale_fill_gradientn(colors=c("yellow", "red", "blue", "black")) +
      theme_void() +
      theme(legend.position="right", strip.text=element_text(size=12),legend.direction='vertical',plot.margin=unit(c(0.5,0.5,0.5,0.5), 'cm'))+
      labs(fill= "wind hours")

# Add land outlines
   plot + 

     geom_path(data = world, aes(x = long, y = lat, group = group), color = "black")+
      xlim(c(min_x, max_x)) +  # Set specific x-axis limits
      ylim(c(min_y, max_y)) +
      coord_fixed(ratio = 1.2)

  ggsave(paste0(args$o,".site",i,".wind_map.png"))
  ggsave(paste0(args$o,".site",i,".wind_map.pdf"))
 
  }
}  
 



if (args$i==1){
 q() 
}  
  
  library(maps)

# read list of locations to consider

coord<- read.csv(args$s, header= TRUE)
#coord<- read.csv("samples.csv", header= TRUE)
coord_reord <- coord[, c("Sample.Name", "Longitude", "Latitude")]

coord_nona<- na.omit(coord_reord)

samples_name <-as.vector(coord_nona[,1])
coord1 <- coord_nona[,-1]

class(samples_name)
print(samples_name)

#rownames(coord1) <- coord$Sample.Name

#print(coord1)

wind_distance_matrix <- matrix(nrow = nrow(coord1), ncol=nrow(coord1))
wind_distance_matrix_average_sym <- matrix(nrow = nrow(coord1), ncol=nrow(coord1))

## calculate distance matrix between all

 for (o in 1:nrow(coord1)){
   print(o)
   for (a in o:nrow(coord1)){
     if (o == a){
       wind_distance_matrix[o,a]=0
       wind_distance_matrix_average_sym[o,a]=0
       }else{
       site1 <- matrix(as.numeric(coord1[o,]), ncol = 2)
       site2 <- matrix(as.numeric(coord1[a,]), ncol = 2)
       #print (site1)
       #print(site2)
       if (site1[1] == site2[1] &site1[2] == site2[2]){
         distance_down=0
         distance_up=0
       }else{
           distance_down <-(build_wind_graph(conductance, "downwind") %>% costDistance(site1,site2))
           distance_up <-(build_wind_graph(conductance, "downwind") %>% costDistance(site2,site1))
       }
       wind_distance_matrix[o,a] <- distance_down
       wind_distance_matrix[a,o] <- distance_up
       wind_distance_matrix_average_sym[o,a]= (distance_down+distance_up)/2
       wind_distance_matrix_average_sym[a,o]= (distance_down+distance_up)/2
     }
   }
 }
 
rownames(wind_distance_matrix) <- samples_name
colnames(wind_distance_matrix) <- samples_name

rownames(wind_distance_matrix_average_sym) <- samples_name
colnames(wind_distance_matrix_average_sym) <- samples_name


write.csv(wind_distance_matrix, file=paste0(args$o,".wind_distance_asym.csv"))
write.csv(wind_distance_matrix_average_sym, file=paste0(args$o,".wind_distance_sym.csv"))

sink(file = paste0(args$o,".analysis.session_info.txt"))
sessionInfo()
sink() 

