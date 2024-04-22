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
library(ggpubr)

# convert longitude from 0-360 to -180 to 180
if (as.integer(args$L[1]) > 180){min_x = as.integer(args$L[1])-360}else{min_x = as.integer(args$L[1])}
if (as.integer(args$R[1]) > 180){max_x = as.integer(args$R[1])-360}else{max_x = as.integer(args$R[1])}

min_y = as.integer(args$S[1])
max_y = as.integer(args$N[1])


# process list of focal sites for plotting

if (!is.null(args$P) & (!is.na(args$P))) {

sites<-as.data.frame(matrix(unlist(strsplit(args$P,",")),ncol=2,byrow=T))

colnames(sites)<-c("Longitude","Latitude")

sites_clean <- filter(sites, Latitude >= min_y & Latitude <= max_y & Longitude >= min_x & Longitude <= max_x)
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
  plot_list <- list()
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
      theme(legend.position="right", strip.text=element_text(size=12),legend.direction='vertical',plot.margin=unit(c(0.25,0.25,0.25,0.25), 'cm'),
        legend.text=element_text(size=7), legend.key.size = unit(0.3, "cm"))+
      labs(fill= "wind hours")

# Add land outlines
plot_list[[i]] <-  plot + 

     geom_path(data = world, aes(x = long, y = lat, group = group), color = "black")+
      xlim(c(min_x, max_x)) +  # Set specific x-axis limits
      ylim(c(min_y, max_y)) +
      coord_fixed(ratio = 1.2)

  }
}  
 

ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],ncol=1,nrow=4)
ggsave("2012-2021_allsites_wind_map.png")
ggsave("2012-2021_allsites_wind_map.pdf")



