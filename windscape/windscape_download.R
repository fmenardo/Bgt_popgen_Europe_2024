#devtools::install_github("matthewkling/windscape")

.libPaths(c("/home/fmenar/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib/R/site-library","/usr/local/lib/R/library" ))

library(argparser)#, lib = "~/data/R_lib")

parser <- arg_parser("run windscape analysis")

parser <- add_argument(parser,"-y", nargs= 1,help="min value for interval years between 1979 and 2010 ")
parser <- add_argument(parser,"-Y", nargs= 1,help="max value for intervalyears between 1979 and 2010 ")
parser <- add_argument(parser,"-m", nargs= 1,help="min value for interval months to be included")
parser <- add_argument(parser,"-M", nargs= 1,help="max value for interval months to be included")
parser <- add_argument(parser,"-d", nargs= 1,help="min value for interval days to be included")
parser <- add_argument(parser,"-D", nargs= 1,help="max value for interval days to be included")

parser <- add_argument(parser,"-L", nargs= 1,help="left longitude in degrees (0-360)")
parser <- add_argument(parser,"-R", nargs= 1,help="right longitude in degrees (0-360)")
parser <- add_argument(parser,"-S", nargs= 1,help="min lat in degrees (0-90). will work out of the box only for N emisphere")             # will work out of the box only for N emisphere
parser <- add_argument(parser,"-N", nargs= 1,help="max lat in degrees (0-90). will work out of the box only for N emisphere")
parser <- add_argument(parser,"-o", nargs= 1,help="output stem")
#parser <- add_argument(parser,"-c", nargs= 1,help=" n cores", default=1)
#parser <- add_argument(parser,"-i", nargs= 1, default=0, help="if = 0 (default) performs complete analysis, if = 1 start from matrix calculation. i=2 analysis start after calculation of pairwise matrix")



args <- parse_args(parser)

#print (args)


library(windscape)
library(raster)
library(tidyverse)
#library(gdistance)
library(abind)
library(dplyr)
 



cfsr_dl_FM <- function(years = 1979, months = 1, days = 1,
                       xlim = c(260, 270), ylim = c(40, 50)){
  require(terra)
  require(purrr)
  q <- expand.grid(day = days, month = months, year = years,
                   #                 xmin = min(xlim), xmax = max(xlim), ymin = min(ylim), ymax = max(ylim))
                   xmin = xlim[1], xmax = xlim[2], ymin = min(ylim), ymax = max(ylim))
  
  
  q1 <- filter(q, day < 29)
  q2 <- filter(q, day == 29 & year%%4 == 0 & month == 2)
  q3 <- filter(q, day == 29 & month != 2)
  q4 <- filter(q, day == 30 & month != 2)
  q5 <- filter(q, day == 31 & month != 4 & month != 6 & month != 9 & month != 11 & month !=2)
  
  q <- rbind(q1,q2,q3,q4,q5)
#  print (q)
  w <- pmap(q, cfsr_dl_day_FM)
#  message("here")
  c(rast(map(w, "u")), rast(map(w, "v")))
}





cfsr_dl_day_FM <- function(year = 1979, month = 1, day = 1,
                           xmin = 260, xmax = 270, ymin = 40, ymax = 50){
  
  if(! year %in% 1979:2023) stop("'year' must be between 1979 and 2010")
  if(! month %in% 1:12) stop("'month' must be between 1 and 12")
  if(min(xmin, xmax) < 0 | max(xmin, xmax) > 360) stop("longitude and latitude must be between 0 and 360")
  require(ncdf4)
  require(lubridate)
  require(terra)
  
  # open connection
  # example url: "https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/1990/wnd10m.gdas.199001.grb2"
  if (year < 2011){
    url <- paste0("https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds093.1/", year, "/wnd10m.gdas.",
                  year, stringr::str_pad(month, 2, "left", 0), ".grb2")
  }
  if (year > 2010){
    url <- paste0("https://thredds.rda.ucar.edu/thredds/dodsC/files/g/ds094.1/", year, "/wnd10mx0.5.cdas1.",
                  year, stringr::str_pad(month, 2, "left", 0), ".grb2")
  }
  ds <- nc_open(url)
  
  # dimensions
  lon <- ncvar_get(ds, "lon")
  lon1 <- lon-180
  lat <- ncvar_get(ds, "lat")
  time <- ncvar_get(ds, "time")
  
  # convert time
  t_units <- ncatt_get(ds, "time", "units")
  startdate <- gsub("T00:00:00Z", "", unlist(strsplit(t_units$value, " "))[3])
  timestamp <- ymd(startdate) + dhours(time-1)
  
  # bounds
  x <- c(xmin, xmax)
  # x[x<0] <- x[x<0] + 360
  y <- c(ymin, ymax)
  date <- paste0(year, "-",
                 stringr::str_pad(month, 2, "left", 0), "-",
                 stringr::str_pad(day, 2, "left", 0))
  t <- as.POSIXct(paste0(date, c(" 00:00:00 UTC", " 23:00:00 UTC")), tz = "UTC")
  
  # indices
  btw <- function(data, z) range(which(data <= max(z) & data >= min(z)))
  btw1<- function(data,z) range(which(data >= max(z) & 360))   
  btw2<- function(data,z) range(1, which(data <= min(z))) 
  
  if (xmin > xmax){
    lon_i <- btw1(lon, x)
    lon_i_2 <- btw2(lon, x)
  }else{  lon_i <- btw(lon, x)}
  lat_i <- btw(lat, y)
  time_i <- btw(timestamp, t)
  lon_count <- lon_i[-1] - lon_i[1] + 1
  lat_count <- lat_i[-1] - lat_i[1] + 1
  time_count <- time_i[-1] - time_i[1] + 1
  start <- c(lon_i[1], lat_i[1], 1, time_i[1]) # x,y,z,t
  count <- c(lon_count, lat_count, 1, time_count)
  
  
  # get data
  message(paste0("... retrieving CFSR data for ", date, " ..."))
  v <- ncvar_get(ds, "v-component_of_wind_height_above_ground", start = start, count = count)
  
  u <- ncvar_get(ds, "u-component_of_wind_height_above_ground", start = start, count = count)
  
  # get second part of data if xmin > xmax
  if (xmin > xmax){
    start <- c(lon_i_2[1], lat_i[1], 1, time_i[1]) # x,y,z,t
    lon_count_2 <- lon_i_2[-1] - lon_i_2[1] + 1
    
    count_2 <- c(lon_count_2, lat_count, 1, time_count)
    v1 <- ncvar_get(ds, "v-component_of_wind_height_above_ground", start = start, count = count_2)
    u1 <- ncvar_get(ds, "u-component_of_wind_height_above_ground", start = start, count = count_2)
    
    # Combine arrays along the first dimension (lon)
    
    V <- abind(v, v1, along = 1)
    U <- abind(u, u1, along = 1)
  }else{
    V <- v
    U <- u
  }  
  
  
  # convert to raster object
  yres <- base::diff(sort(y)) / (lat_count - 1)
  
  if (xmin > xmax){
    
    xres <- (360-xmin+xmax) / (lon_count_2 + lon_count - 2)
    
  }else{xres <- base::diff(sort(x)) / (lon_count - 1)}
  
  
  # if longitide between 180 and 360 transform to negative (-180:0)
  if (x[1] > 180){x1diff = 360}else{x1diff = 0}
  if (x[2] > 180){x2diff = 360}else{x2diff = 0}
  
  
  extent <- c(x[1]-x1diff - 0.5 * xres, x[2]-x2diff + 0.5 * xres,
              y[1] - 0.5 * yres, y[2] + 0.5 * yres)
  
  vr <- rast(aperm(v, c(2, 1, 3)), extent = extent)
  ur <- rast(aperm(u, c(2, 1, 3)), extent = extent)
  names(vr) <- paste("v", timestamp[time_i[1]:time_i[2]])
  names(ur) <- paste("u", timestamp[time_i[1]:time_i[2]])
  return(list(u = ur, v = vr))
}


rm(wind)
wind <- cfsr_dl_FM(years = as.integer(args$y[1]):as.integer(args$Y[1]), months = as.integer(args$m[1]):as.integer(args$M[1]), days = as.integer(args$d[1]):as.integer(args$D[1]),
                   xlim = c(as.integer(args$L[1]), as.integer(args$R[1])), #longitudes to be in [0, 360] range for CFSR
                   ylim = c(as.integer(args$S[1]), as.integer(args$N[1])))

wind <- stack(wind)

#save(wind, file = paste0(args$o,".wind.Rdata"))
saveRDS(wind, file = paste0(args$o,".wind.rds"))

sink(file = paste0(args$o,".session_info.txt"))
sessionInfo()
sink() 
