#
# Clear-sky Radiation Model - Rock Glaciers Wasatch Front, UT
# Matt Olson
# 01-10-2025
#
# R 4.03 Geospatial Packages (CHPC)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# call functions
rm(list=ls())
source('radiation/rg-radiation-func-1.R', echo=FALSE)


# RUN FUNCTIONS FOR TEST
topo.bounds(dpath)
tf_day = topo.sw.day.wrapper(day_string = "2021-09-22", dpath, save_path = NULL, save_day = NULL, single_day = TRUE)
tf_series = topo.sw.multiday.wrapper(month = "09", year = "2021", dpath, start_day = NULL, end_day = NULL, 
                                     save_path = NULL, save_day = NULL, return_grid = TRUE)

# PLOT TEST - quick
plot(tf_day, col = heat.colors(100))
plot(tf_series, col = heat.colors(100), 
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade"))  


# SEASONAL INTEGRATIONS
# integrate full melt season
timp_files <- list.files("~/olson/glacier-data/rg-wasatch/timp/tf_series2/", full.names = T)

# DJF
djf = seasonal_vals3_5(seasonal_path = timp_files, s1=12, s2=1, s3=2)
plot(djf, col = heat.colors(100), 
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade")) 
# MAM
mam = seasonal_vals3_5(seasonal_path = timp_files, s1=3, s2=4, s3=5)
plot(mam, col = heat.colors(100), 
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade")) 
# JJA
jja = seasonal_vals3_5(seasonal_path = timp_files, s1=6, s2=7, s3=8)
plot(jja, col = heat.colors(100), 
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade")) 
# SON
son = seasonal_vals3_5(seasonal_path = timp_files, s1=9, s2=10, s3=11)
plot(son, col = heat.colors(100), 
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade"))
# MJJAS
mjjas = seasonal_vals3_5(seasonal_path = timp_files, s1=5, s2=6, s3=7, s4=8, s5=9)
oma=c(0,0,0,2)
plot(mjjas, col = heat.colors(100), axes=F, box=F, mar = c(1,1,1,1), bg = "gray20",
     main = c("MJJAS Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade"))
# SAVE
seas_name = "SON"
writeRaster(son, paste0("/uufs/chpc.utah.edu/common/home/u1037042/data/timp/rad_results/tf_series/clearsky_tf_season_", seas_name, ".tif") )
list.files("/uufs/chpc.utah.edu/common/home/u1037042/data/timp/rad_results/tf_series/")

#






# ERROR in central point assumption
dem <- raster(dpath) # FULL DEM
# get coordinates in radians
if(gsub(".*proj=([a-z]+)\\s.*", '\\1', crs(dpoly) ) != "longlat"){pdpoly <- spTransform(dpoly, "EPSG:4326") }
lat_lon <- rowMeans(pdpoly@bbox)[2:1]
lat_lon2 = lat_lon * pi/180
lat_range <- seq(pdpoly@bbox[2, 1], pdpoly@bbox[2, 2], length.out = 5) * pi/180
lon_range <- seq(pdpoly@bbox[1, 1], pdpoly@bbox[1, 2], length.out = 5) * pi/180
#
# calculate angular distance
ang_distance = acos( sin(lat_lon2[1])*sin(lat_range) + cos(lat_lon2[1])*cos(lat_range)*cos(lon_range - lat_lon2[2]) )
#
# solar position error
k_high = 0.5 
k_morning = 0.01 # degrees per minute
k_noon = 0.005

# SZ error in meters
k_morning * ang_distance * 6371000

## THOUGHTS
# outputs:
# - Days of sun per year
# - 




