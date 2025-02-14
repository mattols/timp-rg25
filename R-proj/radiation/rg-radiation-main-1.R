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
topo.bounds(dpath, vpath)
tf_day = topo.sw.day.wrapper(day_string = "2021-09-22", dpath, save_path = NULL, save_day = NULL, single_day = TRUE)
tf_series = topo.sw.multiday.wrapper(month = "09", year = "2021", dpath, start_day = NULL, end_day = NULL, 
                                     save_path = NULL, save_day = NULL, return_grid = TRUE)

# PLOT TEST - quick
plot(tf_day, col = heat.colors(100))
plot(tf_series, col = heat.colors(100), 
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade"))  


# # #
# Merge rasters
base_path = "~/olson/glacier-data/rg-wasatch"
sub_fls <- c("central_wasatch", "north_wasatch1", "north_wasatch2", "south_wasatch")
save_path <- file.path(base_path, "wasatch_merge")
if(!dir.exists(save_path)){dir.create(save_path)}
mosaic.month.save(month_list = 1:12, base_path, sub_fls, save_path)

# # # 
# generate seasonal data
seasonal_wrap_standard(fpath = save_path, sub_folder = NULL)


# cp -r ~/olson/glacier-data/rg-wasatch/wasatch_merge ~/data/timp/rad_results
#





#

# SEASONAL INTEGRATIONS
# integrate full melt season
timp_files <- list.files("~/olson/glacier-data/rg-wasatch/timp/tf_series2/", full.names = T)
timp_files <- list.files(save_path, full.names = T)[1:12]
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
seas_name = "DJF"
writeRaster(djf, paste0("~/olson/glacier-data/rg-wasatch/wasatch_merge/clearsky_tf_season_", seas_name, ".tif") )
# writeRaster(djf, paste0("/uufs/chpc.utah.edu/common/home/u1037042/data/timp/rad_results/tf_series/clearsky_tf_season_", seas_name, ".tif") )
list.files("/uufs/chpc.utah.edu/common/home/u1037042/data/timp/rad_results/tf_series/")



# MERGE FINAL RASTER
#
fls <- c("central_wasatch", "north_wasatch1", "north_wasatch2", "south_wasatch")

mfpths <- list.files(file.path("~/olson/glacier-data/rg-wasatch", fls,"tf_series/"),
           pattern = "01", full.names = TRUE)

mr <- do.call(merge, lapply(mfpths, function(r) stack(r)) )
#
month_rast = extend(mr, dem)
names(month_rast) <- c("norm_insol", "tot_insol", "topo_factor", "topo_norm", "shade_hrs")
plot(month_rast[[c(1,2)]])

# NA values
sum(is.na(values(dem)))/ ncell(dem)
# 92% of data is NA
sum(!is.na(values(dem))) / ncell(dem) # 7.5%

# test - only keep non-NA
insol_norm_month <- matrix(NA, nrow = dem_NA_len, ncol = 3)
insol_norm_month[, 1] <- as.array(dem)[dem_NA_idx]
d = array.to.raster2(insol_norm_month[, 1], dem)

# did not work for single day? lines in raster / incorrect values



# view_rast = make.raster(viewf, dem)
# view_rast2 = mask(view_rast, dem)
# writeRaster(view_rast2, "~/data/timp/dem_20m/Wstch_gtr_2400m_20m_svf.tif")





# ERROR in central point assumption
dem <- raster(dpath) # FULL DEM
dpoly <- as(extent(dem), 'SpatialPolygons')
crs(dpoly) <- crs(dem)
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
k_noon * ang_distance * 6371000

# change in length of shadow of a 100 ft pole
100 * ( (1/tan((90 - lat_range*180/pi + declination(160) + 25)*pi/180)) - (1/tan((90 - lat_lon[1] + declination(160) + 25)*pi/180)) )
# -1.6869893 -0.8435355  0.0000000  0.8437373  1.6877965 (ft error in length of shadow)

# How off are shadows across a Landsat scene?
# what is the equation that shows the change in solar position with respect to latitude.



## THOUGHTS
# outputs:
# - Days of sun per year
# - 




