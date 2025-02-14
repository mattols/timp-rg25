#
# Clear-sky Radiation Model - Rock Glaciers Wasatch Front, UT
# Matt Olson
# 01-10-2025
#
# R 4.03 Geospatial Packages (CHPC)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

require(insol);require(raster);require(rgdal);require(rgeos)
require(sp);require(dplyr)
# require(ggplot2)
# require(ncdf4)

# wasatch grid
dpath = "~/data/timp/dem_20m/Wstch_gtr_2400m_20m_grid.tif"
vpath = "~/data/timp/dem_20m/Wstch_gtr_2400m_20m_svf.tif"
# may need to be split?

# WORKSPACE STORAGE
timp_storage = "~/olson/glacier-data/rg-wasatch/timp" # UVU group
rg_storage <- "~/olson/glacier-data/rg-wasatch/central_wasatch"
ftype = ".tif" # ".grd"
sub_folder = "tf_series" # For multiday
# "~/molson/" 
# ls ~/olson/glacier-data/rg-wasatch/timp/tf_series
# cp -a ~/olson/glacier-data/rg-wasatch/timp/tf_series /uufs/chpc.utah.edu/common/home/u1037042/data/timp/rad_results


# TO DO
# Add gauss filter to DEM
# filter shade

#######################################################################
# # FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################################################################

###############################################
# 0. DEM-BOUNDS
# set initial area and conditions (based on DEM)
topo.bounds <- function(dpath, vpath = NULL, subset=TRUE){
  #
  # returns
  # dem, viewf, dpoly
  # read in DEM
  dem <<- raster(dpath)

  if(subset){
    # 
    # Cottonwoods and Timpanogos
    # pts = cbind(445079.2, 4487065) # subset DEM
    # pt = SpatialPoints(coords = pts, proj4string = crs(dem))
    # b1 = buffer(pt, 2.1e4)
    # dem <<- crop(dem, b1)
    
    # Only Timpanogos
    # pts = cbind(445173.4, 4471190) # subset DEM to Timp
    # pt = SpatialPoints(coords = pts, proj4string = crs(dem))
    # # b1 = buffer(pt, 2e3)
    # b1 = buffer(pt, 5.4e3)
    # dem <<- crop(dem, b1)
    
    # Central Wasatch
    # pts = cbind(445079.2, 4479300) # subset DEM
    # pt = SpatialPoints(coords = pts, proj4string = crs(dem))
    # b1 = buffer(pt, 3.1e4)
    # dem <<- crop(dem, b1)
    
    # Southern Wasatch
    # b1 <- extent(c(412372.6, 465152.6, 4403026, 4448300))
    # dem <<- crop(dem, b1)
    
    # Northern Wasatch 1
    b1 <- extent(c(412372.6, 465152.6, 4510300, 4567212))
    dem <<- crop(dem, b1)
    
    # Northern Wasatch 2
    # b1 <- extent(c(412372.6, 465152.6, 4567212, 4617186))
    # dem <<- crop(dem, b1)
    
    # trim to NA values in scene
    non_na_cells <- which(!is.na(values(dem)))
    coords <- xyFromCell(dem, non_na_cells)
    longitudes <- coords[, 1]
    min_longitude <- min(longitudes); max_longitude <- max(longitudes)
    latitudes <- coords[, 2]
    min_latitudes <- min(latitudes); max_latitudes <- max(latitudes)
    b2 = extent(c(min_longitude, max_longitude, min_latitudes, max_latitudes)) # b1@bbox[2, 1], b1@bbox[2, 2])) # xmin, xmax, ymin, ymax
    dem <<- crop(dem, b2)
    
  }
  
  # BOUNDARY - for computational spatial trasformations
  dpoly <<- as(extent(dem), 'SpatialPolygons')
  crs(dpoly) <<- crs(dem)
  
  # slope and aspect
  s_a <<- slope.aspect(dem)
  
  # view factor
  dl <<-  res(dem)[1]
  d_mat <<- as.matrix(dem)
  if(!is.null(vpath)){ # create view factor if none supplied
    view_rast <<- raster(vpath)
    if(subset){
      view_rast <<- crop(view_rast, dem)
    }
    viewf <<- as.matrix(view_rast) 
  }else{viewf <<- view.factor(dem_mat = d_mat, dem = dem, dem_res = dl) }
  
  # non-NA values
  # dem_NA_idx <<- !is.na(values(dem))
  # dem_NA_len <<- sum(dem_NA_idx)
}


###############################################
# CLEARSKY - DAILY
# model insol with grid info
clear.sky.sw.topo <- function(dem, dpoly, local_seq, save_path=NULL, single_day = FALSE){
  # model sw
  # returns:
  # zenith, Ib0, Id0, sh0, cosi 
  #   s - 'dem', 'viewf', 'Ib0', 'Id0', 'cos_inc', 'shade', 'Dh', 'Ib_terrain', 'Id_terrain', 'dswtotal'
  #
  # get coordinates
  if(gsub(".*proj=([a-z]+)\\s.*", '\\1', crs(dpoly) ) != "longlat"){
    pdpoly <- spTransform(dpoly, "EPSG:4326") 
  }
  lat_lon <- rowMeans(pdpoly@bbox)[2:1]
  # ddproj = projectRaster(dem, crs = crs("EPSG:4326"))
  # lat_lon <- c(round((ddproj@extent@ymax + ddproj@extent@ymin)/2,5), round((ddproj@extent@xmax + ddproj@extent@xmin)/2,5))
  # sun position and vector
  jd <- JD(local_seq)
  # jd <-  JD(dtime)
  # timezone
  # zchng = dtime; attr(zchng, "tzone") <- "US/Mountain"
  # tmz = as.numeric(format(zchng, "%H")) - as.numeric(format(dtime, "%H"))
  # tmz = -6
  # sunvector
  # ADD CHECK FOR SUNDOWN VALUES ?
  # daylight hours (zenith <= 90)
  # sp <-  ifelse(sp1[2]>=90, 90, sp1[2]) # sp1[which(sp1[,2]<=90),] 
  # sv <- sv[which(sp1[,2]<=90),]
  # sv <- c(1,-0,-0)

  # for moments in day
  
  # jd = JD(dtime)
  # JD(jd, inverse=T)
  
  # sun position and vector
  # NO TIMEZONE NEEDED DUE TO LOCAL TIME
  sv = sunvector(jd,lat_lon[1],lat_lon[2],0); sp1=sunpos(sv)
  # verified https://www.sunearthtools.com/dp/tools/pos_sun.php
  # data.frame(hour = as.character(format(dtime, "%H:%M:%S")), elevation = 90-sp1[,2], azimuth = sp1[,1])
  # daylight hours (zenith <= 90)
  sp=sp1[which(sp1[,2]<=90),]
  sv=sv[which(sp1[,2]<=90),]
  
  dtime_final = local_seq[which(sp1[,2]<=90)]
  
  zenith <<- sp[,2]
  az_noon = which.min(abs(180-sunpos(sv))[,1])
  azimuth_eq = c(sunpos(sv)[,1][1:az_noon]-180,sunpos(sv)[,1][(az_noon+1):length(sunpos(sv)[,1])]-180)
  
  # # SAVE empty rays
  insol_norm <- matrix(NA, nrow = ncell(dem), ncol = length(zenith))
  insol_topo <- matrix(NA, nrow = ncell(dem), ncol = length(zenith))
  shade_arr <- matrix(NA, nrow = ncell(dem), ncol = length(zenith))
  
  # EMPTY rays (no NA)
  # insol_norm <- matrix(NA, nrow = dem_NA_len, ncol = length(zenith))
  # insol_topo <- matrix(NA, nrow = dem_NA_len, ncol = length(zenith))
  # shade_arr <- matrix(NA, nrow = dem_NA_len, ncol = length(zenith))
  
  # LOOP
  for (i in 1:length(zenith)){
    # print(paste(i, "of", length(zenith)))
    # zenith and incident angles
    cos_inc <- cos.slope(zenith[i], azimuth_eq[i], aspect = s_a[[2]], slope = s_a[[1]], as_mat = TRUE)
    cos_sfc <- cos(radians(zenith[i]))
    # plot(make.raster(cos_inc, dem) )
    
    # test
    # hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
    # plot(hsh, col=grey(c(0:100)/100) )
    # for (i in 1:length(zenith)){
    #   cos_inc <- cos.slope(zenith[i], azimuth_eq[i], aspect = s_a[[2]], slope = s_a[[1]])
    #   sh0 <- doshade(d_mat, sv[i,], dl=30)
    #   plot(hsh, col=grey(c(0:100)/100), main = paste("zenith:", round(zenith[i]), "| time:", JD(jd[i], inverse=T)), legend=F)
    #   plot(cos_inc, add = T, col = terrain.colors(100, 0.4))
    #   plot(make.raster(sh0, dem), add=T, col=c("black",NA), legend=F)
    # }
    # 
    # shh <- doshade(dem, sp[5,])
    
    # specify
    # zenith <<-  sp1[2]
    # azim <-  sp1[1]
    # azimuth_eq <-  azim -180
    # run shade
    sh <- doshade(d_mat, sv[i,], dl=dl) 
    # sh0 <<- make.raster(sh, dem) 
    # zenith and incident angles
    # insolation at hour
    if(!exists('visibility')){insol.params()}
    # with HEIGHT
    Idirdif <-  insolation(zenith[i],jd[i],height,visibility,RH,tempK,O3,alphag)
    Ib <-  matrix(Idirdif[,1],nrow=nrow(dem),ncol=ncol(dem), byrow=T)
    Id <-  matrix(Idirdif[,2],nrow=nrow(dem),ncol=ncol(dem), byrow=T)
    # layers
    #
    # final output - TOPOGRAPHIC SOLAR MODELS
    sw0 <<- Ib * cos_sfc + Id * cos_sfc
    sw1 <- Ib * cos_inc * sh # direct beam component
    sw2 <- Id * cos_sfc * viewf # diffuse component
    sw3 <<- sw1 + sw2 # direct and diffuse
    
    # add to empty dataframe
    insol_norm[, i] <- sw0
    insol_topo[, i] <- sw3
    shade_arr[, i] <- !sh
    
    # add to no-NA dataframe
    # insol_norm[, i] <- as.array(sw0)[dem_NA_idx]
    # insol_topo[, i] <- as.array(sw3)[dem_NA_idx]
    # shade_arr[, i] <- as.array(!sh)[dem_NA_idx]
  }
  
  # CALC topographic factor
  day_norm_sum <-  rowSums(insol_norm)
  day_topo_sum <-  rowSums(insol_topo)
  day_shade_sum <-  rowSums(shade_arr) / 4   # hours of shade
  day_topo_factor <- day_topo_sum / day_norm_sum
  day_topo_f1 <- vnormalize_1(day_topo_factor)
  
  # STACK AND SAVE
  if(!is.null(save_path)){
    # stack
    small_grid = TRUE
    topo_stk <- stack(array.to.raster(day_topo_sum, dem, small=small_grid),
                      array.to.raster(day_topo_factor, dem, small=small_grid),
                      array.to.raster(day_topo_f1, dem, small=small_grid),
                      array.to.raster(day_shade_sum, dem, small=small_grid) )
    names(topo_stk) <- c("tot_insol", "topo_factor", "topo_norm", "shade_hrs")
    topo_stk <- mask(topo_stk, dem)
    # name
    dmoment_char <- paste0("tf_day_", gsub('-','_',as.Date(JD(jd[1], inverse=T))))
    sname = paste0("clearsky_", dmoment_char,  ftype)
    # save
    writeRaster(topo_stk, file.path(save_path, "tf_day", sname))
  }
  
  # RETURN
  if(single_day){
    if(!exists("topo_stk")){
      # stack and return
      small_grid = TRUE
      topo_stk <- stack(array.to.raster(day_topo_sum, dem, small=small_grid),
                        array.to.raster(day_topo_factor, dem, small=small_grid),
                        array.to.raster(day_topo_f1, dem, small=small_grid),
                        array.to.raster(day_shade_sum, dem, small=small_grid) )
      names(topo_stk) <- c("tot_insol", "topo_factor", "topo_norm", "shade_hrs")
      topo_stk <- mask(topo_stk, dem)
    }
    return(topo_stk)
  }else{ # MONLTHLY OR SEASON
    # continue to build monthly total
    return(cbind(day_norm_sum, day_topo_sum, day_shade_sum))
  }
}


###############################################
# 6. TOPO-DAY TIMEWRAP
# execute over single day
topo.sw.day.wrapper <- function(day_string = "2021-03-21", dpath, save_path = NULL, save_day = NULL, single_day = FALSE){
  #
  # run clearsky topo for single day
  # return: tot_insol, topo_factor (0, 1), topo_norm (-1, 1)
  #
  # start
  # read in dpoly
  if(!exists('dpoly')){topo.bounds(dpath)} # returns: dem, viewf, slp, dpoly
  dtime0 <- as.POSIXct(paste0(day_string, " 00:00:00"), tz = "UTC")
  # get coordinates
  if(gsub(".*proj=([a-z]+)\\s.*", '\\1', crs(dpoly) ) != "longlat"){
    pdpoly <- spTransform(dpoly, "EPSG:4326") 
  }
  lat_lon <- rowMeans(pdpoly@bbox)[2:1]
  # SUNRISE AND SET
  # April 1st 2022 Sunrise 6:52 AM MT (13:52 UTC) - Sunset 7:32 PM (19:32 MT or 2:32 UTC +1 day) 
  # MST is UTC-07:00 | MDT is UTC-06:00)
  # MST to MDT at 2 am MST to 3 am MDT on the second Sunday in March and returns at 2 am MDT to 1 am MST on the first Sunday in November
  # March 13 2am is switch
  # close but not as accurate
  # srs2 <- daylength(lat_lon[1], lat_lon[2], JD(dtime), -6) # jdate # JD(as.Date(dtime))
  # better (UTC) & MT
  sunr <- maptools::sunriset(cbind(lat_lon[2], lat_lon[1]), dtime0, direction = "sunrise", POSIXct.out=TRUE)
  suns <- maptools::sunriset(cbind(lat_lon[2], lat_lon[1]), dtime0, direction = "sunset", POSIXct.out=TRUE)
  sr_mt = sunr$time; attr(sr_mt, "tzone") <- "US/Mountain" # .POSIXct(sunr$time, tz="US/Mountain")
  ss_mt = suns$time; attr(ss_mt, "tzone") <- "US/Mountain"
  # find nearest hour
  sun_start = lubridate::floor_date(sunr$time, unit = 'hours')
  sun_end = lubridate::ceiling_date(suns$time, unit = 'hours')
  local_start = lubridate::floor_date(sr_mt, unit = 'hours')
  local_end = lubridate::ceiling_date(ss_mt, unit = 'hours')
  
  # print(paste(" UTC: ", sunr$time, " - ", suns$time))
  # print(paste("MT: ", sr_mt, " - ", ss_mt))
  
  # format(local_start, "%H"):format(local_end, "%H")
  local_seq = seq(local_start, local_end, by = "15 mins")
  
  # check timechange
  tmz = as.numeric(format(local_start, "%H")) - as.numeric(format(sun_start, "%H"))
  
  # run for day
  day_mat <- clear.sky.sw.topo(dem, dpoly, local_seq, save_path=save_day, single_day = single_day)
  
  return(day_mat)
  
}

###############################################
# TOPOT-MULTIDAY TIMEWRAP
# execute over single day - Keep to month?
topo.sw.multiday.wrapper <- function(month = "03", year = "2021", dpath, start_day = NULL, end_day = NULL, 
                                  save_path = NULL, save_day = NULL, return_grid = TRUE){
  #
  # specify start and end days as strings
  #
  # start
  strt = Sys.time()
  
  # month time series
  start_date <- as.Date(paste0("01", month, year), "%d%m%Y")
  end_date <- seq(start_date, length.out = 2, by = "month")[2]
  dates_ls <- seq(start_date, end_date - 1, by  = "day")
  
  # fi specific dates specified
  if(!is.null(start_day)){
    # start_day = "2021-05-01" ; end_day = "2021-09-30" # MJJAS
    # start_day = "2021-10-01" ; end_day = "2022-04-30" # ODJFMA
    dates_ls <- seq(as.Date(start_day), as.Date(end_day), by="day")
  }
  # check length for computation
  if((length(dates_ls)>32) ){stop("! Periods beyond a month are currently discouraged")} # & !is.null(save_day)
  
  # SAVE empty rays (OLD)
  insol_norm_month <- matrix(NA, nrow = ncell(dem), ncol = length(dates_ls))
  insol_topo_month <- matrix(NA, nrow = ncell(dem), ncol = length(dates_ls))
  shade_arr_month <- matrix(NA, nrow = ncell(dem), ncol = length(dates_ls))
  
  # EMPTY rays (no NA)
  # insol_norm_month <- matrix(NA, nrow = dem_NA_len, ncol = length(dates_ls))
  # insol_topo_month <- matrix(NA, nrow = dem_NA_len, ncol = length(dates_ls))
  # shade_arr_month <- matrix(NA, nrow = dem_NA_len, ncol = length(dates_ls))
  
  for(dl in 1:length(dates_ls)){
    print(paste("...processing day", dl, "of", length(dates_ls), " | ", 
                dates_ls[dl]))
    
    # run current day
    current_day = as.character(dates_ls[dl])
    day_mat <- topo.sw.day.wrapper(day_string = current_day, dpath, save_path = save_path, save_day = save_day)
    
    
    # add to empty dataframe
    insol_norm_month[, dl] <- day_mat[, 1]
    insol_topo_month[, dl] <- day_mat[, 2]
    shade_arr_month[, dl] <- day_mat[, 3]
    
    # new 
    # insol_norm_month[, dl] <- day_mat[dem_NA_idx, 1]
    # insol_topo_month[, dl] <- day_mat[dem_NA_idx]
    # shade_arr_month[, dl] <- day_mat[dem_NA_idx, 3]
    
    # remove main vars
    # rm(list=setdiff(ls(), c('except these...') ) )
  }
  
  # CALC topographic factor
  day_norm_sum <-  rowSums(insol_norm_month)
  day_topo_sum <-  rowSums(insol_topo_month)
  day_shade_sum <-  rowSums(shade_arr_month) / 4   # hours of shade
  day_topo_factor <- day_topo_sum / day_norm_sum
  day_topo_f1 <- vnormalize_1(day_topo_factor)
  
  # stack
  small_grid = FALSE
  topo_stk <- stack(array.to.raster(day_norm_sum, dem, small=small_grid),
                    array.to.raster(day_topo_sum, dem, small=small_grid),
                    array.to.raster(day_topo_factor, dem, small=small_grid),
                    array.to.raster(day_topo_f1, dem, small=small_grid),
                    array.to.raster(day_shade_sum, dem, small=small_grid) )
  names(topo_stk) <- c("norm_insol", "tot_insol", "topo_factor", "topo_norm", "shade_hrs")
  topo_stk <- mask(topo_stk, dem)
  
  # ONLY works for smaller dataset
  # dft0 <- t(matrix(day_topo_sum,nrow=nrow(dem),ncol=ncol(dem), byrow = T))
  # dft1 <- t(matrix(day_topo_factor,nrow=nrow(dem),ncol=ncol(dem), byrow = T))
  # dft2 <- t(matrix(day_topo_f1,nrow=nrow(dem),ncol=ncol(dem), byrow = T))
  # dft3 <- t(matrix(day_shade_sum,nrow=nrow(dem),ncol=ncol(dem), byrow = T))
  # topo_stk <- stack(make.raster(dft0, dem), make.raster(dft1, dem), make.raster(dft2, dem), make.raster(dft3, dem))
  # names(topo_stk) <- c("tot_insol", "topo_factor", "topo_norm", "shade_sun_hrs")
  # topo_stk <- mask(topo_stk, dem)
  
  # STACK AND SAVE
  if(!is.null(save_path)){
    # name
    # dmoment_char <- paste0("tf_series_", gsub('-','_', paste0(dates_ls[1], "__", dates_ls[length(dates_ls)])))
    dmoment_char <- paste0("tf_series_m", format(dates_ls[1], "%m"))
    sname = paste0("clearsky_", dmoment_char,  ftype)
    # save
    writeRaster(topo_stk, file.path(save_path, sub_folder, sname))
  }
  
  # end time
  print( paste("Calculated topo factor for month:", dates_ls[1], "to", dates_ls[length(dates_ls)]) )
  Sys.time() - strt 
  
  # RETURN
  if(return_grid){
    return(topo_stk)
  }
  
}


###############################################
# MONTHLY TIMEWRAP
# save monthly values
topo.sw.monthly.save.wrapper <- function(start_month = NULL, end_month = NULL, year = "2021", dpath, 
                                       save_path, save_day = NULL){
  #
  # specify start and end days as strings
  #
  # start
  strt = Sys.time()
  
  # dates list
  if(!is.null(start_month)){
    months_ls <- start_month:end_month
  }else{
    months_ls <- 1:12
  }

  for(ml in 1:length(months_ls)){
    # run current month
    if(months_ls[ml]< 10){
      mostr = paste0("0", months_ls[ml])
    }else{mostr = months_ls[ml]}
    topo.sw.multiday.wrapper(month = mostr, year = year, dpath, start_day = NULL, end_day = NULL, 
                             save_path = save_path, save_day = save_day, return_grid = FALSE)
    
    # remove main vars
    # rm(list=setdiff(ls(), c('???', '?') ) )
  }
  
  print( paste("Completed months:", months_ls[1], "to", months_ls[length(months_ls)]) )
  # end time
  print( Sys.time() - strt )
}



# # # # # # # # # # # # # # # # # # # # # # # #
# Additional HELPER Functions # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # #

###############################################
exo_irr <- function(day_of_year) {
  # Solar constant (W/m^2)
  I_sc <- 1361
  # Calculate extraterrestrial solar irradiance using the formula
  I_0 <- I_sc * (1 + 0.034 * cos((2 * pi * day_of_year) / 365))
  return(I_0)
}

###############################################
insol.params <- function(use_hrrr=FALSE){
  visibility <<- 25 
  O3 <<- 340/1000 # DU = 0.01 mm (CO 325-340 DU)
  # https://disc.gsfc.nasa.gov/datasets?page=1&source=AURA%20OMI
  alphag <<- 0.35
  height <<- array(dem)
  if(use_hrrr){
    # print("! currently unable to solve for rh and temp at sub grid-level")
    # stop()
    # RH <- array(h[['rh']])
    RH <<- array(resample(h[['rh']], dem))
    # tempK <- array(h[['tempC']]) + 273.15 #278.15 # 5 C or 41 F
    tempK <<- array(resample(h[['tempC']], dem)) + 273.15
  }else{
    RH <<- 40 # 60
    tempK <<- 273.15 #278.15 # 5 C or 41 F
  }
}

###############################################
make.raster <- function(matrix, dem){
  raster(matrix,
         xmn=dem@extent@xmin, xmx=dem@extent@xmax,
         ymn=dem@extent@ymin, ymx=dem@extent@ymax, 
         crs=crs(dem))
}

###############################################
array.to.raster <- function(array, dem, small=FALSE){
  # for use in summary columns
  if(small){
    ar_mat <- t(matrix(array,nrow=nrow(dem),ncol=ncol(dem), byrow = T))
  }else{
    ar_mat <- matrix(array,nrow=nrow(dem),ncol=ncol(dem), byrow = F)
  }
  return(make.raster(ar_mat, dem))
}

###############################################
array.to.raster2 <- function(array, dem, small=FALSE){
  # adjusts for raster NA values
  array_d = as.array(dem); array_d[which(dem_NA_idx)] = array  
  if(small){
    ar_mat <- t(matrix(array_d,nrow=nrow(dem),ncol=ncol(dem), byrow = T))
  }else{
    ar_mat <- matrix(array_d,nrow=nrow(dem),ncol=ncol(dem), byrow = F)
  }
  return(make.raster(ar_mat, dem))
}

###############################################
slope.aspect <- function(dem, units = 'radians', neighbor = 8){
  s <- terrain(dem, opt='slope',unit=units,neighbors=neighbor)
  a <- terrain(dem, opt='aspect',unit=units,neighbors=neighbor)
  stk <- stack(s,a)
  return(stk)
}

###############################################
cos.slope <- function(zenith, azimuth_eq, aspect, slope, as_mat=FALSE){
  # returns a matrix of the cosine of the incident angle at a given moment
  exposures  = aspect - radians(180)
  cos_inc = acos((cos(radians(zenith)) * cos(slope)) +
                   (sin(radians(zenith)) * sin(slope) * cos(radians(azimuth_eq) - exposures)))
  if(as_mat){cos_inc = as.matrix(cos_inc)}
  # get rid of self shading values
  cos_inc[cos_inc > radians(90)] = radians(90)
  cos_inc = cos(cos_inc)
  return(cos_inc)
}

###############################################
view.factor <- function(dem_mat, dem, dem_res, elv_interval = 5, az_interval = 15){
  print("____Generating sky view factor_____")
  ptm <- proc.time()
  ELV = rev(seq(0, 90, elv_interval))
  AZI = seq(0, 345, az_interval)
  AZ = matrix(0,nrow=dim(dem_mat)[1]*dim(dem_mat)[2],ncol=length(AZI))
  for (vv in 1:length(AZI)){
    print(paste(vv, "of", length(AZI)))
    Z1 = matrix(0,nrow=dim(dem_mat)[1]*dim(dem_mat)[2],ncol=length(ELV))
    for (mm in 1:length(ELV)){
      sv = normalvector(ELV[mm],AZI[vv])
      sh <- doshade(dem_mat, sv, dl=dem_res)
      Z1[,mm] = as.array(sh)
    }
    AZ[,vv] = rowSums(Z1)/length(ELV)
  }
  VF_mat = matrix(rowMeans(AZ), nrow=nrow(dem), ncol=ncol(dem))
  #VF_dem <- make.raster(VF_mat, dem)
  print(proc.time() - ptm)
  return(VF_mat)
}

###############################################
vnormalize_1 <- function(data_vector, a = -1, b = 1){
  # normalize vector between range a and b
  dmin <- min(data_vector, na.rm = T)
  dmax <- max(data_vector, na.rm = T)
  (b - a) * ( (data_vector - dmin) / ( dmax - dmin ) ) + a
  
}

###############################################
rnormalize_1 <- function(data_raster, a = -1, b = 1){
  # normalize vector between range a and b
  dmin <- cellStats(data_raster, min, na.rm = T)
  dmax <- cellStats(data_raster, max, na.rm = T)
  (b - a) * ( (data_raster - dmin) / ( dmax - dmin ) ) + a
}

###############################################
gauss.window <- function(sigma=2, n=5) { #(spatialEco)
  m <- matrix(ncol=n, nrow=n)
  mcol <- rep(1:n, n)
  mrow <- rep(1:n, each=n)
  x <- mcol - ceiling(n/2)
  y <- mrow - ceiling(n/2)
  m[cbind(mrow, mcol)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
  m / sum(m)
}

###############################################
filter.dem <- function(dem, sigma = 2, wn = 5, funct = mean, methd = 'bilinear'){
  # low pass gaussian filter to prevent aliasing
  gm <- gauss.window(sigma=sigma, n=wn)
  smooth.dem <- focal(dem, w = gm, fun = funct, na.rm=TRUE, pad=FALSE)
  return(smooth.dem)
}

###############################################
filter.shade <- function(sh, dem, shade_size_filter = 40){
  ## SHADE
  # ? convert to 
  shd <- make.raster(sh, dem)
  shd2 = shd;shd2[is.na(dem)] = NA
  #
  # filter pixels based on size
  shclump <-  clump(!shd2, directions=8)
  f <-as.data.frame(freq(shclump))
  exludeShade <- f$value[which(f$count <= shade_size_filter)]
  shfilter <- shclump
  shfilter[shclump %in% exludeShade] <- NA 
  shf <- (!is.na(shfilter))*!is.na(d)
  shdd <- !shf
  return(as.matrix(shdd))
}


# # # # # # # # #
# POST-functions

###############################################
seasonal_vals3_5 <- function(seasonal_path, s1 = 3, s2 = 4, s3 = 5, s4 = NULL, s5 = NULL){
  #
  # default: MAM (others: DJF, MAM, JJA, SON, MJJAS)
  # combine monthly data into seasonal information
  # at minimum 3 seasons must be specified
  # up to 5 seasons may be provided
  #
  # optional seasons
  if(!is.null(s4)){
    s4_norm <- stack(seasonal_path[s4])[[1]]; s4_topo <- stack(seasonal_path[s4])[[2]]; s4_shade <- stack(seasonal_path[s4])[[5]]
  }else{s4_norm <- 0; s4_topo <- 0; s4_shade <- 0}
  # optional seasons
  if(!is.null(s5)){
    s5_norm <- stack(seasonal_path[s5])[[1]]; s5_topo <- stack(seasonal_path[s5])[[2]]; s5_shade <- stack(seasonal_path[s5])[[5]]
  }else{s5_norm <- 0; s5_topo <- 0; s5_shade <- 0}
  #
  # CALCULATE SEASONAL VALUES
  season_norm = stack(seasonal_path[s1])[[1]] + stack(seasonal_path[s2])[[1]] + stack(seasonal_path[s3])[[1]] + s4_norm + s5_norm
  season_topo = stack(seasonal_path[s1])[[2]] + stack(seasonal_path[s2])[[2]] + stack(seasonal_path[s3])[[2]] + s4_topo + s5_topo
  season_shade = stack(seasonal_path[s1])[[5]] + stack(seasonal_path[s2])[[5]] + stack(seasonal_path[s3])[[5]] + s4_shade + s5_shade
  # final values
  seas_shade_sum <-  season_shade / 4   # hours of shade
  seas_topo_factor <- season_topo / season_norm
  seas_topo_f1 <- rnormalize_1(seas_topo_factor)
  # final
  new_stk <- stack(season_topo, seas_topo_factor, seas_topo_f1, seas_shade_sum)
  names(new_stk) <- c("tot_insol", "topo_factor", "topo_norm", "shade_hrs")
  return(new_stk)
}

###############################################
# ISSUES
# Error in .rasterObjectFromFile(x, objecttype = "RasterBrick", ...) : 
#   Cannot create a RasterLayer object from this file.
seasonal_wrap_standard <- function(fpath, sub_folder){
  #
  # standard months: DJF, MAM, JJA, SON, MJJAS
  # run and save monthly data
  #
  s_files <- list.files(fpath, full.names = T)
  # standard seasons
  djf = seasonal_vals3_5(seasonal_path = s_files, s1=12, s2=1, s3=2)
  mam = seasonal_vals3_5(seasonal_path = s_files, s1=3, s2=4, s3=5)
  jja = seasonal_vals3_5(seasonal_path = s_files, s1=6, s2=7, s3=8)
  son = seasonal_vals3_5(seasonal_path = s_files, s1=9, s2=10, s3=11)
  mjjas = seasonal_vals3_5(seasonal_path = s_files, s1=5, s2=6, s3=7, s4=8, s5=9)
  # SAVE
  seas_name = c("DJF","MAM","JJA","SON","MJJAS")
  ss = list(djf, mam, jja, son, mjjas)
  for(s in 1:length(seas_name)){
    writeRaster(ss[[s]], paste0(file.path(fpath, "clearsky_tf_season_"), seas_name[s], ".tif") )
    # writeRaster(son, paste0(file.path(fpath, sub_folder, "clearsky_tf_season_"), seas_name, ".tif") )
  }
  print("Seasonal files saved!")
}

###############################################
mosaic.month.save <- function(month_list, base_path, sub_fls, save_path){
  #
  # save monthly merge
  for (m in month_list){
    #
    ml = nchar( trunc( abs(m) ) )
    if (ml<2){mchar = paste0("0", m)}else{mchar = as.character(m)}
    print(paste("Merging month", mchar))
    # merge month
    month_rast <- mosaic.month.region(mchar, base_path, sub_fls)
    names(month_rast) <- c("norm_insol", "tot_insol", "topo_factor", "topo_norm", "shade_hrs")
    #
    writeRaster(month_rast, file.path(save_path, paste0("clearsky_tf_m", mchar, ".tif")) )
  }
}

###############################################
mosaic.month.region <- function(month_char, base_path, sub_fls){
  # 
  # merge all rasters from same month
  #
  # list files for given month
  mfpths <- list.files(file.path(base_path, fls,"tf_series/"),
                       pattern = month_char, full.names = TRUE)
  # perform merge
  mr <- do.call(merge, lapply(mfpths, function(r) stack(r)) )
  # extend and name
  month_rast = extend(mr, dem)
  names(month_rast) <- c("norm_insol", "tot_insol", "topo_factor", "topo_norm", "shade_hrs")
  
  return(month_rast)
}

