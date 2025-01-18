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


# 15th of every month


# integrate full melt season



## THOUGHTS
# outputs:
# - Days of sun per year
# - 




