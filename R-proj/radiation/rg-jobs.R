


# call functions
rm(list=ls())
source('rg-radiation-func-1.R', echo=FALSE)

# save_path
# timp_storage = "~/olson/glacier-data/rg-wasatch/timp"
#
# central wasatch
# dir.create("~/olson/glacier-data/rg-wasatch/central_wasatch")
# dir.create("~/olson/glacier-data/rg-wasatch/central_wasatch/tf_series")
rg_storage <- "~/olson/glacier-data/rg-wasatch/central_wasatch"
# ymin       : 4448316 
# ymax       : 4510296 

# funs
topo.bounds(dpath, vpath)
# topo.sw.monthly.save.wrapper(dpath = dpath, save_path = rg_storage)
topo.sw.monthly.save.wrapper(start_month = 6, end_month = 6,
                             dpath = dpath, save_path = rg_storage)

list.files("~/olson/glacier-data/rg-wasatch/central_wasatch/tf_series/")

# [1] "Completed months: 1 to 12"
# Time difference of 21.56584 mins (For TIMP extent)

# For central Wasatch
# Error in .local(.Object, ...) : Dataset copy failed
# Calls: sourceWithProgress ... <Anonymous> -> new -> initialize -> initialize -> .local
# Execution halted

# # 2 hrs 23 mins
# ?