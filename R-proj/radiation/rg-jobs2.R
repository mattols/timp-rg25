

# call functions
rm(list=ls())
source('rg-radiation-func-1.R', echo=FALSE)

# save_path
# timp_storage = "~/olson/glacier-data/rg-wasatch/timp"
#
# central wasatch
# dir.create("~/olson/glacier-data/rg-wasatch/north_wasatch1")
# dir.create("~/olson/glacier-data/rg-wasatch/north_wasatch1/tf_series")
rg_storage <- "~/olson/glacier-data/rg-wasatch/north_wasatch1"

# funs
topo.bounds(dpath, vpath)
# topo.sw.monthly.save.wrapper(dpath = dpath, save_path = rg_storage)
topo.sw.monthly.save.wrapper(start_month = 7, end_month = 12,
                             dpath = dpath, save_path = rg_storage)