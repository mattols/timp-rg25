

# call functions
rm(list=ls())
source('rg-radiation-func-1.R', echo=FALSE)

# save_path
# timp_storage = "~/olson/glacier-data/rg-wasatch/timp"
#
# central wasatch
# dir.create("~/olson/glacier-data/rg-wasatch/south_wasatch")
# dir.create("~/olson/glacier-data/rg-wasatch/south_wasatch/tf_series")
rg_storage <- "~/olson/glacier-data/rg-wasatch/south_wasatch"

# funs
topo.bounds(dpath, vpath)
topo.sw.monthly.save.wrapper(dpath = dpath, save_path = rg_storage)