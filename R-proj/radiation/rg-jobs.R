


# call functions
rm(list=ls())
source('rg-radiation-func-1.R', echo=FALSE)

# funs
topo.bounds(dpath)
topo.sw.monthly.save.wrapper(dpath = dpath, save_path = timp_storage)