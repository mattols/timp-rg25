#
# Clear-sky Radiation Model - Rock Glaciers Wasatch Front, UT
# Matt Olson
# 01-10-2025
#
# R 4.03 Geospatial Packages (CHPC)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #







# individual layers
# topofactor
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100), legend=F, axes=F, box=F,  main = "MJJAS insolation topography factor" )
# plot(tf_series[[2]], add = T, col = heat.colors(100, 0.4))
plot(mjjas[[2]], add = T, col = heat.colors(100, 0.6))

plot(hsh, col=grey(c(0:100)/100), legend=F, axes=F, box=F,  main = "MJJAS insol topography factor", mar = c(2,2, 1, 4))
plot(mjjas[[2]], add = T, legend=F, axes=F, box=F, col = heat.colors(100, 0.6))
plot(mjjas[[2]], legend.only=TRUE, col=heat.colors(100, 0.6),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=" ", side=3, font=2, line=2.5, cex=1.2))


plot(tf_series[[3]], legend=F, axes=F, box=F, col = c(grey(c(0:100)/100), heat.colors(100)),
     main = "Normalized topographic factor")
plot(tf_series[[3]], legend.only=TRUE, col=c(grey(c(0:100)/100), heat.colors(100)),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text="?", side=3, font=2, line=2.5, cex=1.2))

# normalized
plot(tf_series[[3]], legend=F, axes=F, box=F, col = c(grey(c(0:100)/100), heat.colors(100)),
     main = "Normalized topographic factor")
plot(tf_series[[3]], legend.only=TRUE, col=c(grey(c(0:100)/100), heat.colors(100)),
     legend.width=2.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text="?", side=3, font=2, line=2.5, cex=1.2))
# legend("bottomleft", "Normalized variation in insolation due to topography" )

# shade
plot(tf_series[[4]], col = grey(c(100:0)/100), main = "March solar hours of shade")

plot(mjjas[[4]], col = grey(c(100:0)/100), legend=F, axes=F, box=F, 
     main = "MJJAS solar hours of shade")
plot(mjjas[[4]], legend.only=TRUE, col=grey(c(100:0)/100),
     legend.width=1.5, legend.shrink=0.7,
     axis.args=list(cex.axis=1.2),
     legend.args=list(text=" ", side=3, font=2, line=2.5, cex=1.2))

# test plot
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100), legend=F, main = paste("TopoFactor:", as.Date(JD(jd[1], inverse=T))) )
plot(make.raster(dft1, dem), add = T, col = heat.colors(100, 0.4))
#
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100), legend=F, main = paste("TopoFactor:", as.Date(JD(jd[1], inverse=T))) )
plot(make.raster(dft2, dem), add = T, col = heat.colors(100, 0.4))



#
#
#
library(raster)
mpath <- list.files("~/olson/glacier-data/rg-wasatch/timp/tf_series", full.names = T)
plot(stack(mpath[1]), col = heat.colors(100),
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade"))  
plot(stack(mpath[9]), col = heat.colors(100),
     main = c("Total insolation", "Topographic insol factor", "Topo factor norm", "Solar hours of shade"))  


