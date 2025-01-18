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
plot(hsh, col=grey(c(0:100)/100), legend=F, main = "March insolation topography factor" )
plot(tf_series[[2]], add = T, col = heat.colors(100, 0.4))

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


# test plot
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100), legend=F, main = paste("TopoFactor:", as.Date(JD(jd[1], inverse=T))) )
plot(make.raster(dft1, dem), add = T, col = heat.colors(100, 0.4))
#
hsh = hillShade(s_a[[1]], s_a[[2]], angle = 25, direction = 145)
plot(hsh, col=grey(c(0:100)/100), legend=F, main = paste("TopoFactor:", as.Date(JD(jd[1], inverse=T))) )
plot(make.raster(dft2, dem), add = T, col = heat.colors(100, 0.4))
