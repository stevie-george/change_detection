library(rgdal)
library(raster)
library(rgeos)


### change area
magn_bkp <- raster("/home/stevie/Change_detection/KIUIC/tifs/magn_bkp_15-all_Kiuic.tif")
magnthresh <- magn_bkp /10000
magnthresh[magnthresh > -0.061] <- NA #### set threshold value here 
magn_areasieve <- areaSieve(magnthresh)

par(mfrow=c(1,1))

plot(magnthresh, main ="magnitude all deforestatiom", col="blue")
plot(magn_areasieve, main = "magnitude areasieve", col="red", add=T)

change_area_pol <- rasterToPolygons(magn_areasieve)
change_area <- gArea(change_area_pol)/10000

