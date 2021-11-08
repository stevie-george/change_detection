library(rgdal)
library(raster)
library(rgeos)
  
  
### change area
change_area_year <- raster("/home/stevie/Change_detection/KIUIC/yearmap15.tif")
plot(change_area_year, col ="dark orange")
change_area_pol <- readOGR("/home/stevie/Change_detection/KIUIC/shapes/val_poly_15.shp")
change_area <- gArea(change_area_pol)/10000 # layer in m, convert to ha
  
total_pol <- readOGR("/home/stevie/Change_detection/KIUIC/shapes/Kiuic_square_poly_01_utm16.shp")
total_area <- gArea(total_pol)/10000 # layer in m, convert to ha
  #roads <- readOGR("/home/stevie/Change_detection/shapes/roads/conjunto_de_datos/red_vial.shp")
  # forest <- raster('/home/stevie/Change_detection/KIUIC/tifs/Kiuic_fnf.ti f')
  # forest_all <- raster('/home/stevie/Change_detection/KIUIC/tifs/Kiuic_fnf.tif')
  # plot(forest)
  # forest[forest==0] <- NA
  # forest_proj <- projectRaster(forest, crs =crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # plot(forest)
  # 
  # #get sizes of all cells in raster [km2]
  # cell_size<-area(forest_proj, na.rm=TRUE, weights=FALSE)
  # #delete NAs from vector of all raster cells
  # ##NAs lie outside of the rastered region, can thus be omitted
  # cell_size<-cell_size[!is.na(cell_size)]
  # 
  # #compute area [km2] of all cells in geo_raster
  # raster_area<-round((length(cell_size)*median(cell_size))*100, digits = 2)
  # 
  # ## subtract change area to obtain total forest area
  # forest_area <- raster_area # layer in m, convert to ha dividing over 10000
  
forest <- readOGR('/home/stevie/kiuic_non_forest_poly.shp')
forest_area <- gArea(forest)/10000
no_change_area <- round(forest_area - change_area)
# print results
print(paste("total change area =", change_area, "ha"))
print(paste("total no change area =",no_change_area, "ha"))
print(paste("total area =",round(total_area, digits = 2), "ha"))
print(paste("total forest area =",round(forest_area, digits = 2), "ha"))
  
  
# calculate weights
wi_c <- change_area/forest_area
  
wi_nc <-no_change_area/forest_area
  
wi_c + wi_nc
  
cochrans <- function(wi_c,  si_c,wi_nc, si_nc, so) {
  (((wi_c*si_c) + (wi_nc*si_nc))/so)^2
}
  
si_c  = sqrt(.90*(1-.90))
si_nc = sqrt(.95*(1-.95))
so = .01  
  
  
N = cochrans(wi_c, si_c, wi_nc, si_nc, so)
N
  
#--------------------------------------   
# plot(change_area_year, add = T)
# forest_nochange_1 <- mask(nonforest_all, change_area_year, inverse=T)
#   
# plot(forest_nochange_1)
#   
# sample_change <- sampleRandom(change_area_year, size=300, cells=TRUE, sp=TRUE)
# sample_noch <- sampleRandom(forest_nochange_1, size = 700, cells = TRUE, sp=TRUE)
# plot(sample_change, add =T, pch =20, col= "red")
# plot(sample_noch, add=T, pch=20, col ="yellow")
# shapefile(sample_change, "Val_change.shp", overwrite=T)
# shapefile(sample_noch, "Val_nochange.shp", overwrite=T)
# getwd()

