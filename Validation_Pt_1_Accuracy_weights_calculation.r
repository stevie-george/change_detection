library(rgdal)
library(raster)
library(rgeos)
  
  
### change area
magn_bkp <- raster("/home/stevie/Change_detection/PALMAR/tifs/magn_bkp_08-all_Palmar.tif")
plot(magn_bkp)
magnthresh <- magn_bkp /10000
plot(magnthresh)
magnthresh[magnthresh > -0.061] <- NA #### set threshold value here

#magnthresh <- raster('/home/stevie/Change_detection/KIUIC/tifs/sieves/magn15areasieve.tif')
#plot(magnthresh)

change_area_pol <- rasterToPolygons(magnthresh)
change_area <- gArea(change_area_pol)/10000 # layer in m, convert to ha
  
total_pol <- readOGR("/home/stevie/Change_detection/PALMAR/shapes_old/shapes/Palmar_square_poly6_01.shp") # Palmar
#total_pol <- readOGR("/home/stevie/Change_detection/FCP/shapes/FCP_square_poly.shp") # FCP
#total_pol <- readOGR("/home/stevie/Change_detection/KIUIC/shapes_old/Kiuic_square_poly_01_utm16.shp") # Kiuic

total_area <- gArea(total_pol)/10000 # layer in m, convert to ha
forest <- raster('/home/stevie/Change_detection/PALMAR/tifs/Palmar_fnf.tif')
plot(forest)
forest[forest==0] <- NA
#plot(forest)

no_change_area <- round(total_area - change_area)
# print results
print(paste("total change area =", change_area, "ha"))
print(paste("total no change area =",no_change_area, "ha"))
print(paste("total area =",round(total_area, digits = 2), "ha"))
#print(paste("total forest area =",round(forest_area, digits = 2), "ha"))
  
  
# calculate weights
wi_c <- change_area/total_area
  
wi_nc <-no_change_area/total_area
  
wi_c + wi_nc
  
cochrans <- function(wi_c,  si_c,wi_nc, si_nc, so) {
  (((wi_c*si_c) + (wi_nc*si_nc))/so)^2
}
    
si_c  = sqrt(.90*(1-.90))
si_nc = sqrt(.95*(1-.95))
so = .01  
  
N = cochrans(wi_c, si_c, wi_nc, si_nc, so)
N
  
#---------------- Get samples in change and no change area  ----------------------     
# plot(forest, add=T)
# 
# change_area_year <- raster("/home/stevie/Change_detection/KIUIC/tifs/sieves/magn15areasieve.tif")
# 
# # change_cropped <- crop(change_area_year, total_pol)
# 
# forest <- raster('/home/stevie/Change_detection/KIUIC/tifs/KIUIC_fnf.tif')
# 
# # forest_cropped <- crop(forest, total_pol)
# 
# plot(forest)
#         
# plot(change_area_year, add=T, col="dark orange")
#         
# forest[change_area_year] <- NA
#         
# forest[forest ==0] <- NA
# plot(forest)
#         
# plot(magnthresh, col ="dark orange", add=T)
#         
# sample_change <- sampleRandom(change_area_year, size=500, cells=TRUE, sp=TRUE)
#         
# sample_noch <- sampleRandom(forest, size = 700, cells = TRUE, sp=TRUE)
#         
# plot(sample_change, add =T, pch =20, col= "blue")
#         
# plot(sample_noch, add=T, pch=20, col ="yellow")
# 
# shapefile(sample_change, "/home/stevie/Change_detection/KIUIC/Val_change_KIUIC_10.shp", overwrite=T)
#  
# shapefile(sample_noch, "/home/stevie/Change_detection/KIUIC/Val_nochange_KIUIC_10.shp", overwrite=T)
# 
#   getwd()
  

# forest_proj <- projectRaster(forest, crs =crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(forest)
# 
