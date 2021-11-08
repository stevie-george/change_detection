#### Change detection validation 
####  11/06/2020    (1/2)  
#### Accuracy assessment of change maps following good practices established by Oloffson et al. 2014
# Based in Devries et al. 2015 and Wagenigen University Advanced raster processing materials and
# Introduction to BfastSpatial http://www.loicdutrieux.net/bfastSpatial/
# This script loads a change raster (magn bkpt) previously produced by bfastSpatial()
# then clumps change areas into 11 pixels of change corresponding to 1 ha clusters.
      
## load libraries
library(bfastSpatial)
library(raster)
library(rgeos)
    
## set working directory
setwd('/home/stevie/Change_detection/KIUIC/')
    
## define useful paths
path <- getwd()
path_to_tiff <- file.path(path, "tifs")
#dir.create(path_to_tiff)

### useful project info  ## change this stuff
year = "20" ## two digits please
site = "Kiuic"
ext = ".tif"
ext.csv = ".csv"
ext.shp <- ".shp"

############### Option 1.  ######################################################################################
# Get results one by one -- or -- skip to Option 2. to get all results

  # # ### read in raster data
  # rastername <- paste("magn_bkp_", year,"-all_", site, ext, sep="")
  # magn_bkp <- raster(paste(path_to_tiff,rastername, sep = "/"))
  # plot(magn_bkp)
  # plot(magn_bkp, add=T, col="orange")
  # #
  # 
  # # extract and rescale magnitude and apply a -600 threshold
  # magnthresh <- magn_bkp /10000
  # plot(magnthresh)
  # magnthresh[magnthresh > -0.061] <- NA #### set threshold value here
  # # compare
  # op <- par(mfrow=c(1, 2))
  # plot(magnthresh, main="magnitude < -0.06", col ="blue")
  # plot(magn_bkp, main="magnitude all", col ="dark green")
  # 
  # ## applying sieve masks to leave out small changes
  # # two pixel sieve
  # magn_sieve <- areaSieve(magnthresh, thresh=1800)
  # 
  # # apply sieve to filter out change areas of 0.5 ha
  # magn_areasieve <- areaSieve(magnthresh)
  # plot(magn_areasieve)
  # 
  # 
  # ##### Rook's case change area calculation
  # magn_as_rook <- areaSieve(magnthresh, directions=4)
  # 
  # 
  # # compare all magn rasters
  # op <- par(mfrow=c(2, 2))
  # plot(magnthresh, main="magnitude")
  # plot(magn_sieve, main="pixel sieve")
  # plot(magn_areasieve, main="0.5ha sieve")
  # plot(magn_as_rook, main="0.5ha sieve, rook's case")
  # 
  # 
  # changeSize_queen <- clumpSize(magn_areasieve)
  # changeSize_rook <- clumpSize(magn_areasieve, directions=4)
  # # compare
  # op <- par(mfrow=c(1, 2))
  # plot(changeSize_queen, col=bpy.colors(50), main="Clump size: Queen's case")
  # plot(changeSize_rook, col=bpy.colors(50), main="Clump size: Rook's case")
  # 
  # ### choosing areasieve of half hectare
  # changeSize <- clumpSize(magn_areasieve, f=900/10000)
  # plot(changeSize, col=bpy.colors(50), main="Clump size (hectares)")
  # 
  # 
  # 
  # 
  # ### plot stats
  # changeSize <- clumpSize(magn_areasieve, f=900/10000, stats=TRUE)
  # print(changeSize$stats)
  # write.csv(changeSize$stats, file = paste("changeSize_stats", year,".csv"))
  # 
  # 
  #   ##### NEXT
  # #### IN QGIS
  # ### 1. Generate 300 Random points inside boundary layer ***
  # ### 2. Extract by location (select pixel for generating validation sites): val_poly*randompoints
  # ### 3. Generate centroids
  # ### 4. Generate 45 m pixel buffer (9 pixels)
  # ### 5. Mask magnitude layer (clip)
  # ### 6. Raster pixels to points   (extracts value)
  # ### 8. Add columns rownum as id (field calculator), xy, REF, EST, notes.
  # 
  # 
  # ##
  # formask <- paste(site,"_fnf",ext, sep = "")
  # ## read in forest mask
  # fnf <- raster(paste(path_to_tiff, formask, sep = "/"))
  # plot(fnf)
  # 
  # ## mask 0 as NA's
  # fnf[fnf == 0] <- NA
  # 
  # #
  # # ## if mask is inverted
  # # fnf[fnf == 1] <- NA
  # # fnf[fnf == 0] <- 1
  # 
  # magn_areasieve[is.na(fnf)] <- NA
  # sievetype = "areasieve"
  # outfile = paste("magn",year, sievetype, ext, sep = "")
  # 
  # writeRaster(magn_areasieve, filename=paste(path_to_tiff,"sieves",outfile, sep = "/"))
  # 
  # # #writeRaster(magn_areasieve, filename="magn_areasieve_all_15.tif")
  # # magn_areasieve<- raster(paste(path_to_tiff,"sieves",outfile, sep = "/"))
  # 
  # 
  # ### convert raster to polygon
  # pol <- rasterToPolygons(magn_areasieve)
  # #plot(pol)
  # #pol2 <-rasterToPolygons(magnthresh)
  # # gArea(pol)/10000
  # # gArea(pol2)/10000
  # 
  # ### write polygon shapefile
  # path_to_shapes <- paste(path,"shapes", sep = "/")
  # #dir.create(path_to_shapes)
  # 
  # fn <- paste("val_poly_",year,ext.shp, sep="")
  # shapefile(pol,file = paste(path_to_shapes,fn, sep = "/"))

########## Option 2.################################################################################################################ 
## use a loop to produce all change results for rasters in folder

years <- c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20")

get_results <- function(years, site, ext, ext.csv, ext.shp) {for (i in 1:length(years)){
  ### read in raster data
  rastername <- paste("magn_bkp_", years[i],"-all_", site, ext, sep="")
  magn_bkp <- raster(paste(path_to_tiff,rastername, sep = "/"))
  plot(magn_bkp)
  plot(magn_bkp, add=T, col="orange")
  # 
  # extract and rescale magnitude and apply a -600 threshold
  magnthresh <- magn_bkp /10000
  plot(magnthresh)
  magnthresh[magnthresh > -0.061] <- NA #### set threshold value here
  # compare
  op <- par(mfrow=c(1, 2))
  plot(magnthresh, main="magnitude < -0.06", col ="blue")
  plot(magn_bkp, main="magnitude all", col ="dark green")
  
  ## applying sieve masks to leave out small changes
  # two pixel sieve 
  magn_sieve <- areaSieve(magnthresh, thresh=1800)
  
  # apply sieve to filter out change areas of 0.5 ha
  magn_areasieve <- areaSieve(magnthresh)
  plot(magn_areasieve)
  
  ##### Rook's case change area calculation 
  magn_as_rook <- areaSieve(magnthresh, directions=4)
  
  # compare all magn rasters
  #op <- par(mfrow=c(2, 2))
  #plot(magnthresh, main="magnitude")
  #plot(magn_sieve, main="pixel sieve")
  #plot(magn_areasieve, main="0.5ha sieve")
  #plot(magn_as_rook, main="0.5ha sieve, rook's case")
  
  changeSize_queen <- clumpSize(magn_areasieve)
  changeSize_rook <- clumpSize(magn_areasieve, directions=4)
  # compare
  #op <- par(mfrow=c(1, 2))
  #plot(changeSize_queen, col=bpy.colors(50), main="Clump size: Queen's case")
  #plot(changeSize_rook, col=bpy.colors(50), main="Clump size: Rook's case")
  
  ### choosing areasieve of half hectare
  changeSize <- clumpSize(magn_areasieve, f=900/10000)
  #plot(changeSize, col=bpy.colors(50), main="Clump size (hectares)")
  
  ### plot stats  
  changeSize <- clumpSize(magn_areasieve, f=900/10000, stats=TRUE)
  print(changeSize$stats)
  write.csv(changeSize$stats, file = paste("changeSize_stats", years[i],".csv"))
  ## 
  formask <- paste(site,"_fnf",ext, sep = "")
  ## read in forest mask
  fnf <- raster(paste(path_to_tiff, formask, sep = "/"))
  # plot(fnf)
  ## mask 0 as NA's 
  #fnf[fnf == 0] <- NA
  # 
  # ## if mask is inverted
  fnf[fnf == 1] <- NA
  fnf[fnf == 0] <- 1
  
  magn_areasieve[is.na(fnf)] <- NA
  sievetype = "areasieve"
  outfile = paste("magn",years[i], sievetype, ext, sep = "")
  
  writeRaster(magn_areasieve, filename=paste(path_to_tiff,"sieves",outfile, sep = "/"), overwrite=T)
  
  # #writeRaster(magn_areasieve, filename="magn_areasieve_all_15.tif")
  # magn_areasieve<- raster(paste(path_to_tiff,"sieves",outfile, sep = "/"))
  
  ### convert raster to polygon
  pol <- rasterToPolygons(magn_areasieve)
  #plot(pol)
  #pol2 <-rasterToPolygons(magnthresh)
  # gArea(pol)/10000
  # gArea(pol2)/10000
  
  ### write polygon shapefile
  path_to_shapes <- paste(path,"shapes", sep = "/")
  #dir.create(path_to_shapes)
  
  fn <- paste("val_poly_",years[i],ext.shp, sep="")
  shapefile(pol,file = paste(path_to_shapes,fn, sep = "/"))
}
}

get_results(years, site='Palmar', ext, ext.csv, ext.shp)
  