## load libraries
#install.packages("rasterVis")
library(bfastSpatial)
library(raster)
library(rasterVis)
## set working directory
setwd('/home/stevie/Change_detection/FCP')

## define useful paths
path <- getwd()
path_to_tiff <- file.path(path, "tifs")
#dir.create(path_to_tiff)
path_to_sieve <- file.path(path_to_tiff, "sieves")
### useful project info  ## change this stuff
year = "00" ## two digits please
year_of_change = 2000
site = "FCP"
ext = ".tif"
ext.csv = ".csv"
ext.shp <- ".shp"


### read in raster data
rastername <- paste("magn_bkp_", year,"-all_", site, ext, sep="")
magn_bkp <- raster(paste(path_to_tiff,rastername, sep = "/"))
plot(magn_bkp)
rescaled <- magn_bkp/10000

### read in raster sieve
rastername <- paste("magn", year,"areasieve", ext, sep="")
areasieve <- raster(paste(path_to_sieve,rastername, sep = "/"))


### plot comparison
par(mfrow=c(1,2))
plot(rescaled, main = "Magnitudes by breakpoint")
plot(areasieve, main = "0.5 ha areasieve")

## produce binary mask
mask = areasieve/areasieve

## produce year map
year_map = mask * year_of_change
plot(year_map)

writeRaster(year_map, filename = paste("yearmap",year, ext, sep = ""))

### inmport agemap 
#agemapfn <- paste(path,"/agemap/",site,"_agemap",ext,sep="")
#agemap <- raster(paste(path,"/agemap/",site,"_agemap",ext,sep=""))
#plot(agemap)

#plot(agemap, col=topo.colors(21))
    

