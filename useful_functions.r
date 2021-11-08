#### Auxiliary Functions  Change detection validation 
#### Accuracy assessment of change maps following good practices established by Oloffson et al. 2014
# This script loads a change raster (magn bkpt) previously produced by bfastSpatial()
# then clumps change areas into 11 pixels of change corresponding to 1 ha clusters.

### load necessary libraries
#install.packages('spatstat')
#install.packages('igraph')
#remotes::install_github("raff-k/Lslide")

library(sf)
library(raster)
library(rgdal)
library(Lslide)

## define useful paths
path <- file.path("/home/stevie/Change_detection")
## set working directory
setwd(path)
path_to_tiff <- file.path(path, "tifs")

## set year of change to be evaluated
year <- '2015'

### define filename prefix
def<- 'DEF_nf_'

### define extension 
ext <- '.tif'

## get filename object
def_raster <- file.path(path_to_tiff, paste(def,year,ext, sep = ""))

### Generate random points inside data pixels of raster file
### load raster 
#def_raster <- file.path(path_to_tiff, paste(def,year,ext, sep = ""))

def_2016 <- raster(def_raster)



## Group raster cells into clumps based on the Rook's Case

forestclumps <- clump(def_2016, directions=4, filename=fn)
formask <- raster(def_raster)
plot(forestclumps)


clumpFreq <- freq(forestclumps)
head(clumpFreq)
tail(clumpFreq)


## Coerce freq table to data.frame
clumpFreq <- as.data.frame(clumpFreq)

## which rows of the data.frame are only represented by one cell?
str(which(clumpFreq$count<11))

## which values do these correspond to?
str(clumpFreq$value[which(clumpFreq$count<11)])

## Put these into a vector of clump ID's to be removed
excludeID <- clumpFreq$value[which(clumpFreq$count<11)]

## Make a new forest mask to be sieved
formaskSieve <- def_2016

## Assign NA to all clumps whose IDs are found in excludeID
formaskSieve[forestclumps %in% excludeID] <- NA

## Zoom in to a small extent to check the results
# Note: you can define your own zoom by using e <- drawExtent()

opar <- par(mfrow=c(1, 2)) # allow 2 plots side-by-side
plot(formask, col="dark green", legend=FALSE)
plot(formaskSieve, col="dark green", legend=FALSE)
poly <- readOGR(poly_path)    


#### IN QGIS  
### 1. Generate 300 Random points inside boundary layer ***
### 2. Extract by location (select pixel for generating validation sites)
### 3. Generate centroids
### 4. Generate 45 m pixel buffer (9 pixels)


points <- readOGR("path/to/points_generated_in_qgis.shp")

### plot resulting points over change raster
plot(def_2016)
points(points) 

##### Section UNDER CONSTRUCTION  
##### sample random points at minimum distance

#pol_sf <- st_read("/home/stevie/Change_detection/shapes/val_poly_16.shp")
#points <- genRandomPnts(pol_sf, out="samplepnts.shp", n = 100, dist = 500, cores = 4)
#st_write(points, dsn = "points.shp")
#pts <- readOGR("points.shp")
#library(sp)

#dmat <- spDists(pts)

#min.dist <- 500 
#dmat[dmat <= min.dist] <- NA
#head(pts)
#m = gDistance(points, byid = TRUE)
# #samples <- data.frame(pts$ID, kNN=NA)
# for(i in 1:nrow(dmat) ) {
#   x <- as.vector( dmat[,i] )
#   names(x) <- samples$ID
#   x <- x[!is.na(x)]
#   if(!length(x) == 0) {
#     samples[i,][2] <- names(x)[sample(1:length(x), 1)]
#   } else {
#     samples[i,][2] <- NA
#   }   
# }


### select a number of polygons for validation

# require(sp)
# select 20 random polygons:
#sample = pol[sample(length(pol),100),]
#plot(sample)
# sample a single random point inside each polygon:
#pts = do.call(rbind,
 #             lapply(sample@polygons, spsample, n=1, type="random"))
###

#coords <- coordinates(points)
#head(coords)
#df <- data.frame(points, coords)
#head(df)
#head(coords)
      
###### Under construction
### mask to non forest areas

# fnf <- raster('../tifs/Kiuic_FNF.tif')
# magn_bkp <- raster('magn_bkp_11-all_Kiuic.tif')
# plot(fnf)
# plot(magn_bkp, add =T)
# r1 <- crop(magn_bkp, extent(e))
# r2 <- crop(fnf, extent(e))
# #writeRaster(r2, filename = "kiuic_fnf.tif")
# getwd()
# fnf[fnf == 0] <- NA
# 
# plot(r2)
# plot(fnf)
# magn_bkp_1 <- mask(r1, r2, maskvalue = 1)
# magn_bkp[is.na(fnf)] <- NA
# plot(fnf)
# plot(magn_bkp, add=T)
# 
# magn_bkp = r1 * r2

# extract and rescale magnitude
#magn <- raster(bfm, 2) / 10000
# remove all non-breakpoints (we already stored this object as magn_bkp )
#magn[is.na(change)] <- NA

  