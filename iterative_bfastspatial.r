## Iterative application of 'bfmspatial' for detecting change in Landsat timeseries
## This script was developed from the workflow proposed by Verbesselt et al., 2015 and DeVries et al-. 2015
## Stephanie P. George  e mail = stephanie.p.george@gmail.com 

# This script contains the following analyses
# 1. Preprocessing of landsat imagery to create timestack (main input for bfmspatial)
# 2. Application of bfmspatial iteratively for the length of the time series. This method applies the algorithm for a single year using all available information as the historical period.
# This script will produce: NDVI data for all timeseries, breakpoint rasters containing the year of deforestation, magnitude rasters containing the difference between observed NDVI and estimated NDVI from the harmonic model, and breakpoint-magnitude layer obtaining the magnitude for each breakpoint where observed magnitude differs from the model with 95% confidence.

# set working directory to site
setwd('/media/stevie/HD710 PRO/ESPA/BFast/KIUIC/')

# Load bfastSpatial 
library(bfastSpatial)
library(parallel)

## Make sure all directories are created at the correct locations
# set path of working directory
path <- getwd()

# Set temporary directory for raster processing
tmpDir <- rasterOptions()$tmpdir

# Directory where data is stored 
inDir <- file.path(path, 'data')

# Directory where intermediate outputs are stored
stepDir <- file.path(inDir, 'datastep')

# Directory of Landsat data (copy all data in tar.gz format, without stats and other files. Only "L0720046589956T1" files).
# landsatDir <- file.path(stepDir, 'landsat')
## directory for FCP images path/row 20/47
# directory for KIUIC and PALMAR path/row 20/46

# Directory where raw landsat timeseries is stored
landsatDir57 <- file.path("/media/stevie/HD710 PRO/ESPA/BFast/KIUIC/L57")
landsatDir8 <- file.path("/media/stevie/HD710 PRO/ESPA/BFast/KIUIC/L08") 

# Directory where NDVI imagery is stored
ndviDir <- file.path(stepDir, 'ndvi')

# Output directory
outDir <- file.path(inDir, 'out')


# 1. Preprocessing of Landsat data
# create an extent object in order to clip processing extent to 60 x 60 km window (UTM 16 N coords).
# It is extremely important for all the images to have the same extent, therefore we have to specify via 
# an extent object with four coordinates in UTM (xmin, xmax, ymin, ymax).

## This extent must correspond to the four corners of the 60 x 60 km sites in Yucatan 
# this must match the site
e <- extent(201753.52, 262685.59, 2170910.51, 2231840.28) # extent in Kaxil Kiuic
#e <- extent(159727,221807.4,2269844,2332127) # extent in Palmar 
#e <- extent(331533.8, 391974.4, 2099930, 2160368) #extent in FCP

# Warning: Due to the change in USGS ESPA file naming convention when using developers version of bfastSpatial use the following to apply the cloud mask: 

# for landsat 7 keep = c(66, 130))
# for landsat 8 keep = c(322, 386))

# Process landsat unzips Landsat files, applies cloud mask, and calculates VI when file is not present in the director.
# If the file is present in the directory, a loop will read in the file. 

# processLandsatBatch(x = landsatDir8, outdir = ndviDir,
#                     delete = TRUE, e=e, overwrite = FALSE, mask = 'pixel_qa', vi = 'ndvi',
#                     keep = c(322, 386))
# 

# Stack ndvi or process NDVI from landsat imagery 

# IMPORTANT: Run separately  for Landsat 7 and 8 and compile all output in same directory prior to stacking
if (!file.exists(file.path(inDir, 'ndvi_stack.grd'))) {
                    processLandsatBatch(x = landsatDir8, outdir = ndviDir,
                    delete = TRUE, e=e, overwrite = FALSE, mask = 'pixel_qa', vi = 'ndvi',
                    keep = c(322, 386))
                # make temporal ndvi stack
                    ndviStack <- timeStack(x = ndviDir, pattern = glob2rx('*.grd'),
                    filename = file.path(inDir, 'ndvi_stack.grd'),
                    datatype = 'INT2S')
                    } else {
                    ndviStack <- brick(file.path(inDir, 'ndvi_stack.grd'))
                    }

# if (!file.exists(file.path(inDir, 'ndvi_stack.grd'))) {
#                     processLandsatBatch(x = landsatDir57, outdir = ndviDir,
#                     delete = TRUE, e=e, overwrite = FALSE, mask = 'pixel_qa', vi = 'ndvi',
#                     keep = c(66, 130))
#                 # make temporal ndvi stack
#                     ndviStack <- timeStack(x = ndviDir, pattern = glob2rx('*.grd'),
#                     filename = file.path(inDir, 'ndvi_stack.grd'),
#                     datatype = 'INT2S')
#                     } else {
#                     ndviStack <- brick(file.path(inDir, 'ndvi_stack.grd'))
#                     }

## Prepare stack with all NDVI layers
# when running only ndviStack for 1st time use:
# ndviStack <- timeStack(x = ndviDir, pattern = glob2rx('*.grd'), filename = '/home/stevie/Change_detection/KIUIC/out/ndvi_stack.grd', datatype = 'INT2S')

## Preprocessed data exploration
# pass ndvistack to object
x <- ndviStack
# check nomenclature for landsat is correct
# show names
names(x)
s <- getSceneinfo(names(x))
s
## Observation of NDVI
# filter the scene to show pixels above a certain threshold (NDVI >.7 in this case)
checkThresh <- function(x){
  # first, get rid of NA's
  x <- x[!is.na(x)]
  # if there still values left, count how many are above the threshold
  # otherwise, return a 0
  if(length(x) > 0){
    y <- length(x[x > 7000])
  } else {
    y <- 0
  }
  # return ndvivalue
  return(y)
}

# Pass function to summarybrick
customStat <- summaryBrick(x, fun=checkThresh)
plot(customStat, main = "# of observations where NDVI > 0.7")

# Create a histogram to show available scenes per year
s$year <- as.numeric(substr(s$date, 1, 4))
hist(s$year, main="Scenes per Year", 
     xlab="year", ylab="# of scenes")

# Count valid observations
obs <- countObs(x)
plot(obs)

# Summary of valid observations
obs <- countObs(x, as.perc=TRUE)
summary(obs)
hist(obs, main ="percent obs")
# % NA per pixel 
percNA <- 100 - countObs(x, as.perc=TRUE)
hist(percNA, main = "percent NA per pixel", xlab = "% NA per pixel")
plot(percNA, main="percent NA per pixel")
summary(percNA)


## 2. Bfm implementation
# Define site according to extent
site = "kiuic"

# Make year object with timeseries extent
years <- 1999 + seq(1: 20)


## Loop for application of bfmspatial algorithm iteratively

for (i in 1:length(years)) {
  #### Apply bfmspatial to timeseries
  # file = paste("bfm_",years[i],"-harmon-",site,".grd", sep = "")
  if (!file.exists(fn <- file.path(outDir,paste("bfm_",years[i],"-harmon-",site,".grd", sep = "")
  ))) {
    # make sure the filename is correct and the folder exists
    bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                      start = years[i], history = c("all"), h=0.25, monend=years[i+1], level=0.05,
                      formula = response ~ harmon, type = "OLS-MOSUM",
                      order = 1, mc.cores = 1, returnLayers = c("breakpoint", "magnitude", "error",
                                                                "history", "r.squared", "adj.r.squared", "coefficients"),
                      filename = file.path(outDir, paste("bfm_",years[i],"-harmon-",site,".grd", sep = ""))
    )
  } else {
    bfm <- brick(fn)
  }
  
  # separate filename in object for further processing
  fn <- file.path(outDir, filename)
  
  bfm <- brick(fn)
  
  # the first layer corresponds to 'breakpoints' for the analysis year, which will correspond to 'change' location for specific point.
  change <- raster(bfm, 1)
  
  # write change layer to raster  
  writeRaster( 
    change, 
    filename=paste("bkp_",years[i],"_",site,".tif", sep = ""), 
    format="GTiff", 
    overwrite=TRUE)  
  
  # Make month object
  months <- changeMonth(change)
  # Define color labels for month object
  monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
                 "jul", "aug", "sep", "oct", "nov", "dec")
  cols <- rainbow(12)
  
  # Produce plot with estimated change pixels per month
  plot(months, col=cols, breaks=c(1:12), legend=FALSE)
  # apply legend
  legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)
  
  # Second layer corresponds to magnitude i.e. the difference between observed and predicted ndvi from the harmonic model 
  magn <- raster(bfm, 2)
  
  # Write magnitude layer to output
  writeRaster( 
    magn, 
    filename=paste("magn_",years[i],"_",site,".tif", sep = ""), 
    format="GTiff", 
    overwrite=TRUE)  
  
  # Make a magnitude-breakpoint layer containing only the magnitudes where breakpoints have been identified
  magn_bkp <- magn
  magn_bkp[is.na(change)] <- NA
  op <- par(mfrow=c(1, 2))
  plot(magn_bkp, main="Magnitudes for breakpoints")
  plot(magn, main="Magnitudes for all pixels")
  
  # Write magnitude bkpts output in separate raster
  writeRaster( 
    magn_bkp, 
    filename=paste("magn_bkp_",years[i],"_",site,".tif", sep = ""), 
    format="GTiff", 
    overwrite=TRUE)
}








