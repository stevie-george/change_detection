# Change detection : BFast implementation
# Stephanie P. George stephanie.p.george@gmail.com 31/08/2018
# Edited 27/04/2920 
# This code is based in Smith et al. 2019

# Install necessary packages (uncomment section 'ctrl + shift + c' and run)

# install.packages("snow") # for multicore processing only
# install.packages("devtools")
# ref = 'develop' # install developers version of BFast spatial
# devtools::install_github('loicdtx/bfastSpatial', ref = 'develop')

  
# For windows 
# set working directory 
# setwd("E:/ESPA/BFast/KIUIC/")

# Set directory path. This is the path where all landsat data is stored and where the following directories have to be created 
# 1) C:/ESPA/BFast/data/datastep/ 
# 2) C:/ESPA/BFast/data/out 
# 3) C:/ESPA/BFast/data/datastep/landsat 
# 4)C:/ESPA/BFast/data/datastep/ndvi

# for UNIX  change all working directories to adjust folder architecture
# setwd('/Volumes/ADATA\ HD710/')
setwd('/media/stevie/HD710 PRO/ESPA/BFast/KIUIC/')

# make sure all directories have been created in the correct order
# set working directory to the directory where "data" is stored 
# set path for reading and saving files
path <- getwd()

# load bfastSpatial 
library(bfastSpatial)
library(parallel)

# and set tmpdir
tmpDir <- rasterOptions()$tmpdir

# set the path to the location of data 
inDir <- file.path(path, 'data')

# stepDir is where intermediary outputs are stored
stepDir <- file.path(inDir, 'datastep')

# directory for Landsat data (copy all data in tar.gz format, without stats and other files. Only "L0720046589956T1" files).
# landsatDir <- file.path(stepDir, 'landsat')
# landsatDir57 <- file.path("/media/stevie/HD710 PRO/ESPA/BFast/FCP/L57")
# landsatDir8 <- file.path("/media/stevie/HD710 PRO/ESPA/BFast/FCP/L08_2020") 

# directory for KIUIC and PALMAR path/row 20/46
# directory for FCP images path/row 20/47
# landsatDir57 <- file.path("E:/ESPA/BFast/FCP/L57") 
# landsatDir8 <- file.path("E:/ESPA/BFast/KIUIC/L08")#directory for FCP
# where individual VI layers are stored prior to being stacked; ndviDir, eviDir, etc. are subdirectories of stepDir (Here we only used NDVI)
ndviDir <- file.path(stepDir, 'ndvi')
#ndviDir <- 'direccion/de/carpeta'
# eviDir <- file.path(stepDir, 'evi')
# msaviDir <- file.path(stepDir, 'msavi')
# ndmiDir <- file.path(stepDir, 'ndmi')
# nbrDir <- file.path(stepDir, 'nbr')
# nbr2Dir <- file.path(stepDir, 'nbr2')

#----- Set tempdir to create directory 
#srdir <- dirout <- file.path(rasterOptions()$tmpdir, 'bfmspatial')
#dir.create(dirout, showWarning=FALSE) 


# outDir is where outputs are stored
outDir <- file.path(inDir, 'out')
#outDir <- 'direccion/de/carpeta'
# create an extent object in order to clip processing extent to 60 x 60 km window (UTM 16 N coords).
# It is extremely important for all the images to have the same extent, therefore we have to specify via 
# an extent object with four coordinates in UTM (xmin, xmax, ymin, ymax).

e <- extent(201753.52, 262685.59, 2170910.51, 2231840.28) # extent in Kaxil Kiuic
#e <- extent(159727,221807.4,2269844,2332127) # extent in Palmar 
#e <- extent(331533.8, 391974.4, 2099930, 2160368) #extent in FCP

# processLandsatBatch is variable due to the change in USGS ESPA file naming convention. If using developers version of bfastSpatial use the following to apply the cloud mask: keep = c(322, 386) applies to Landsat 8 data. Change to: keep = c(66, 130) for Landsat 5-7 data
# for landsat 7 keep = c(66, 130))
# for landsat 8 keep = c(322, 386))
# unzip Landsat files, apply cloud mask, and calculate VI if not available, 
# loop is created to read files when they have already created, if not function processLandsatBatch 
# will calculate Vi and apply cloud mask.

# processLandsatBatch(x = landsatDir8, outdir = ndviDir,
#                     delete = TRUE, e=e, overwrite = FALSE, mask = 'pixel_qa', vi = 'ndvi',
#                     keep = c(322, 386))
# 

# 
# if (!file.exists(file.path(inDir, 'ndvi_stack.grd'))) {
#                     processLandsatBatch(x = landsatDir8, outdir = ndviDir,
#                     delete = TRUE, e=e, overwrite = FALSE, mask = 'pixel_qa', vi = 'ndvi',
#                     keep = c(322, 386))
#                 # make temporal ndvi stack
#                     ndviStack <- timeStack(x = ndviDir, pattern = glob2rx('*.grd'),
#                     filename = file.path(inDir, 'ndvi_stack.grd'),
#                     datatype = 'INT2S')
#                     } else {
#                     ndviStack <- brick(file.path(inDir, 'ndvi_stack.grd'))
#                     }
# 
#     
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


# filename = '/home/stevie/Change_detection/FCP/out/ndvi_stack.grd'
# filename = '/home/stevie/Change_detection/PALMAR/out/ndvi_stack.grd'

filename = '/home/stevie/Change_detection/KIUIC/out/ndvi_stack.grd'


ndviStack <- brick(filename)
# when running only ndviStack for 1st time use:
ndviStack <- timeStack(x = ndviDir, pattern = glob2rx('*.grd'), filename = '/home/stevie/Change_detection/KIUIC/out/ndvi_stack.grd', datatype = 'INT2S')



# set ndviStack to x to prepare to run through bfmSpatial
x <- ndviStack

# Section for checking ndviStack has been created  correctly and have a visual overview of data
# show scene info from  ndviStack layers. Check for duplicate dates and incorrect names, which are a problem for bfmSpatial.
# show the layer names
names(x)
  s <- getSceneinfo(names(x))
  s
  
  
  ### ADD THIS FUNCTION
  # filter the scene to show pixels above a certain threshold ( NDVI >.7 in this case)
  # define a function that takes a vector as an argument
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
    # return the value
    return(y)
}

# pass this function summaryBrick
customStat <- summaryBrick(x, fun=checkThresh)
plot(customStat, main = "# of observations where NDVI > 0.7")

# IMPORTANT: create histogram for number of scenes per year
# add a column for years and plot # of scenes per year one must change start:end dates
s$year <- as.numeric(substr(s$date, 1, 4))
hist(s$year, main="Scenes per Year", 
     xlab="year", ylab="# of scenes")

# plot count of valid observations 
obs <- countObs(x)
plot(obs)

# obtain summary of valid observations
obs <- countObs(x, as.perc=TRUE)
summary(obs)
hist(obs, main ="percent obs")
# % NA per pixel 
percNA <- 100 - countObs(x, as.perc=TRUE)
hist(percNA, main = "percent NA per pixel", xlab = "% NA per pixel")
plot(percNA, main="percent NA per pixel")
summary(percNA)

### ADD THIS FUNCTION
### Calculates annual summary statistics of a Raster Brick

annualSummary <- function(x, fun, dates=NULL, years=NULL, sensor=NULL, na.rm=NULL, ...){

    # if sensor is given (!is.null(sensor)), then limit the analysis to a particular sensor
    if(!is.null(sensor)){
        if ("ETM+" %in% sensor) {
            sensor <- unique(c(sensor, "ETM+ SLC-on", "ETM+ SLC-off"))
        }
        if(!.isLandsatSceneID(x)){
            warning("Scene IDs should be supplied as names(x) to subset by sensor. Ignoring...\n")
            scenes <- NULL
        } else {
            # 'allowed' scenes
            scenes <- which(getSceneinfo(names(x))$sensor %in% sensor)
        }
    } else {
        scenes <- NULL
    }
        
    # get dates (if is.null(dates))
    if(is.null(dates)) {
        if(is.null(getZ(x))) {
            if(!.isLandsatSceneID(x)){ # Check if dates can be extracted from layernames
                stop('A date vector must be supplied, either via the date argument, the z dimension of x or comprised in names(x)')
            } else {
                dates <- as.Date(getSceneinfo(names(x))$date)
            }
        } else {
            dates <- getZ(x)
        }
    } else {
        if(length(dates) != nlayers(x)){
            stop("dates should be of same length as nlayers(x)")
        }
    }
    
    # trim dates if a sensor had been supplied
    if(!is.null(scenes)){
        dates <- dates[scenes]
    }
    
    # extract years
    y <- substr(dates, 1, 4)
    
    # vector of years over which to process
    yrs <- sort(unique(y))
    
    # limit to user-defined period
    if(!is.null(years))
        yrs <- yrs[yrs %in% years]
    
    # function to be applied over each pixel in the RasterBrickStack
    pixStat <- function(b){
        if(!is.null(scenes))
            b <- b[scenes]
        ps <- vector("numeric", length(yrs))
        for(i in 1:length(yrs)){
            args <- list(b[which(y == yrs[i])])
            if(is.logical(na.rm))
                args$na.rm <- na.rm
            ps[i] <- do.call(fun, args)
        }
        
        names(ps) <- yrs
        return(ps)
    }
    
    out <- mc.calc(x, fun=pixStat, ...)
    
    return(out)
}


#Calculate Statistics pixel wise per year function median
annualMed <- annualSummary(x, fun=median, na.rm=TRUE)
plot(annualMed, main = "Annual median")

# calculate mean and standard deviation values per year
annualMean <- annualSummary(x, fun=mean, na.rm=TRUE)
plot(annualMean, main = "Annual mean")
hist(annualMean, main = "Annual mean")

annualSD <- annualSummary(x, fun=sd, na.rm=TRUE)
plot(annualSD, main = "Annual SD")
hist(annualSD, main "Annual SD")
# custom function to calculate # of non-NA values per pixel per year (similar to countObs())
ff <- function(x)
length(x[!is.na(x)])
annualObs <- annualSummary(x,fun=ff)


######################################################################################################
# BFM INTERACTIVE
# in order to use bfmPixel in interactive mode it is important to first plot a layer with visible data.
plot(x, 322)

# run bfmPixel() in interactive mode with a monitoring period 
# starting @ the 1st day in 1990
bfm <- bfmPixel(x, start=2016, interactive=TRUE)
# choose the pixel whose time series you want to see by clicking on the map you plotted a moment ago
bfm$cell
# assign same pixel as a target cell
targcell <- bfm$cell
# inspect and plot the $bfm output
bfm <- bfmPixel(x, cell=targcell, start=c(2016, 1))
bfm
plot(bfm$bfm)

# trial pixel run for research parameter adjustment
bfm2 <- bfmPixel(x, cell=targcell, start=1995, 
                 formula=response~trend+harmon, h=0.25, history=c(1991, 4),type="OLS-MOSUM", order=1, plot=TRUE)
				 
				 # only trend
#bfm3 <- bfmPixel(x, cell=targcell, start=c(1990, 1), 
             #    formula=response~trend, plot=TRUE)
				 
# bfmPixel using a 1-year monitoring period
#bfm4 <- bfmPixel(x, cell=targcell, start=c(1990, 1), 
            #     monend=c(1991, 1), plot=TRUE)
				 
# it is possible to filter the analysis by sensor 
# e.g. apply bfmPixel only on ETM+ data 
#bfm5 <- bfmPixel(x, cell=targcell, start=c(2009, 1), 
           #      sensor="ETM+", plot=TRUE)


# get MEDIAN value for timestack, it is also possible to obtain mean and other statistics and plot them
if (!file.exists(fn <- file.path(outDir, 'medianVI.grd')))
  # median values for all layers
  medVI <- summaryBrick(ndviStack, fun=median, na.rm=TRUE, 
                        filename = fn) else {
                          medVI <- brick(fn)
                        }
# produce median plot
plot(medVI/10000)


#######################################################################################################
# BFMSPATIAL
# run bfmSpatial through the entire timeStack (all pixels)
# This is a computationally demanding process, multicore processing is preferred and can be specified (mc.cores= ) a set of parameters need to be set: formula (harmonic or trend or both), history, start of monitoring period (start) level (significance), (h=) it is possible to retrieve the following layers via returnLayers = c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared", "coefficients"). For editing parameters start= , monend, history


 
# Write breakpoint, yearly break month product, and breakpoint magnitude raster layers to GeoTiff files as well as the raster brick to a .grd file.
writeRaster(out[[1]], filename = "Kiuic1_NDVI_breaks.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth1995, filename = "Kiuic1_NDVI_breaksmos95.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth1996, filename = "Kiuic1_NDVI_breaksmos96.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth1997, filename = "Kiuic1_NDVI_breaksmos97.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth1998, filename = "Kiuic1_NDVI_breaksmos98.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth1999, filename = "Kiuic1_NDVI_breaksmos99.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2000, filename = "Kiuic1_NDVI_breaksmos00.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2001, filename = "Kiuic1_NDVI_breaksmos01.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2002, filename = "Kiuic1_NDVI_breaksmos02.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2003, filename = "Kiuic1_NDVI_breaksmos03.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2004, filename = "Kiuic1_NDVI_breaksmos04.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2005, filename = "Kiuic1_NDVI_breaksmos05.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2006, filename = "Kiuic1_NDVI_breaksmos06.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2007, filename = "Kiuic1_NDVI_breaksmos07.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2008, filename = "Kiuic1_NDVI_breaksmos08.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2009, filename = "Kiuic1_NDVI_breaksmos09.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2010, filename = "Kiuic1_NDVI_breaksmos10.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2011, filename = "Kiuic1_NDVI_breaksmos11.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2012, filename = "Kiuic1_NDVI_breaksmos12.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2013, filename = "Kiuic1_NDVI_breaksmos13.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2014, filename = "Kiuic1_NDVI_breaksmos14.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2015, filename = "Kiuic1_NDVI_breaksmos15.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2016, filename = "Kiuic1_NDVI_breaksmos16.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2017, filename = "Kiuic1_NDVI_breaksmos17.tif", format = "GTiff", overwrite = TRUE)
writeRaster(months$changeMonth2018, filename = "Kiuic1_NDVI_breaksmos18.tif", format = "GTiff", overwrite = TRUE)

# Test breakpoints
plot(ndviStack[[322]], col = grey.colors(255), legend = F)
plot(bfm[[1]], add=TRUE)
# Test months product
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)
# Test magnitudes
plot(magn_bkp, main="Magnitude of a breakpoint")
plot(magn, main="Magnitude: all pixels") 
 
rsq_map <- raster(bfm, 5) 
 
# Extract R2 from map

writeRaster( 
  rsq_map, 
  filename="rsq_95-4y_Kiuic.tif", 
  format="GTiff", 
  overwrite=TRUE)  
  
#Apply ndviStack function to all images

x <- ndviStack	

#####  Predicting one year with all previous history #######################################
############################################################################################

  
  ######################################### 2020  ###########################################
      
        
        if (!file.exists(fn <- file.path(outDir, 'bfm_20-all_harmon_1.grd'))) {
          # make sure the filename is correct and the folder exists
          bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                            start = 2020, history = c("all"), h=0.25, monend=NULL, level=0.05,
                            formula = response ~ harmon, type = "OLS-MOSUM",
                            order = 1, mc.cores = 1, returnLayers = c("breakpoint", "magnitude", "error",
                                                                      "history", "r.squared", "adj.r.squared", "coefficients"),
                            filename = file.path(outDir, 'bfm_20-all_harmon_1.grd')
          )
        } else {
          bfm <- brick(fn)
        }
        
        
  
        
  
  
  
######################################### 2019  ###########################################
            
            
            
            if (!file.exists(fn <- file.path(outDir, 'bfm_19-all_harmon_1.grd'))) {
              # make sure the filename is correct and the folder exists
              bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                                start = 2019, history = c("all"), h=0.25, monend=2020, level=0.05,
                                formula = response ~ harmon, type = "OLS-MOSUM",
                                order = 1, mc.cores = 1, returnLayers = c("breakpoint", "magnitude", "error",
                                                                          "history", "r.squared", "adj.r.squared", "coefficients"),
                                filename = file.path(outDir, 'bfm_19-all_harmon_1.grd')
              )
            } else {
              bfm <- brick(fn)
            }
            
            

######################################### 2018  ###########################################
        
        
        if (!file.exists(fn <- file.path(outDir, 'bfm_18-all_harmon_1.grd'))) {
          # make sure the filename is correct and the folder exists
          bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                            start = 2018, history = c("all"), h=0.25, monend=2019, level=0.05,
                            formula = response ~ harmon, type = "OLS-MOSUM",
                            order = 1, mc.cores = 1, returnLayers = c("breakpoint", "magnitude", "error",
                                                                      "history", "r.squared", "adj.r.squared", "coefficients"),
                            filename = file.path(outDir, 'bfm_18-all_harmon_1.grd')
          )
        } else {
          bfm <- brick(fn)
        }
  
######################################### 2017 ###########################################
  
  
  if (!file.exists(fn <- file.path(outDir, 'bfm_17-all_harmon_1.grd'))) {
    # make sure the filename is correct and the folder exists
    bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                      start = 2017, history = c("all"), h=0.25, monend=2018, level=0.05,
                      formula = response ~ harmon, type = "OLS-MOSUM",
                      order = 1, mc.cores = 1, returnLayers = c("breakpoint", "magnitude", "error",
                                                                "history", "r.squared", "adj.r.squared", "coefficients"),
                      filename = file.path(outDir, 'bfm_17-all_harmon_1.grd')
    )
  } else {
    bfm <- brick(fn)
  }

#########################################   2016   ###########################################


if (!file.exists(fn <- file.path(outDir, 'bfm_16-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2016, history = c("all"), h=0.25, monend=2017, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_16-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

########################################   2015   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_15-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2015, history = c("all"), h=0.25, monend=2016, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_15-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}
#########################################   2014   ###########################################


if (!file.exists(fn <- file.path(outDir, 'bfm_14-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2014, history = c("all"), h=0.25, monend=2015, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_14-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2013   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_13-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2013, history = c("all"), h=0.25, monend=2014, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_13-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2012   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_12-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2012, history = c("all"), h=0.25, monend=2013, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_12-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2011   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_11-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2011, history = c("all"), h=0.25, monend=2012, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_11-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2010   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_10-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2010, history = c("all"), h=0.25, monend=2011, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_10-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2009   ###########################################


if (!file.exists(fn <- file.path(outDir, 'bfm_09-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2009, history = c("all"), h=0.25, monend=2010, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_09-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2008   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_08-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2008, history = c("all"), h=0.25, monend=2009, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_08-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2007   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_07-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2007, history = c("all"), h=0.25, monend=2008, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_07-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2006   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_06-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2006, history = c("all"), h=0.25, monend=2007, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_06-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2005   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_05-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2005, history = c("all"), h=0.25, monend=2006, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_05-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2004   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_04-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2004, history = c("all"), h=0.25, monend=2005, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_04-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2003   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_03-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2003, history = c("all"), h=0.25, monend=2004, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_03-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2002   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_02-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2002, history = c("all"), h=0.25, monend=2003, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_02-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2001   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_01-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2001, history = c("all"), h=0.25, monend=2002, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_01-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   2000   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_00-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 2000, history = c("all"), h=0.25, monend=2001, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_00-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   1999   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_99-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 1999, history = c("all"), h=0.25, monend=2000, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_99-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   1998   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_98-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 1998, history = c("all"), h=0.25, monend=1999, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_98-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   1997   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_97-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 1997, history = c("all"), h=0.25, monend=1998, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_97-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   1996   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_96-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 1996, history = c("all"), h=0.25, monend=1997, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_96-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}

#########################################   1995   ###########################################

if (!file.exists(fn <- file.path(outDir, 'bfm_95-all_harmon_1.grd'))) {
  # make sure the filename is correct and the folder exists
  bfm <- bfmSpatial(x = ndviStack, pptype = 'irregular',
                    start = 1995, history = c("all"), h=0.25, monend=1996, level=0.05,
                    formula = response ~ harmon, type = "OLS-MOSUM",
                    order = 1, mc.cores = 4, returnLayers = c("breakpoint", "magnitude", "error",
                                                              "history", "r.squared", "adj.r.squared", "coefficients"),
                    filename = file.path(outDir, 'bfm_95-all_harmon_1.grd')
  )
} else {
  bfm <- brick(fn)
}



# Extraction of year segment output
#------------------------ 2020 ------------------------------------------------------

# filename: bfm_20-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_20-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

plot(change, main="Change: Breakpoints 2020")

writeRaster( 
  change, 
  filename="breaks_20-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)

plot(magn, main="Magnitude")
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_20-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_20-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  


#------------------------ 2019 ------------------------------------------------------

# filename: bfm_19-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_19-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_19-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_19-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_19-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  







#------------------------ 2018 ------------------------------------------------------

# filename: bfm_18-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_18-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_18-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_18-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_18-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  





#------------------------ 2017 ------------------------------------------------------

# filename: bfm_17-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_17-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_17-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_17-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_17-all_FCP.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2016 ------------------------------------------------------

# filename: bfm_16-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_16-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_16-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_16-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_16-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2015 ------------------------------------------------------

# filename: bfm_15-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_15-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_15-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_15-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_15-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2014 ------------------------------------------------------

# filename: bfm_14-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_14-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_14-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_14-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_14-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  


#------------------------ 2013 ------------------------------------------------------

# filename: bfm_13-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_13-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_13-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_13-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_13-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2012 ------------------------------------------------------

# filename: bfm_12-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_12-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_12-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_12-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_12-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2011 ------------------------------------------------------

# filename: bfm_11-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_11-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_11-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_11-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_11-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  


#------------------------ 2010 ------------------------------------------------------

# filename: bfm_10-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_10-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_10-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_10-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_10-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2009 ------------------------------------------------------

# filename: bfm_09-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_09-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_09-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_09-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_09-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2008 ------------------------------------------------------

# filename: bfm_08-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_08-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_08-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_08-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_08-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)

#------------------------ 2007 ------------------------------------------------------

# filename: bfm_07-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_07-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_07-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_07-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_07-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  


#------------------------ 2006 ------------------------------------------------------

# filename: bfm_06-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_06-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_06-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_06-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_06-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2005 ------------------------------------------------------

# filename: bfm_05-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_05-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_05-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_05-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_05-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  




#------------------------ 2004 ------------------------------------------------------

# filename: bfm_04-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_04-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_04-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_04-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_04-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2003 ------------------------------------------------------

# filename: bfm_03-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_03-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_03-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_03-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_03-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2002 ------------------------------------------------------

# filename: bfm_02-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_02-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_02-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_02-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_02-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  


#------------------------ 2001 ------------------------------------------------------

# filename: bfm_01-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_01-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_01-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_01-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_01-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 2000 ------------------------------------------------------

# filename: bfm_00-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_00-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_00-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_00-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_00-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 1999 ------------------------------------------------------

# filename: bfm_99-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_99-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_99-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_99-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_99-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 1998 ------------------------------------------------------

# filename: bfm_98-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_98-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_98-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_98-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_98-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 1997 ------------------------------------------------------

# filename: bfm_97-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_97-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_97-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_97-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_97-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

#------------------------ 1996 ------------------------------------------------------

# filename: bfm_96-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_96-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_96-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_96-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_96-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)

#------------------------ 1995 ------------------------------------------------------

# filename: bfm_95-all_harmon_1.grd

fn <- file.path(outDir, 'bfm_95-all_harmon_1.grd')

bfm <- brick(fn)

# the first layer corresponds to breakpoints per year, which will correspond to change.
# write a raster to extract change layer

change <- raster(bfm, 1)

writeRaster( 
  change, 
  filename="breaks_95-all_harmon_1.tif", 
  format="GTiff", 
  overwrite=TRUE)  


# make months object
months <- changeMonth(change)
# set up labels and color map for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)

# produce plot with changes per month
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# second layer corresponds to magnitude 
magn <- raster(bfm, 2)
# store magnitude output in separate layer
writeRaster( 
  magn, 
  filename="magn_95-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)  

# make a version showing only breakpoint pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")


# store magnitude bkpts output in separate layer
writeRaster( 
  magn_bkp, 
  filename="magn_bkp_95-all_Palmar.tif", 
  format="GTiff", 
  overwrite=TRUE)    