#### Validation of time series using timeSyncR 

# library(devtools)
# install_github("bendv/timeSyncR")
library(timeSyncR)

r <- brick("test/MDD_red.grd")
g <- brick("test/MDD_green.grd")
b <- brick("test/MDD_blue.grd")

r <- brick(b4_stack)
g <- brick(b3_stack)
b <- brick(b2_stack)
myRaster <- writeRaster(b,"myStack.grd", format="raster")
raster_brick <- brick(myRaster)
nlayers(r)
plot(g, 1:9)

# layer names correspond to Landsat ID's
names(raster_brick)
# see ?getSceneinfo to get more info out of these
s <- getSceneinfo(names(b))
print(b)

# each brick has a z-dimension containing date info
getZ(b)
# also found in the date column of the getSceneinfo() data.frame
s$date
all(getZ(b) == s$date)

xy <- c(251908.5, 2211539)
# save original (default) plotting parameters to workspace
op <- par()

tsChipsRGB(xr = r, xg = g, xb = b, loc = xy, start = "2000-01-01")
# reset plotting parameters (plotRGB() changes parameters internally)
par(op)

r$names <- as.numeric(rownames(s))

library(stringr)

library(stringr)
g$names <- str_sub(g$names, end=-2)
head(r)
pix <- pixelToPolygon(x = r, cell = xy)
plot(pix)
tsChipsRGB(xr = r, xg = g, xb = b, loc = pix, start = "2000-01-01") 
head(g)
