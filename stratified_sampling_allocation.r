### Accuracy sampling
## stephanie.p.george@gmail.com

library(bfastSpatial)
library(rgeos)


## estimating proportions

hi_buffer <- readOGR("/home/stevie/Change_detection/shapes/strata/high_impact_buffer.shp")
med_buffer <- readOGR("/home/stevie/Change_detection/shapes/strata/medium_impact_buffer.shp")
lo_buffer <- readOGR("/home/stevie/Change_detection/shapes/strata/low_impact_buffer.shp")

plot(hi_buffer, col = "orange")
plot(med_buffer, add =T, col = "yellow")
plot(lo_buffer, add = T, col = "green")

total_poly <- readOGR("/home/stevie/Change_detection/KIUIC/shapes/Kiuic_square_poly_01_utm16.shp")


hi_buffer_site <- crop(hi_buffer, total_poly)
med_buffer_site <- crop(med_buffer, total_poly)
lo_buffer_site <- crop(lo_buffer, total_poly)



plot(hi_buffer_site, col = "orange")
plot(med_buffer_site,add=T, col = "yellow")
plot(lo_buffer_site, add = T, col = "green")


#change_year <- raster("/home/stevie/Change_detection/KIUIC/yearmap15.tif")
change_year <- raster("/home/stevie/Change_detection/KIUIC/magn_areasieve_all_15.tif")



hi_change <- mask(change_year, hi_buffer_site)
plot(hi_change, add = T, col="red")
med_change <- mask(change_year, med_buffer_site)
plot(med_change, add = T, col="red")
lo_change <- mask(change_year, lo_buffer_site)
plot(lo_change, add = T, col="red")

hi_change_poly <- rasterToPolygons(hi_change)
med_change_poly <- rasterToPolygons(med_change)
lo_change_poly <- rasterToPolygons(lo_change)

area_change/total_area



## area calculation

total_area <- gArea(total_poly)
total_area_ha <- total_area/10000
total_area_ha

hi_change_area <- gArea(hi_change_poly)
med_change_area <- gArea(med_change_poly)
lo_change_area <- gArea(lo_change_poly)

total_change_area <- 
hi_nochange <- (gArea(hi_buffer)) - hi_change_area
med_nochange <- (gArea(med_buffer)) - med_change_area
lo_nochange <- (gArea(lo_buffer)) - lo_change_area


(hi_change_area + med_change_area + lo_change_area)*100/total_area


### weights calculation
hi_w_c <- hi_change_area/total_area
hi_w_nc <- hi_nochange/total_area

med_w_c <- med_change_area/total_area
med_w_nc <- med_nochange/total_area

lo_w_c <- lo_change_area/total_area
lo_w_nc <- lo_nochange/total_area
  
  
  si_c  = sqrt(.70*(1-.70))
  si_nc = sqrt(.95*(1-.95))
  so = .02 
  
  
  # get total N according to Cochran's 1977 formula
  N = (((hi_w_c * si_c) + 
          (hi_w_nc * si_nc) + 
          (med_w_c * si_c) + 
          (med_w_nc * si_nc) + 
          (lo_w_c * si_c) + 
          (lo_w_nc * si_nc))/so)**2
  
  N = round(N)
  N
  n_equal = round(N/6)
  n_equal


#--------------------------------------

# Neyman's allocation 
n_hi_c =  N * ((hi_w_c * si_c)/((hi_w_c * si_c) + 
                              (hi_w_nc * si_nc) + 
                              (med_w_c * si_c) + 
                              (med_w_nc * si_nc) + 
                              (lo_w_c * si_c) + 
                              (lo_w_nc * si_nc)))
n_hi_c
n_hi_nc =  N * ((hi_w_nc * si_nc)/((hi_w_c * si_c) + 
                                  (hi_w_nc * si_nc) + 
                                  (med_w_c * si_c) + 
                                  (med_w_nc * si_nc) + 
                                  (lo_w_c * si_c) + 
                                  (lo_w_nc * si_nc)))
n_hi_nc
n_med_nc =  N * ((med_w_nc * si_nc)/((hi_w_c * si_c) + 
                                     (hi_w_nc * si_nc) + 
                                     (med_w_c * si_c) + 
                                     (med_w_nc * si_nc) + 
                                     (lo_w_c * si_c) + 
                                     (lo_w_nc * si_nc)))
n_med_nc
n_med_c =  N * ((med_w_c * si_c)/((hi_w_c * si_c) + 
                                       (hi_w_nc * si_nc) + 
                                       (med_w_c * si_c) + 
                                       (med_w_nc * si_nc) + 
                                       (lo_w_c * si_c) + 
                                       (lo_w_nc * si_nc)))
n_med_c
n_lo_c =  N * ((lo_w_c * si_c)/((hi_w_c * si_c) + 
                                    (hi_w_nc * si_nc) + 
                                    (med_w_c * si_c) + 
                                    (med_w_nc * si_nc) + 
                                    (lo_w_c * si_c) + 
                                    (lo_w_nc * si_nc)))
n_lo_c
n_lo_nc =  N * ((lo_w_nc * si_nc)/((hi_w_c * si_c) + 
                                  (hi_w_nc * si_nc) + 
                                  (med_w_c * si_c) + 
                                  (med_w_nc * si_nc) + 
                                  (lo_w_c * si_c) + 
                                  (lo_w_nc * si_nc)))
n_lo_nc  

n_total_change <- n_lo_c + n_hi_c + n_med_c
n_total_n_change <- n_lo_nc + n_hi_nc + n_med_nc

    


#get sizes of all cells in raster [km2]
cell_size<-area(change_year, na.rm=TRUE, weights=FALSE)
#delete NAs from vector of all raster cells
##NAs lie outside of the rastered region, can thus be omitted
cell_size<-cell_size[!is.na(cell_size)]
#compute area [km2] of all cells in geo_raster
raster_area<-length(cell_size)*median(cell_size)
#print area of Georgia according to raster object
print(paste("Area of Kiuic (raster):",round(raster_area, digits=1),"km2"))



change_poly <- rasterToPolygons(change_year)
total_change_area <- gArea(change_poly)
total_area <- gArea(total_poly)
wic <- total_change_area/total_area
winc <- 1 - wic  
  