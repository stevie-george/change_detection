## Plotting overlap Deforestation/nodeforestation area
# Stephanie P. George 
## 25/01/2021

library(rgdal)
library(raster)
library(rgeos)

#magn_bkp_15-all_Kiuic.tif

setwd('/home/stevie/Change_detection/KIUIC/tifs/')

years <- c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20")
rasternames <- paste("magn_bkp_", years,"-all_Kiuic.tif", sep="")
rasternames

r00 <- raster(rasternames[1])/10000
r00[r00 > -0.061] <- NA
raster00 <- gArea(rasterToPolygons(r00))/10000 

r01 <- raster(rasternames[2])/10000
r01[r01 > -0.061] <- NA
raster01 <- gArea(rasterToPolygons(r01))/10000 

r02 <- raster(rasternames[3])/10000
r02[r02 > -0.061] <- NA
raster02 <- gArea(rasterToPolygons(r02))/10000 

r03 <- raster(rasternames[4])/10000
r03[r03 > -0.061] <- NA
raster03 <- gArea(rasterToPolygons(r03))/10000 

r04 <- raster(rasternames[5])/10000
r04[r04 > -0.061] <- NA
raster04 <- gArea(rasterToPolygons(r04))/10000 

r05 <- raster(rasternames[6])/10000
r05[r05 > -0.061] <- NA
raster05 <- gArea(rasterToPolygons(r05))/10000 


r06 <- raster(rasternames[7])/10000
r06[r06 > -0.061] <- NA
raster06 <- gArea(rasterToPolygons(r06))/10000 

r07 <- raster(rasternames[8])/10000
r07[r07 > -0.061] <- NA
raster07 <- gArea(rasterToPolygons(r07))/10000 


r08 <- raster(rasternames[9])/10000
r08[r08 > -0.061] <- NA
raster08 <- gArea(rasterToPolygons(r08))/10000 


r09 <- raster(rasternames[10])/10000
r09[r09 > -0.061] <- NA
raster09 <- gArea(rasterToPolygons(r09))/10000 


r10 <- raster(rasternames[11])/10000
r10[r10 > -0.061] <- NA
raster10 <- gArea(rasterToPolygons(r10))/10000 


r11 <- raster(rasternames[12])/10000
r11[r11 > -0.061] <- NA
raster11 <- gArea(rasterToPolygons(r11))/10000 


r12 <- raster(rasternames[13])/10000
r12[r12 > -0.061] <- NA
raster12 <- gArea(rasterToPolygons(r12))/10000 


r13 <- raster(rasternames[14])/10000
r13[r13 > -0.061] <- NA
raster13 <- gArea(rasterToPolygons(r13))/10000 


r14 <- raster(rasternames[15])/10000
r14[r14 > -0.061] <- NA
raster14 <- gArea(rasterToPolygons(r14))/10000 

r15 <- raster(rasternames[16])/10000
r15[r15 > -0.061] <- NA
raster15 <- gArea(rasterToPolygons(r15))/10000 


r16 <- raster(rasternames[17])/10000
r16[r16 > -0.061] <- NA
raster16 <- gArea(rasterToPolygons(r16))/10000 

r17 <- raster(rasternames[18])/10000
r17[r17 > -0.061] <- NA
raster17 <- gArea(rasterToPolygons(r17))/10000 

r18 <- raster(rasternames[19])/10000
r18[r18 > -0.061] <- NA
raster18 <- gArea(rasterToPolygons(r18))/10000 


r19 <- raster(rasternames[20])/10000
r19[r19 > -0.061] <- NA
raster19 <- gArea(rasterToPolygons(r19))/10000 


r20 <- raster(rasternames[21])/10000
r20[r20 > -0.061] <- NA
raster20 <- gArea(rasterToPolygons(r20))/10000 


rasters <- rbind(raster00, raster01, raster02, raster03, raster04, raster05, raster06,
                 raster07, raster08, raster09, raster10, raster11, raster12, raster13, 
                 raster14, raster15, raster16, raster17, raster18, raster19, raster20)
write.csv(rasters, "total_change_kiuic")
year <- 2000:2020

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(year, rasters, type="l", col = "blue",lwd =3, ylab = "Forest cover loss (ha)", ylim =c(0,50000), xlab="Year", cex.lab=1.5, cex.axis = 1.2) # kiuic

##### PALMAR
setwd('/home/stevie/Change_detection/PALMAR/tifs/')

years <- c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20")
#rasternames <- paste("magn", years,"areasieve.tif", sep="")
rasternames <- paste("magn_bkp_", years,"-all_Palmar.tif", sep="")
rasternames

# raster00 <- gArea(rasterToPolygons(raster(rasternames[1])))/10000 
# raster01 <- gArea(rasterToPolygons(raster(rasternames[2])))/10000 
# raster02 <- gArea(rasterToPolygons(raster(rasternames[3])))/10000 
# raster03 <- gArea(rasterToPolygons(raster(rasternames[4])))/10000 
# raster04 <- gArea(rasterToPolygons(raster(rasternames[5])))/10000 
# raster05 <- gArea(rasterToPolygons(raster(rasternames[6])))/10000 
# raster06 <- gArea(rasterToPolygons(raster(rasternames[7])))/10000 
# raster07 <- gArea(rasterToPolygons(raster(rasternames[8])))/10000 
# raster08 <- gArea(rasterToPolygons(raster(rasternames[8])))/10000 
# raster09 <- gArea(rasterToPolygons(raster(rasternames[9])))/10000 
# raster10 <- gArea(rasterToPolygons(raster(rasternames[10])))/10000 
# raster11 <- gArea(rasterToPolygons(raster(rasternames[11])))/10000 
# raster12 <- gArea(rasterToPolygons(raster(rasternames[12])))/10000 
# raster13 <- gArea(rasterToPolygons(raster(rasternames[13])))/10000 
# raster14 <- gArea(rasterToPolygons(raster(rasternames[14])))/10000 
# raster15 <- gArea(rasterToPolygons(raster(rasternames[15])))/10000 
# raster16 <- gArea(rasterToPolygons(raster(rasternames[16])))/10000 
# raster17 <- gArea(rasterToPolygons(raster(rasternames[17])))/10000 
# raster18 <- gArea(rasterToPolygons(raster(rasternames[18])))/10000 
# raster19 <- gArea(rasterToPolygons(raster(rasternames[19])))/10000 
# raster20 <- gArea(rasterToPolygons(raster(rasternames[20])))/10000 

r00 <- raster(rasternames[1])/10000
r00[r00 > -0.061] <- NA
raster00 <- gArea(rasterToPolygons(r00))/10000 

r01 <- raster(rasternames[2])/10000
r01[r01 > -0.061] <- NA
raster01 <- gArea(rasterToPolygons(r01))/10000 

r02 <- raster(rasternames[3])/10000
r02[r02 > -0.061] <- NA
raster02 <- gArea(rasterToPolygons(r02))/10000 

r03 <- raster(rasternames[4])/10000
r03[r03 > -0.061] <- NA
raster03 <- gArea(rasterToPolygons(r03))/10000 

r04 <- raster(rasternames[5])/10000
r04[r04 > -0.061] <- NA
raster04 <- gArea(rasterToPolygons(r04))/10000 

r05 <- raster(rasternames[6])/10000
r05[r05 > -0.061] <- NA
raster05 <- gArea(rasterToPolygons(r05))/10000 


r06 <- raster(rasternames[7])/10000
r06[r06 > -0.061] <- NA
raster06 <- gArea(rasterToPolygons(r06))/10000 

r07 <- raster(rasternames[8])/10000
r07[r07 > -0.061] <- NA
raster07 <- gArea(rasterToPolygons(r07))/10000 


r08 <- raster(rasternames[9])/10000
r08[r08 > -0.061] <- NA
raster08 <- gArea(rasterToPolygons(r08))/10000 


r09 <- raster(rasternames[10])/10000
r09[r09 > -0.061] <- NA
raster09 <- gArea(rasterToPolygons(r09))/10000 


r10 <- raster(rasternames[11])/10000
r10[r10 > -0.061] <- NA
raster10 <- gArea(rasterToPolygons(r10))/10000 


r11 <- raster(rasternames[12])/10000
r11[r11 > -0.061] <- NA
raster11 <- gArea(rasterToPolygons(r11))/10000 


r12 <- raster(rasternames[13])/10000
r12[r12 > -0.061] <- NA
raster12 <- gArea(rasterToPolygons(r12))/10000 


r13 <- raster(rasternames[14])/10000
r13[r13 > -0.061] <- NA
raster13 <- gArea(rasterToPolygons(r13))/10000 


r14 <- raster(rasternames[15])/10000
r14[r14 > -0.061] <- NA
raster14 <- gArea(rasterToPolygons(r14))/10000 

r15 <- raster(rasternames[16])/10000
r15[r15 > -0.061] <- NA
raster15 <- gArea(rasterToPolygons(r15))/10000 


r16 <- raster(rasternames[17])/10000
r16[r16 > -0.061] <- NA
raster16 <- gArea(rasterToPolygons(r16))/10000 

r17 <- raster(rasternames[18])/10000
r17[r17 > -0.061] <- NA
raster17 <- gArea(rasterToPolygons(r17))/10000 

r18 <- raster(rasternames[19])/10000
r18[r18 > -0.061] <- NA
raster18 <- gArea(rasterToPolygons(r18))/10000 


r19 <- raster(rasternames[20])/10000
r19[r19 > -0.061] <- NA
raster19 <- gArea(rasterToPolygons(r19))/10000 


r20 <- raster(rasternames[21])/10000
r20[r20 > -0.061] <- NA
raster20 <- gArea(rasterToPolygons(r20))/10000 




rasters <- rbind(raster00, raster01, raster02, raster03, raster04, raster05, raster06,
                 raster07, raster08, raster09, raster10, raster11, raster12, raster13, 
                 raster14, raster15, raster16, raster17, raster18, raster19, raster20)
write.csv(rasters, "total_change_palmar")
year <- 2000:2020

lines(year, rasters, type="l", col = "red", lwd=3) # palmar


setwd('/home/stevie/Change_detection/FCP/tifs/')

years <- c("00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20")
# rasternames <- paste("magn", years,"areasieve.tif", sep="")
rasternames <- paste("magn_bkp_", years,"-all_FCP.tif", sep="")
rasternames

# raster00 <- gArea(rasterToPolygons(raster(rasternames[1])))/10000 
# raster01 <- gArea(rasterToPolygons(raster(rasternames[2])))/10000 
# raster02 <- gArea(rasterToPolygons(raster(rasternames[3])))/10000 
# raster03 <- gArea(rasterToPolygons(raster(rasternames[4])))/10000 
# raster04 <- gArea(rasterToPolygons(raster(rasternames[5])))/10000 
# raster05 <- gArea(rasterToPolygons(raster(rasternames[6])))/10000 
# raster06 <- gArea(rasterToPolygons(raster(rasternames[7])))/10000 
# raster07 <- gArea(rasterToPolygons(raster(rasternames[8])))/10000 
# raster08 <- gArea(rasterToPolygons(raster(rasternames[8])))/10000 
# raster09 <- gArea(rasterToPolygons(raster(rasternames[9])))/10000 
# raster10 <- gArea(rasterToPolygons(raster(rasternames[10])))/10000 
# raster11 <- gArea(rasterToPolygons(raster(rasternames[11])))/10000 
# raster12 <- gArea(rasterToPolygons(raster(rasternames[12])))/10000 
# raster13 <- gArea(rasterToPolygons(raster(rasternames[13])))/10000 
# raster14 <- gArea(rasterToPolygons(raster(rasternames[14])))/10000 
# raster15 <- gArea(rasterToPolygons(raster(rasternames[15])))/10000 
# raster16 <- gArea(rasterToPolygons(raster(rasternames[16])))/10000 
# raster17 <- gArea(rasterToPolygons(raster(rasternames[17])))/10000 
# raster18 <- gArea(rasterToPolygons(raster(rasternames[18])))/10000 
# raster19 <- gArea(rasterToPolygons(raster(rasternames[19])))/10000 
# raster20 <- gArea(rasterToPolygons(raster(rasternames[20])))/10000 

r00 <- raster(rasternames[1])/10000
r00[r00 > -0.061] <- NA
raster00 <- gArea(rasterToPolygons(r00))/10000 

r01 <- raster(rasternames[2])/10000
r01[r01 > -0.061] <- NA
raster01 <- gArea(rasterToPolygons(r01))/10000 

r02 <- raster(rasternames[3])/10000
r02[r02 > -0.061] <- NA
raster02 <- gArea(rasterToPolygons(r02))/10000 

r03 <- raster(rasternames[4])/10000
r03[r03 > -0.061] <- NA
raster03 <- gArea(rasterToPolygons(r03))/10000 

r04 <- raster(rasternames[5])/10000
r04[r04 > -0.061] <- NA
raster04 <- gArea(rasterToPolygons(r04))/10000 

r05 <- raster(rasternames[6])/10000
r05[r05 > -0.061] <- NA
raster05 <- gArea(rasterToPolygons(r05))/10000 


r06 <- raster(rasternames[7])/10000
r06[r06 > -0.061] <- NA
raster06 <- gArea(rasterToPolygons(r06))/10000 

r07 <- raster(rasternames[8])/10000
r07[r07 > -0.061] <- NA
raster07 <- gArea(rasterToPolygons(r07))/10000 


r08 <- raster(rasternames[9])/10000
r08[r08 > -0.061] <- NA
raster08 <- gArea(rasterToPolygons(r08))/10000 


r09 <- raster(rasternames[10])/10000
r09[r09 > -0.061] <- NA
raster09 <- gArea(rasterToPolygons(r09))/10000 


r10 <- raster(rasternames[11])/10000
r10[r10 > -0.061] <- NA
raster10 <- gArea(rasterToPolygons(r10))/10000 


r11 <- raster(rasternames[12])/10000
r11[r11 > -0.061] <- NA
raster11 <- gArea(rasterToPolygons(r11))/10000 


r12 <- raster(rasternames[13])/10000
r12[r12 > -0.061] <- NA
raster12 <- gArea(rasterToPolygons(r12))/10000 


r13 <- raster(rasternames[14])/10000
r13[r13 > -0.061] <- NA
raster13 <- gArea(rasterToPolygons(r13))/10000 


r14 <- raster(rasternames[15])/10000
r14[r14 > -0.061] <- NA
raster14 <- gArea(rasterToPolygons(r14))/10000 

r15 <- raster(rasternames[16])/10000
r15[r15 > -0.061] <- NA
raster15 <- gArea(rasterToPolygons(r15))/10000 


r16 <- raster(rasternames[17])/10000
r16[r16 > -0.061] <- NA
raster16 <- gArea(rasterToPolygons(r16))/10000 

r17 <- raster(rasternames[18])/10000
r17[r17 > -0.061] <- NA
raster17 <- gArea(rasterToPolygons(r17))/10000 

r18 <- raster(rasternames[19])/10000
r18[r18 > -0.061] <- NA
raster18 <- gArea(rasterToPolygons(r18))/10000 


r19 <- raster(rasternames[20])/10000
r19[r19 > -0.061] <- NA
raster19 <- gArea(rasterToPolygons(r19))/10000 


r20 <- raster(rasternames[21])/10000
r20[r20 > -0.061] <- NA
raster20 <- gArea(rasterToPolygons(r20))/10000 

rasters <- rbind(raster00, raster01, raster02, raster03, raster04, raster05, raster06,
                 raster07, raster08, raster09, raster10, raster11, raster12, raster13, 
                 raster14, raster15, raster16, raster17, raster18, raster19, raster20)

write.csv(rasters, "total_change_fcp")
year <- 2000:2020

lines(year, rasters, type="l", col = "green", lwd = 3) # fcp
legend("topleft", legend=c("Semi-deciduous", "Semi-evergreen", "Deciduous"), col=c("blue", "green", "red"), lty=1, lwd = 3)

# Draw confidence intervals



