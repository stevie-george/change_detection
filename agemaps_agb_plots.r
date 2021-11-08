# plotting
library(raster)
par(mfrow=c(1,3))

FCP <- raster('/home/stevie/Change_detection/FCP/tifs/sieves/FCP_AGB.tif')
Kiuic <- raster('/home/stevie/Change_detection/KIUIC//tifs/sieves/Kiuic_AGB.tif')
Palmar <- raster('/home/stevie/Change_detection/PALMAR/tifs/palmar_agb.tif')


plot(FCP, main = "FCP")
plot(Kiuic, main = "KK")
plot(Palmar, main = "EP")


age_FCP <- raster('/home/stevie/Change_detection/FCP/agemap/FCP_agemap_2020.tif')
age_Kiuic <- raster('/home/stevie/Change_detection/KIUIC/agemap/Kiuic_agemap_2020.tif')
age_Palmar <- raster('/home/stevie/Change_detection/PALMAR/agemap/Palmar_agemap_2020.tif')

plot(age_FCP, main = "FCP")
plot(age_Kiuic, main = "KK")
plot(age_Palmar, main = "EP")



plot(seq(0, 100), log(seq(0,100)), pch=20)
  
