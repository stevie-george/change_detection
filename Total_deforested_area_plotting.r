library(raster)

######## Total area calculation
path_to_sieves = '/home/stevie/Change_detection/FCP/tifs/sieves/'
setwd(path_to_sieves)

magn_00 <- raster('magn00areasieve.tif')
mask_00 = magn_00/magn_00
total_00 = (sum(mask_00[mask_00==1])*900)/10000 


magn_01 <- raster('magn01areasieve.tif')
mask_01 = magn_01/magn_01
total_01 = (sum(mask_01[mask_01==1])*900)/10000 


magn_02 <- raster('magn02areasieve.tif')
mask_02 = magn_02/magn_02
total_02 = (sum(mask_02[mask_02==1])*900)/10000 


magn_03 <- raster('magn03areasieve.tif')
mask_03 = magn_03/magn_03
total_03 = (sum(mask_03[mask_03==1])*900)/10000 


magn_04 <- raster('magn04areasieve.tif')
mask_04 = magn_04/magn_04
total_04 = (sum(mask_04[mask_04==1])*900)/10000 


magn_05 <- raster('magn05areasieve.tif')
mask_05 = magn_05/magn_05
total_05 = (sum(mask_05[mask_05==1])*900)/10000 


magn_06 <- raster('magn06areasieve.tif')
mask_06 = magn_06/magn_06
total_06 = (sum(mask_06[mask_06==1])*900)/10000 


magn_07 <- raster('magn07areasieve.tif')
mask_07 = magn_07/magn_07
total_07 = (sum(mask_07[mask_07==1])*900)/10000 


magn_08 <- raster('magn08areasieve.tif')
mask_08 = magn_08/magn_08
total_08 = (sum(mask_08[mask_08==1])*900)/10000 


magn_09 <- raster('magn09areasieve.tif')
mask_09 = magn_09/magn_09
total_09 = (sum(mask_09[mask_09==1])*900)/10000 


magn_10 <- raster('magn10areasieve.tif')
mask_10 = magn_10/magn_10
total_10 = (sum(mask_10[mask_10==1])*900)/10000 


magn_11 <- raster('magn11areasieve.tif')
mask_11 = magn_11/magn_11
total_11 = (sum(mask_11[mask_11==1])*900)/10000 


magn_12 <- raster('magn12areasieve.tif')
mask_12 = magn_12/magn_12
total_12 = (sum(mask_12[mask_12==1])*900)/10000 


magn_13 <- raster('magn13areasieve.tif')
mask_13 = magn_13/magn_13
total_13 = (sum(mask_13[mask_13==1])*900)/10000 


magn_14 <- raster('magn14areasieve.tif')
mask_14 = magn_14/magn_14
total_14 = (sum(mask_14[mask_14==1])*900)/10000 


magn_15 <- raster('magn15areasieve.tif')
mask_15 = magn_15/magn_15
total_15 = (sum(mask_15[mask_15==1])*900)/10000 


magn_16 <- raster('magn16areasieve.tif')
mask_16 = magn_16/magn_16
total_16 = (sum(mask_16[mask_16==1])*900)/10000 


magn_17 <- raster('magn17areasieve.tif')
mask_17 = magn_17/magn_17
total_17 = (sum(mask_17[mask_17==1])*900)/10000 


magn_18 <- raster('magn18areasieve.tif')
mask_18 = magn_18/magn_18
total_18 = (sum(mask_18[mask_18==1])*900)/10000 
  

magn_19 <- raster('magn19areasieve.tif')
mask_19 = magn_19/magn_19
total_19 = (sum(mask_19[mask_19==1])*900)/10000 


magn_20 <- raster('magn20areasieve.tif')
mask_20 = magn_20/magn_20
total_20 = (sum(mask_20[mask_20==1])*900)/10000 


time = c(2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
         2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020)
total <- c(total_00, total_01, total_02, total_03, total_04,  total_05,  total_06,  total_07,  total_08,
           total_09,  total_10,  total_11,  total_12,  total_13,  total_14,  total_15,  total_16,
           total_17,  total_18,  total_19,  total_20)

df <- data.frame(time, total)
library(ggplot2)



#total_change_fcp <- df

#total_change_kiuic <- df
total_change_palmar <- df
total_change_palmar <- df[-c(21),]





## plot change area per site
 ggplot(total_change_palmar, aes(x = time, y = total)) +  labs(x="Time", y="Total change area (ha)") +                          # bar plot
   geom_col(size= 1, colour = "black") +
   geom_line(size=1, colour ="dark blue")

 ggplot(total_change_kiuic, aes(x = time, y = total)) +  labs(x="Time", y="Total change area (ha)") +                          # bar plot
   geom_col(size= 1, colour = "black") +
   geom_line(size=1, colour ="dark red")

 ggplot(total_change_fcp, aes(x = time, y = total)) +  labs(x="Time", y="Total change area (ha)") +                          # bar plot
   geom_col(size= 1, colour = "black") +
   geom_line(size=1, colour ="purple")
# 



## plot change area for all sites


  

 ggplot() + geom_line(data = total_change_kiuic, aes(x = time, y = total, colour = "KK"))  +
   geom_line(data =total_change_palmar, aes(x =time, y = total, colour = "EP")) +
   geom_line(data =total_change_fcp, aes(x=time, y=total, colour = "FCP"), na.rm = T, stat="identity") + 
   labs(colour = "",x ="Year of change", y ="Total area of change (ha)") +
   theme_light()
 

# ### cumsum
#   
# cumsum_kiuic <- cumsum(total_change_kiuic)
# cumsum_fcp <- cumsum(total_change_fcp)
# cumsum_palmar <- cumsum(total_change_palmar)
# 
# cdf_kiuic <- data.frame(cumsum_kiuic, time) 
# cdf_fcp <-data.frame(cumsum_fco, time)
# cdf_palmar <- data.frame(cumsum_palmar, time[1:20]) 
# 
# 
# 
# ggplot() + geom_line(data = cumsum_kiuic, aes(x = time, y = total, colour = "KK"))  +
#   geom_line(data =cumsum_palmar, aes(x =time, y = total, colour = "EP")) +
#   geom_line(data =cumsum_fcp, aes(x=time, y=total, colour = "FCP"), na.rm = T, stat="identity") + 
#   labs(colour = "",x ="Year of change", y ="Total area of change (ha)") +
#   theme_light()
# 
# 
# 

