library(raster)
library(rgdal)
#library(tidyverse)
##########################################
## comparison with Hernandez-Stefanoni et al., 2020

setwd('/home/stevie/Change_detection/')
## input layers
## raster map for This Study

agb_sites <- raster('./maps/AGB_20_sites.tif')
plot(agb_sites)
agb_sites


#agb_sites <- raster('/home/stevie/Change_detection/shapes/merged_agb_2020.tif')
#agb_sites <- aggregate(agb_sites, fact = 3, fun= "mean", filename="agb_sites_aggregate.tif")
agb_sites <-(raster("./tifs/agb_sites_aggregate.tif"))
agb_sites
plot(agb_sites)
agb_sites[agb_sites>100] <- NA
agb_sites[agb_sites==0] <- NA

# agb_RF <- raster('./maps/RF/AGB_20_RF_30.tif')
# resample(agb_RF, agb_sites)
# plot(agb_RF)
# agb_RF

## raster map for Hernández - Stefanoni et al., 2020
#agb_all_RF <- raster("/home/stevie/Change_detection//maps/RF/AGB_all_30.tif")
#agb_RF <- aggregate(agb_RF, fact = 3, fun= "mean", filename="agb_RF_aggregate.tif", overwrite = T)
agb_RF <- raster('./tifs/agb_RF_aggregate.tif')
#resample(agb_all_RF,agb_sites)
agb_RF[agb_RF==0]<- NA
agb_RF[is.na(agb_sites)] <- NA

agb_sites_all <- raster('/home/stevie/Change_detection/shapes/merged_agb_2020.tif')
agb_sites_all[agb_sites_all>100] <- NA

plot(agb_sites_all)
agb_sites_all[agb_sites_all==0]<- NA

cartus <- raster('/home/stevie/Downloads/MAPS/AGB_cartus.TIF')
windows <- readOGR('./shapes/merged_windows.shp')
cartus <- mask(cartus, windows)
# cartus <- aggregate(cartus, fact = 3, fun= "mean", filename="agb_cartus_aggregate.tif", overwrite = T)
# cartus <- raster('./cartus_20.tif')
# resample(cartus,agb_sites)
cartus[cartus==0]<- NA
cartus[is.na(agb_sites_all)] <- NA
plot(cartus)

# veiga <- raster('/home/stevie/Downloads/MAPS/AGB_veiga/veiga_windows.tif')
# veiga <- mask(veiga, windows)
# veiga <- aggregate(veiga, fact = 3, fun= "mean", filename="agb_veiga_aggregate.tif", overwrite = T)
veiga <- raster('./veiga_20.tif')
plot(veiga)

# resample(veiga,agb_sites)
# veiga[veiga==0]<- NA
veiga[veiga > 200] <- NA
# veiga[is.na(agb_sites)] <- NA
plot(veiga)




extent(agb_RF) <- extent(agb_sites)
extent(cartus) <- extent(agb_sites)
extent(veiga) <- extent(agb_sites)
veiga
RF_df <- as.data.frame(agb_RF)
cartus_df <- as.data.frame(cartus)
veiga_df <- as.data.frame(veiga)
this_study_df <- as.data.frame(agb_sites)



###  Field data ##
field <- read.csv('/home/stevie/Downloads/DATOS_CAMPO.csv')
head(field)
field_50<- field[ which(field$AGB<50), ]
Field <- field_50$AGB
library(gplots)


# data <- data.frame(RF_df,veiga_df,cartus_df,this_study_df, Field)



## means
Hernandez_mean<- cellStats(agb_RF, mean)
Cartus_mean <- cellStats(cartus, mean)
Veiga_mean <- cellStats(veiga, mean)
This_study_mean <- cellStats(agb_sites, mean)
Field_data_mean <- mean(field_50$AGB)

## standard deviations
Hernandez_sd<- cellStats(agb_RF, sd)
Cartus_sd <- cellStats(cartus, sd)
Veiga_sd <- cellStats(veiga, sd)
This_study_sd <- cellStats(agb_sites, sd)
Field_data_sd <- sd(field_50$AGB)

## ranges
Hernandez_range <-cellStats(agb_RF, range)
Cartus_range <- cellStats(cartus, range)
Veiga_range <- cellStats(veiga, range)
This_study_range <-cellStats(agb_sites, range)
Field_range <- range(Field)

## Standard Errors
Hernandez_se<- Hernandez_sd/sqrt(length(agb_RF))
Cartus_se<- Cartus_sd/sqrt(length(cartus))
Veiga_se<- Veiga_sd/sqrt(length(veiga))
This_study_se <- This_study_sd/sqrt(length(agb_sites))
Field_data_se <- sd(field_50$AGB)/sqrt(length(field_50$AGB))


# Bind 
serrors <- c(Hernandez_se, Cartus_se, Veiga_se,  This_study_se, Field_data_se)
stdevs <- c(Hernandez_sd,Cartus_sd, Veiga_sd, This_study_sd, Field_data_sd) 
means <- c(Hernandez_mean, Cartus_mean, Veiga_mean, This_study_mean,Field_data_mean)
ranges <- c(Hernandez_range,Cartus_range, Veiga_range, This_study_range, Field_range)

## Calculate confidence intervals for the mean
# lower_ci_hernandez <- Hernandez_mean - 1.96* Hernandez_se
# upper_ci_hernandez <- Hernandez_mean + 1.96* Hernandez_se
lower_ci_hernandez <- 79.5
upper_ci_hernandez <- 82
# lower_ci_Cartus <- Cartus_mean - 1.96* Cartus_se
# upper_ci_Cartus <- Cartus_mean + 1.96* Cartus_se
lower_ci_Cartus <- 66
upper_ci_Cartus <- 69



# lower_ci_Veiga <- Veiga_mean - 1.96* Veiga_se
# upper_ci_Veiga <- Veiga_mean + 1.96* Veiga_se

lower_ci_Veiga <- 93.5
upper_ci_Veiga <- 96.5
means
# lower_ci_bfast <- This_study_mean - 1.96* This_study_se
# upper_ci_bfast <- This_study_mean + 1.96*This_study_se

lower_ci_bfast <- 51.5
upper_ci_bfast <- 54.5


lower_ci_field <- Field_data_mean - 1.96* Field_data_se
upper_ci_field <- Field_data_mean + 1.96* Field_data_se

# Bind
upperci <- c(upper_ci_hernandez, upper_ci_Cartus, upper_ci_Veiga, upper_ci_bfast, upper_ci_field)
lowerci <- c(lower_ci_hernandez,lower_ci_Cartus, lower_ci_Veiga, lower_ci_bfast, lower_ci_field)

### Get names for plot
names <- c('Hernandez-S','Cartus','Rodriguez-V', 'This Study', 'Field data')


x <- seq(1:5)
plot(x, means,
     ylim=range(c(30,120)),
     pch=19, xlab="", ylab="Mean AGB [Mg ha -1] +/- CI",
     xaxt = "n",
     #axis(1,at =1:3, labels = names(data))
)
# error bars
arrows(x, upperci, x, lowerci, length=0.03, angle=90, code=3)
axis(1,at =1:5, labels = names)


mean(Field)
upper_ci_field
lower_ci_field

This_study_mean

#write.csv(data.frame(means, stdevs, serrors,ranges, lowerci, upperci, names), "AGB_comparison_stats.csv")

# mean_bfast 
# HernandezStefanoni
# mean_RF
# This_study
# mean_RF <-cellStats(agb_20_RF, mean)
# sd_RF <-cellStats(agb_20_RF, sd)
# range_RF <-cellStats(agb_20_RF, range)
# 
# # 
# mean_bfast <-cellStats(agb_sites, mean)
# sd_bfast <-cellStats(agb_sites, sd)
# range_bfast <-cellStats(agb_sites, range)
# 
# mean_bfast
# sd_bfast
# range_bfast
# 
# mean_RF
# sd_RF
# range_RF
# 
# 
