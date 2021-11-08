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
resample(agb_all_RF,agb_sites)
agb_RF[agb_RF==0]<- NA
agb_RF[is.na(agb_sites)] <- NA

agb_RF
agb_sites

extent(agb_RF) <- extent(agb_sites)
maps <- stack(agb_RF, agb_sites)
maps_df <- as.data.frame(maps)

###  Field data ##
field <- read.csv('/home/stevie/Downloads/DATOS_CAMPO.csv')
head(field)
field_50<- field[ which(field$AGB<50), ]
Field <- field_50$AGB
library(gplots)

data <- data.frame(maps_df, Field)

### Get names for plot
names(data) <- c('Hernandez-Stefanoni', 'This Study', 'Field data')

## Boxplot
boxplot(data,col=c('lightsteelblue2', 'lightsteelblue4', 'gray'), ylab='AGB [Mg ha -1]' )


## means
Hernandez_mean<- cellStats(agb_RF, mean)
This_study_mean <- cellStats(agb_sites, mean)
Field_data_mean <- mean(field_50$AGB)

## standard deviations
Hernandez_sd<- cellStats(agb_RF, sd)
This_study_sd <- cellStats(agb_sites, sd)
Field_data_sd <- sd(field_50$AGB)

## ranges
Hernandez_range <-cellStats(agb_RF, range)
This_study_range <-cellStats(agb_sites, range)
Field_range <- range(Field)

## Standard Errors
Hernandez_se<- Hernandez_sd/sqrt(length(agb_RF))
This_study_se <- This_study_sd/sqrt(length(agb_sites))
Field_data_se <- sd(field_50$AGB)/sqrt(length(field_50$AGB))

  
# Bind 
serrors <- c(Hernandez_se, This_study_se, Field_data_se)
stdevs <- c(Hernandez_sd, This_study_sd, Field_data_sd) 
means <- c(Hernandez_mean, This_study_mean,Field_data_mean)
ranges <- c(Hernandez_range, This_study_range, Field_range)

## Calculate confidence intervals for the mean
lower_ci_hernandez <- Hernandez_mean - 1.96* Hernandez_se
upper_ci_hernandez <- Hernandez_mean + 1.96* Hernandez_se

lower_ci_bfast <- This_study_mean - 1.96* This_study_se
upper_ci_bfast <- This_study_mean + 1.96*This_study_se

lower_ci_field <- Field_data_mean - 1.96* Field_data_se
upper_ci_field <- Field_data_mean + 1.96* Field_data_se

# Bind
upperci <- c(upper_ci_hernandez, upper_ci_bfast, upper_ci_field)
lowerci <- c(lower_ci_hernandez, lower_ci_bfast, lower_ci_field)


names(data) <- c("Hernandez-Stefanoni", "This study", "Field data")
x <- 1:length(data)
plot(x, means,
     ylim=range(c(30, 100)),
     pch=19, xlab="", ylab="Mean AGB [Mg ha -1] +/- CI",
     xaxt = "n",
     #axis(1,at =1:3, labels = names(data))
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, upperci, x, lowerci, length=0.1, angle=90, code=3)
axis(1,at =1:3, labels = names(data))

write.csv(data.frame(means, stdevs, serrors,ranges, lowerci, upperci), "AGB_comparison_stats.csv")

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
