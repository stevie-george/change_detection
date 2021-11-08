######### Age/biomass relationships 
# install.packages("spatstat")
# library(spatstat)
library(raster)
library(rgdal)
site = "FCP"
setwd('/home/stevie/Change_detection/FCP/')

#### Part 1. Get yearmaps 
### get yearmaps ####

# setwd("./tifs/sieves/")
# files <- list.files(getwd())
# 
# years <- seq(2000,2020)
# get_yearmaps <- function(files, years) {
#   for (i in 1:(length(files))) {
#     mask1 <- raster(files[i])
#     mask1 <- (mask1/mask1) + years[i]-1
#     plot(mask1)
#     writeRaster(mask1, filename = paste("yearmap",years[i], sep = "_"), format="GTiff")
#   }}
# 
# get_yearmaps(files, years)

# get oldest year from FNF mask

# year21 <- raster('/home/stevie/Change_detection/FCP/tifs/FCP_fnf.tif')
# year21[year21==0] <- NA
# year21 = year21 + 1998
# plot(year21)
# writeRaster(year21, filename = "year_21", format ="GTiff")
# 



### Part 2. stitch maps together
## In QGIS use r.patch function and patch years 2020, 2019, 2018 to fill in NA's from younguest to oldest



### Part 3. Get AGE - BIOMASS model from chronosequence
## read in chronosequence data for obtaining model
chrono_data_kiuic<- read.csv("/home/stevie/CARDAMOM/chrono/Kiuic_Chrono_2010_DALEC.csv")

chrono_data_fcp <- read.csv("/home/stevie/Change_detection/FCP/QROO_oct_13_BIOMASA.csv")

## Get AGB from age - biomass relationship obtained from chronosequence mode.
# FCP
chrono_data_fcp <- chrono_data_fcp[ which(chrono_data_fcp$Edad < 100), ]
# plot(chrono_data_fcp$Edad, chrono_data_fcp$Biomasa, pch = 20, xlab = "Stand age", ylab ="AGB")

lm1 <- lm((chrono_data_fcp$Biomasa)~ log(chrono_data_fcp$Edad))
summary(lm1)
lm2 <- lm(log(chrono_data_fcp$Biomasa)~ log(chrono_data_fcp$Edad))
summary(lm2)
lm3 <- lm(log(chrono_data_fcp$Biomasa)~ (chrono_data_fcp$Edad))
summary(lm3)
lm4 <- lm((chrono_data_fcp$Biomasa)~ log(chrono_data_fcp$Edad))
summary(lm4)
lm5 <- lm(sqrt(chrono_data_fcp$Biomasa)~ sqrt(chrono_data_fcp$Edad))
summary(lm5)

##### Get AGB from age - biomass relationship obtained from chronosequence
### function to obtain AGB from age using a log-transformed agb-age relationship

get_agb <- function(age) {
  exp(lm2$coefficients[1] + log(age)*lm2$coefficients[2])
}

# Plot site 1.
par(mfrow=c(1,2), mar=c(4.5,4.5,4,1.5))

plot(chrono_data_fcp$Edad, (chrono_data_fcp$Biomasa),main="Semi-evergreen tropical forest", pch=20, xlab = "Stand age", ylab ="AGB [Mg ha-1]", cex.lab = 1.5, cex.axis = 1.2)
mtext("FCP", cex = 1.5)
mtext("Rsq = 0.76", side = 1, line = -1)
age_sim <- seq(0,150)
agb_sim <- get_agb(seq(0,150))
df <- data.frame(age_sim, agb_sim)
lines(df$age_sim, df$agb_sim, lty="dashed", col="blue", lwd = 3)


# kiuic & Palmar - read in Kiuic chrono
head(chrono_data_kiuic)

lm1 <- lm(chrono_data_kiuic$ALL~ log(chrono_data_kiuic$Edad.Estimada))
summary(lm1)
lm2 <- lm(log(chrono_data_kiuic$ALL)~ log(chrono_data_kiuic$Edad.Estimada))
summary(lm2)
lm3 <- lm(log(chrono_data_kiuic$ALL)~ (chrono_data_kiuic$Edad.Estimada))
summary(lm3)
lm4 <- lm((chrono_data_kiuic$ALL)~ log(chrono_data_kiuic$Edad.Estimada))
summary(lm4)
lm5 <- lm(sqrt(chrono_data_kiuic$ALL)~ sqrt(chrono_data_kiuic$Edad.Estimada))
summary(lm5)


### plot site 2. 
plot(chrono_data_kiuic$Edad.Estimada,chrono_data_kiuic$ALL ,main="Semi-deciduous tropical forest", pch=20, xlab= "Stand age", ylab= "", cex.axis = 1.2, cex.lab = 1.5)
mtext("KK, EP", cex = 1.5)
mtext("Rsq = 0.52", side = 1, line = -1)
age_sim <- seq(0,150)
agb_sim <- get_agb(seq(0,150))
df <- data.frame(age_sim, agb_sim)
lines(df$age_sim, df$agb_sim, lty="dashed", col="green", lwd = 3)
#########################################################################################################
## Calculate AGB from age raster from agb-age relationship obtained from chronosequence data 

# ### function to obtain AGB from age using a log-transformed agb-age relationship
# 
# get_agb <- function(age) {
#   exp(lm2$coefficients[1] + log(age)*lm2$coefficients[2])
# }


# # introduce the year of forest loss raster
# age_year <- raster('../../agemap/FCP_agemap_2020.tif')
#     plot(age_year)
# age <- 2020 - age_year 
# plot(age)
# 
# ### apply function to age raster
# age_agb <- get_agb(age)
# plot(age_agb, main = site)

# write AGB raster
#writeRaster(age_agb, filename=paste(site,"AGB", sep = "_"), format ="GTiff", overwrite =T)


##########################################
## comparison with Hernandez-Stefanoni et al., 2020
agb_sites
agb_20_RF
## input layers
par(mfrow=c(1,2))
agb_sites <- raster('/home/stevie/Change_detection/shapes/merged_agb_2020.tif')
agb_sites <- aggregate(agb_sites, fact = 3, fun= "mean", filename="agb_sites_aggregate.tif")
plot(agb_sites)
agb_sites[agb_sites>100] <- NA
agb_sites[agb_sites==0] <- NA

agb_all_RF <- raster("/home/stevie/Change_detection//maps/RF/AGB_all_30.tif")
agb_all_RF <- aggregate(agb_all_RF, fact = 3, fun= "mean", filename="agb_RF_aggregate.tif", overwrite = T)
agb_all_RF[agb_all_RF==0]<- NA
agb_20_RF <- agb_all_RF * (agb_sites/agb_sites)
extent(agb_20_RF) <- extent(agb_sites)
plot(agb_20_RF)

# plotting 
extent(agb_20_RF) <- extent(agb_sites)
s <- stack(agb_20_RF, agb_sites)
stack.csv <- as.data.frame(s)


names(s) <- c('Random Forest', 'BFast')
par(mfrow=c(1,1))
extent(agb_20_RF) <- extent(agb_sites)



###  Field data ##
field <- read.csv('/home/stevie/Downloads/ModCamp_pred.csv')
head(field)
field_50<- field[ which(field$obs<50), ]
field_50$study <- 'Hernandez-Stefanoni'

field_50

Field <- field_50$obs
library(gplots)
## NO
#boxplot(s,col=c('lightsteelblue2', 'lightsteelblue4'), ylab='AGB [Mg ha -1]' )
Field_50 <- data.frame(stack.csv,Field)
names(Field_50) <- c('Hernandez-Stefanoni', 'This Study', 'Field data')

boxplot(Field_50,col=c('lightsteelblue2', 'lightsteelblue4', 'gray'), ylab='AGB [Mg ha -1]' )
head(Field_50)
tail(Field_50)

HernandezStefanoni<- mean(Field_50$layer)
This_study <- mean(Field_50$merged_agb_2020)
Field_data <- mean(Field_50$Field)

HernandezStefanoni_sd<- sd(Field_50$layer)
This_study_sd <- sd(Field_50$merged_agb_2020)
Field_data_sd <- sd(Field_50$Field)

HernandezStefanoni_se<- sd_RF/sqrt(length(Field_50$layer))
This_study_se <- sd_bfast/sqrt(length(Field_50$merged_agb_2020))
Field_data_se <- sd(Field_50$Field)/sqrt(length(Field))

serrors <- c(HernandezStefanoni_se, This_study_se, Field_data_se)
stdevs <- c(sd_RF, sd_bfast, Field_data_sd) 
means <- c(mean_RF, mean_bfast,Field_data)
means

upper_ci_hernandez
mean_RF
lower_ci_hernandez

upper_ci_bfast
mean_bfast
lower_ci_bfast

lower_ci_field
mean(field_50$obs)
upper_ci_field

lower_ci_field
lower_ci_hernandez <- mean_RF - 1.96*(sd_RF/sqrt(length(agb_20_RF)))
upper_ci_hernandez <- mean_RF + 1.96*(sd_RF/sqrt(length(agb_20_RF)))

lower_ci_bfast <- mean_bfast - 1.96*(sd_bfast/sqrt(length(agb_sites)))
upper_ci_bfast <- mean_bfast + 1.96*(sd_bfast/sqrt(length(agb_sites)))

lower_ci_field <- Field_data - 1.96*(Field_data_sd/sqrt(length(Field_50)))
upper_ci_field <- Field_data + 1.96*(Field_data_sd/sqrt(length(Field_50)))

lowerci
upperci <- c(upper_ci_hernandez, upper_ci_bfast, upper_ci_field)
lowerci <- c(lower_ci_hernandez, lower_ci_bfast, lower_ci_field)

names(Field_50) <- c("Hernandez-Stefanoni", "This study", "Field data")
x <- 1:length(Field_50)
plot(x, means,
     ylim=range(c(min(lowerci), max(upperci))),
     pch=19, xlab="", ylab="Mean AGB [Mg ha -1] +/- CI",
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, upperci, x, lowerci, length=0.1, angle=90, code=3)




mean_RF <-cellStats(agb_20_RF, mean)
sd_RF <-cellStats(agb_20_RF, sd)
range_RF <-cellStats(agb_20_RF, range)

# 
mean_bfast <-cellStats(agb_sites, mean)
sd_bfast <-cellStats(agb_sites, sd)
range_bfast <-cellStats(agb_sites, range)

mean_bfast
sd_bfast
range_bfast

mean_RF
sd_RF
range_RF



library(raster)  
agb_palmar <- raster('/home/stevie/Change_detection/PALMAR/tifs/palmar_agb.tif')
agb_fcp <- raster('/home/stevie/Change_detection/FCP/tifs/fcp_agb.tif')
agb_kiuic <- raster('/home/stevie/Change_detection/KIUIC/tifs/kiuic_agb.tif')


agb_palmar[agb_palmar>100] <- NA
agb_palmar[agb_palmar==0] <- NA

agb_kiuic[agb_kiuic>100] <- NA
agb_kiuic[agb_kiuic==0] <- NA

agb_fcp[agb_fcp>120] <- NA
agb_fcp[agb_fcp==0] <- NA



round(cellStats(agb_fcp, stat = sd), digits = 2)

round(cellStats(agb_kiuic, stat = sd), digits = 2)

round(cellStats(agb_palmar, stat = sd), digits = 2)


round(cellStats(agb_fcp, stat = mean), digits = 2)

round(cellStats(agb_kiuic, stat = mean), digits = 2)

round(cellStats(agb_palmar, stat = mean), digits = 2)





# 
# writeRaster(agb_sites, "/home/stevie/Change_detection//maps/AGB_20_sites", format="GTiff")
# writeRaster(agb_20_RF, "/home/stevie/Change_detection/maps/RF/AGB_20_RF_30.tif", format="GTiff")
  
###############################################
### Age simulations 

get_agb(seq(0,100)) # example

# Simulate agb for ages 0 - 150 years
age_sim <- seq(0,350)
agb_sim <- get_agb(seq(0,350))

df <- data.frame(age_sim, agb_sim)

plot(df$age_sim, df$agb_sim,type="l", col="red")


#### in qgis patch from younger to older r.patch = 2020, 2019, 2018, etc
# 
# get_agemap <- function(files) {
#   for (i in 1:(length(files))) {
#     tif1 <- files[i]
#     raster1 <- raster(tif1)
#     tif2 <- files[i]
#     raster2 <- raster(tif2)
#     rast <- cover(raster1, raster2)
#     writeRaster(rast,filename = "test2", format="GTiff")
#   }
# }


# ### validation
# 
# agb.db <- read.csv("/home/stevie/Change_detection/agb_validation_infys.csv")
# head(agb.db)
# 
# fcp.db <- agb.db[ which(agb.db$fcp_agb>0), ]
# palmar.db <- agb.db[ which(agb.db$palmar_agb>0), ]
# kiuic.db <- agb.db[ which(agb.db$kiuic_agb>0), ]
# 
# par(mfrow=c(1,3))
# 
# 
# fcp_age_val<-lm(fcp.db$AGB ~fcp.db$fcp_agb)
# summary(fcp_age_val)
# 

