######### Age/biomass relationships 
# install.packages("spatstat")
# library(spatstat)
library(raster)
library(rgdal)
site = "Kiuic"
setwd('/home/stevie/Change_detection/KIUIC/')

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

# get_agb <- function(age) {
#   exp(lm2$coefficients[1] + log(age)*lm2$coefficients[2])
# }

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

#########################################################################################################

# ## Produce error maps
# 
library(raster)  
agb_palmar <- raster('/home/stevie/Change_detection/PALMAR/tifs/palmar_agb.tif')
agb_fcp <- raster('/home/stevie/Change_detection/FCP/tifs/fcp_agb.tif')
agb_kiuic <- raster('/home/stevie/Change_detection/KIUIC/tifs/kiuic_agb.tif')

palmar_fnf <-  raster('/home/stevie/Change_detection/PALMAR/tifs/Palmar_fnf.tif')
kiuic_fnf <-  raster('/home/stevie/Change_detection/KIUIC/tifs/Kiuic_fnf.tif')
FCP_fnf <- raster('/home/stevie/Change_detection/FCP/tifs/FCP_fnf.tif')
# 
# 
agb_palmar[agb_palmar>100] <- NA
agb_palmar[agb_palmar==0] <- NA
# 
agb_kiuic[agb_kiuic>100] <- NA
agb_kiuic[agb_kiuic==0] <- NA
# 
agb_fcp[agb_fcp>120] <- NA
agb_fcp[agb_fcp==0] <- NA
# 


age_kiuic <- raster('/home/stevie/Change_detection/KIUIC/agemap/Kiuic_agemap_2020.tif')
age_kiuic_upper <- (2020 - age_kiuic) + 1.5
agb_palmar[agb_palmar==0] <- NA
lm2 <- lm(log(chrono_data_kiuic$ALL)~ log(chrono_data_kiuic$Edad.Estimada))
summary(lm2)

agb_kiuic_upper <- get_agb(age_kiuic_upper)
plot(agb_kiuic_upper)

extent(agb_kiuic_upper) <- extent(agb_kiuic)
error_kiuic <- ((agb_kiuic_upper - agb_kiuic)*100)/agb_kiuic
error_kiuic[error_kiuic< 0] <- NA

error_kiuic_ton <- agb_kiuic_upper - agb_kiuic
error_kiuic_ton[error_kiuic_ton< 0] <- NA
plot(error_kiuic_ton)


age_fcp <- raster('/home/stevie/Change_detection/FCP/agemap/FCP_agemap_2020.tif')
plot(age_fcp)
age_fcp_upper <- (2020 - age_fcp) + 1.5
plot(age_fcp_upper)
lm2 <- lm(log(chrono_data_fcp$Biomasa)~ log(chrono_data_fcp$Edad))

agb_fcp_upper <- get_agb(age_fcp_upper)
plot(agb_fcp_upper)

extent(agb_fcp_upper) <- extent(agb_fcp)
error_fcp <- ((agb_fcp_upper - agb_fcp)*100)/agb_fcp
error_fcp[error_fcp< 0] <- NA

error_fcp_ton <- agb_fcp_upper - agb_fcp
error_fcp_ton[error_fcp_ton< 0] <- NA

plot(error_fcp_ton)


plot(agb_fcp)
age_palmar <- raster('/home/stevie/Change_detection/PALMAR/agemap/Palmar_agemap_2020.tif')

age_palmar_upper <- (2020 - age_palmar) + 1.5
lm2 <- lm(log(chrono_data_kiuic$ALL)~ log(chrono_data_kiuic$Edad.Estimada))

agb_palmar_upper <- get_agb(age_palmar_upper)
plot(agb_palmar_upper)

extent(agb_palmar_upper) <- extent(agb_palmar)
error_palmar <- ((agb_palmar_upper - agb_palmar)*100)/agb_palmar
error_palmar[error_palmar< 0] <- NA
error_palmar_ton <- agb_palmar_upper - agb_palmar
error_palmar_ton[error_palmar_ton< 0] <- NA

plot(error_palmar_ton)

#time_lag_agb_fcp <- get_agb(1.5) # example





### Map figures
library(RColorBrewer)
par(mar=c(2.5,3,1,.5),mfrow=c(2,2))

## Palmar

#plot(lower_palmar,col =brewer.pal(9,"Blues"), colNA = "light gray", cex.axis =1.5, legend= T, cex = 1.5)
#plot(palmar_fnf, col = c("black","transparent"), add = T, legend= F)
#plot(agb_palmar, colNA = "light gray", cex.axis =1.5)
#plot(palmar_fnf, col = c("black","transparent"), add = T, legend= F)
# plot(error_palmar,col =brewer.pal(9,"Reds"),colNA = "light gray", cex.axis =1.5)
# plot(palmar_fnf, col = c("black","transparent"), add = T, legend= F)
# mtext("a)",side = 3, adj = 0.02, line = -1.5, cex =1.5)

## error in ton/ha

plot(error_palmar_ton,colNA = "light gray", cex.axis =1.5)
plot(palmar_fnf, col = c("black","transparent"), add = T, legend= F)
mtext("a)",side = 3, adj = 0.02, line = -3, cex =1.5)


# Kiuic
#plot(lower_kiuic,col =brewer.pal(9,"Reds"), colNA = "light gray", cex.axis =1.5)
#plot(kiuic_fnf, col = c("black","transparent"), add = T, legend= F)
#plot(agb_kiuic, colNA = "light gray", cex.axis =1.5)
#plot(kiuic_fnf, col = c("black","transparent"), add = T, legend= F)
# plot(error_kiuic,col =brewer.pal(9,"Reds"), colNA = "light gray", cex.axis =1.5)
# plot(kiuic_fnf, col = c("black","transparent"), add = T, legend= F)
# mtext("b)",side = 3, adj = 0.02, line = -1.5, cex =1.5)

# error in ton/ha

plot(error_kiuic_ton, colNA = "light gray", cex.axis =1.5)
plot(kiuic_fnf, col = c("black","transparent"), add = T, legend= F)
mtext("b)",side = 3, adj = 0.02, line = -3, cex =1.5)


## FCP

#plot(lower_fcp,col =brewer.pal(9,"Reds"), colNA = "light gray", cex.axis =1.5)
#plot(FCP_fnf, col = c("black","transparent"), add = T, legend= F)
#plot(agb_fcp, colNA = "light gray", cex =1.5, cex.axis =1.5)
#plot(FCP_fnf, col = c("black","transparent"), add = T, legend= F)
# plot(error_fcp,col =brewer.pal(9,"Reds"), colNA = "light gray", cex.axis =1.5)
# plot(FCP_fnf, col = c("black","transparent"), add = T, legend= F)
# mtext("c)",side = 3, adj = 0.02, line = -1.5, cex =1.5)

## error in ton/ha

plot(error_fcp_ton, colNA = "light gray", cex.axis =1.5)
plot(FCP_fnf, col = c("black","transparent"), add = T, legend= F)
mtext("c)",side = 3, adj = 0.02, line = -3, cex =1.5)


#Histograms
par(mar=c(4,4,1,1),mfrow=c(2,2))
hist(error_palmar_ton, main="Deciduous", xlab="AGB Ton/ha")
hist(error_kiuic_ton, main = "Semi-deciduous", xlab="AGB Ton/ha")
hist(error_fcp_ton, main="Semi-evergreen", xlab="AGB Ton/ha")



sd_error_fcp <- cellStats(error_fcp_ton, stat = sd)

sd_error_kiuic <- cellStats(error_kiuic_ton, stat = sd)

sd_error_palmar <- cellStats(error_palmar_ton, stat = sd)


mean_error_fcp<- cellStats(error_fcp_ton, stat = mean)

mean_error_kiuic<- cellStats(error_kiuic_ton, stat = mean)

mean_error_palmar <- cellStats(error_palmar_ton, stat = mean)

# Cartus_se<- Cartus_sd/sqrt(length(cartus))
se_error_kiuic <- sd_error_kiuic/sqrt(length(error_kiuic_ton))

se_error_palmar <- sd_error_palmar/sqrt(length(error_palmar_ton))

se_error_fcp <- sd_error_fcp/sqrt(length(error_fcp_ton))


upper_ci_error_kiuic <- mean_error_kiuic + 1.96 * se_error_kiuic
lower_ci_error_kiuic  <- mean_error_kiuic - 1.96 * se_error_kiuic

upper_ci_error_palmar <- mean_error_palmar + 1.96 * se_error_palmar
lower_ci_error_palmar  <- mean_error_palmar - 1.96 * se_error_palmar

upper_ci_error_fcp <- mean_error_fcp + 1.96 * se_error_fcp
lower_ci_error_fcp  <- mean_error_fcp - 1.96 * se_error_fcp


# Palmar
upper_ci_error_palmar
mean_error_palmar
lower_ci_error_palmar


# Kiuic
upper_ci_error_kiuic
mean_error_kiuic
lower_ci_error_kiuic

# FCP
upper_ci_error_fcp
mean_error_fcp
lower_ci_error_fcp




# # 
# # writeRaster(agb_sites, "/home/stevie/Change_detection//maps/AGB_20_sites", format="GTiff")
# # writeRaster(agb_20_RF, "/home/stevie/Change_detection/maps/RF/AGB_20_RF_30.tif", format="GTiff")
# 
# ###############################################
# ### Age simulations 
# 
# # Simulate agb for ages 0 - 150 years
# age_sim <- seq(0,350)
# agb_sim <- get_agb(seq(0,350))
# 
# df <- data.frame(age_sim, agb_sim)
# 
# plot(df$age_sim, df$agb_sim,type="l", col="red")


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

