### This script performs counts of change categories and calculates Overall, Producer's and User's accuracy
### 12/06/2020   (2/2)


### Table must be set up according to following categories
## 0 = no change
## 1 = change (disturbance)     
## 3 = agriculture
## 4 = validation imagery unavailable (no imagery available for location/year, clouds, haze, etc)

          
#install.packages("plyr")
library(plyr)
library(dplyr)
library(rgdal)
  
setwd('/home/stevie/Change_detection/PALMAR/')
path <- getwd()
year = 2008
site = "Palmar"
### extension for files
ext = ".csv"
        
# 
# ### make filenames
val_data_ch <- read.csv(paste(path,"/","Val_change_", year,"_",site,".csv", sep=""))
val_data_nch <- read.csv(paste(path,"/","Val_nochange_", year,"_",site,".csv", sep=""))
# val_data_nch2 <- read.csv(paste(path,"/","Val_nochange2_", year,"_",site,".csv", sep=""))

# val_data_nch <- rbind(val_data_nch, val_data_nch2)

# filename= paste(path, infile, sep = "/")
# 
# ### read in data
# val_data <- read.csv(filename) 
# val_data_ch <- read.csv("/home/stevie/Change_detection/Kiuic/Val_change_Kiuic.csv")
# val_data_nch <- read.csv("/home/stevie/Change_detection/Kiuic/Val_nochange_Kiuic_15.csv")


#### subset data according to reference and estimated deforestation events
# 
positive_positive <- subset(val_data_ch, REF==1 & EST == 1)
false_positive <- subset(val_data_ch, REF==0 & EST == 1)
false_negative <- subset(val_data_nch, REF==1 & EST == 0)
negative_negative <- subset(val_data_nch, REF==0 & EST == 0)

# head(positive_positive)
# df <- rbind(positive_positive, false_positive, negative_negative, false_negative)
# df$ID <- 1:nrow(df)
# 
# change_class <- rbind(positive_positive,false_positive)
# nochange_class <- rbind(negative_negative, false_negative)
# 
# sample_change <- sample_n(df, 150)
# 
# valdf <- rbind(sample_change[-c(9)], nochange_class)

tru_positive <- subset(val_data_ch, REF==1 & EST == 1)
f_positive <- subset(val_data_ch, REF==0 & EST == 1)
f_negative <- subset(val_data_nch, REF==1 & EST == 0)
tru_negative <- subset(val_data_nch, REF==0 & EST == 0)


ch <- nrow(tru_positive)
ch_f <- nrow(f_positive)
nch <- nrow(tru_negative)
nch_f <- nrow(f_negative)

#wi_c <- round(wi_c, digits = 2)
#wi_nc <- round(wi_nc, digits = 2)

pi_c <- wi_c*(ch/(ch + ch_f))
pi_c_f <- wi_c*(ch_f/(ch + ch_f))
pi_nc <- wi_nc*(nch/(nch + nch_f))
pic_nc_f <-wi_nc*(nch_f/(nch + nch_f))

## Overall accuracy

Overall = pi_c + pi_nc

Users = pi_c /(pi_c + pi_c_f) 
Users_nc = pi_nc/(pi_nc + pic_nc_f)

Producers = pi_c/ (pi_c + pic_nc_f)
Producers_nc = pi_nc/(pi_c_f + pi_nc)



#####################################################################################################
### Error corrected estimates of ferorestation, error and confidence intervals
o_o <- negative_negative 
o_i <- false_positive
i_o <- false_negative
i_i  <- positive_positive

# 
defo_area <- total_area * pi_c
#defo_area_km2 <- defo_area/1000
#defo_area / forest_area
# 
change_area - defo_area
difference_ha <- (defo_area - change_area)/10000
# 
change_pos <- nrow(i_i)
change_false_pos <- nrow(o_i)
nchange <- change_pos + change_false_pos
nochange_pos <- nrow(o_o)
nochange_false_pos <- nrow(i_o)
nnochange <- nochange_pos + nochange_false_pos
# 
S = sqrt(wi_c^2  *(((change_pos/nchange) * (1- ((change_pos/nchange)))/(nchange-1))) + wi_nc^2*(((nochange_pos/nnochange) * (1- (nochange_pos/nnochange)))/(nnochange-1)))
# 
# #### compr
#S = sqrt(sum((.01^2*(((138/145)* (1-(138/145)))/144)),(.99^2*(((492/507)* (1-(492/507)))/506))))

SArea_ha = total_area * S
ic95 <- 1.96 * SArea_ha

# 
print(paste("error corrected change area =",round(defo_area, digits = 2), "+-", round(ic95, digits = 2), "ha"))
# install.packages('gfcanalysis')   

#difference = loweric - change_area

loweric <- defo_area - ic95
upperic <- defo_area + ic95

##
loweric

change_area

upperic


errors <- data.frame(loweric, defo_area, upperic)
write.csv(errors,paste(year, "errors", sep="_"))

# #------------------------------------------------------------------------------------------------#
# # Calculate variance and confidence intervals for O, U, P accuracies
# 
# get_overall_variance <- function (wi, users, pos_pos, f_pos) {VO = ((wi**2 *users) *(1 - users))/((pos_pos+f_pos) - 1)
# }
#   
# 
# change_variance <- get_overall_variance(wi_c, Users, nrow(positive_positive), nrow(false_positive))
# nochange_variance <- get_overall_variance(wi_nc, Users_nc,nrow(negative_negative), nrow(false_negative))
# variance_overall <- change_variance + nochange_variance
# 
# overall_ci =1.96*sqrt(variance_overall)
# SE <- sqrt(variance_overall)
# 
# Users_variance_c <- Users *(1 - Users)/((nrow(positive_positive) + nrow(false_positive)) - 1)
# users_ci_c =1.96*sqrt(Users_variance_c)
# SE_users_c <- sqrt(Users_variance_c)
# 
# 
# Users_variance_nc <- Users *(1 - Users)/((nrow(negative_negative) + nrow(false_negative)) - 1)
# users_ci_nc =1.96*sqrt(Users_variance_nc)
# SE_users_nc <- sqrt(Users_variance_nc)
# 
# nj_c = nrow(positive_positive) + nrow(false_negative)
# nj_nc = nrow(false_positive) + nrow(negative_negative)
# 
# ni_c = nrow(positive_positive) + nrow(false_positive)
# ni_nc = nrow(false_negative) + nrow(negative_negative)
# 
# Nj_c = round(defo_area/.09)
# Nj_nc = round((total_area - defo_area)/.09)
# 
# Ni_c = round(defo_area/.09)
# Ni_nc = round((total_area - defo_area)/.09)
# 
# 
# 
# #Njhat_c = round(Ni_c/ni_c*(nrow(positive_positive)))
# #Njhat_nc = round(Ni_nc/ni_nc*(nrow(negative_negative)))
# 
# ### Get total deforestation from GLOBAL TREE COVER PRODUCT (Sexton et al. 2013)
# ## add in Hansen's Gloval Tree Cover Product
# Kiuic_TC_2015 <-  raster('/home/stevie/Change_detection/KIUIC/tifs/Kiuic_TC_2015.tif')
# plot(Kiuic_TC_2015)
# Kiuic_TC_2015[Kiuic_TC_2015 < 1] <- NA
# 
# ## Deforestation = forest cover < 10 percent threshold (FAO, 2020) 
# Defo_Kiuic_15 <- Kiuic_TC_2015< 10
# plot(Defo_Kiuic_15)
# # mask out NF area
# Kiuic_FNF <- raster('/home/stevie/Change_detection/KIUIC/tifs/Kiuic_fnf.tif')
# plot(Kiuic_FNF)
# 
# Defo_Kiuic_15[Kiuic_FNF==0]<-NA
# Defo_Kiuic_15[Defo_Kiuic_15 ==0] <- NA
# 
# Defo_Kiuic_2015_poly <- rasterToPolygons(Defo_Kiuic_15)
# Defo_Kiuic_15_N <- length(Defo_Kiuic_2015_poly)
# 
# No_change_Kiuic_2015 <- Kiuic_FNF
# No_change_Kiuic_2015[Defo_Kiuic_15] <- NA
# 
# No_change_Kiuic_15_N <- length(No_change_Kiuic_2015)
# 
# plot(No_change_Kiuic_2015)
# plot(Defo_Kiuic_15, col = 'red', add = T)
# 
# Defo_area_Kiuic_2015 <-(length(Defo_Kiuic_2015_poly) * 900)/10000 ## Get deforestation area in ha
# 
# gArea(Kiuic_30_2000_poly)/10000 - change_area
# plot(magnthresh, add = T)
# 
# ### Get N of change and no change for reference
# 
# Njhat_c = Defo_Kiuic_15_N
# Njhat_nc = No_change_Kiuic_15_N
# 
# 
# 
# Producers_variance_c <- 
#   1/Njhat_c**2*(Nj_c**2*((1-Producers)**2)*(Users*(1-Users))
#         + Producers**2 *(
#           (Ni_c**2*(nrow(positive_positive)/ni_c)*(1-(nrow(positive_positive)/ni_c))/ni_c-1) 
#                           +(Ni_nc**2*(nrow(negative_negative)/ni_nc)*(1-(nrow(negative_negative)/ni_nc))/ni_nc-1)))
# 
# SE_producers_c = sqrt(Producers_variance_c)
# Producers_ci_c =1.96*sqrt(Producers_variance_c)
# 
# Producers_variance_nc <- 
#   1/Njhat_nc**2*(Nj_nc**2*((1-Producers_nc)**2)*(Users_nc*(1-Users_nc))
#                 + Producers_nc**2 *(
#                   (Ni_c**2*(nrow(positive_positive)/ni_c)*(1-(nrow(positive_positive)/ni_c))/ni_c-1) 
#                   +(Ni_nc**2*(nrow(negative_negative)/ni_nc)*(1-(nrow(negative_negative)/ni_nc))/ni_nc-1)))
# 
# SE_producers_nc = sqrt(Producers_variance_nc)
# Producers_ci_nc =1.96*sqrt(Producers_variance_nc)
# 
# accuracy <- data.frame(Overall,wi_c, 
#                        overall_ci,
#                        SE, 
#                        Users,
#                        SE_users_c, 
#                        users_ci_c,
#                        wi_nc, 
#                        Users_nc,
#                        SE_users_nc, 
#                        users_ci_nc, 
#                        Producers,
#                        SE_producers_c, 
#                        Producers_ci_c, 
#                        Producers_nc,
#                        SE_producers_nc, 
#                        Producers_ci_nc)
# 
# confu <- data.frame(pi_c, pi_c_f, pic_nc_f, pi_nc)


# write.csv(accuracy,paste("accuracy", year, site,ext, sep = "_"))
# write.csv(confu,paste("confu", year, site,ext, sep = "_"))

# # ## Overall Accuracy
# total = nrow(positive_positive) + nrow(false_positive) + nrow(false_negative) + nrow(negative_negative)
# overall = (nrow(positive_positive) + nrow(negative_negative))/ total
# #
# # ## disturbance
# producers = nrow(positive_positive) / (nrow(positive_positive) + nrow(false_positive))
# users = nrow(positive_positive) /(nrow(false_negative) + nrow(positive_positive))
# 
# 
# producers = nrow(positive_positive) / (nrow(positive_positive) + nrow(false_positive))
# users = nrow(positive_positive) /(nrow(false_negative) + nrow(positive_positive))
# 
# ###  no change
# producers_nochange = nrow(negative_negative) / (nrow(negative_negative) + nrow(false_negative))
# users_nochange = nrow(negative_negative) /(nrow(false_positive) + nrow(negative_negative))
# Accuracy <- data.frame(overall, producers, users, producers_nochange, users_nochange, year, site)
# 
# ### outfile
# outfile = paste("accuracy", year, site, sep = "_")
# outfile = paste(outfile, ext, sep = "")
# ## write results in table
# write.csv(Accuracy, file = outfile)
# #table_accuracy <- data.frame(Accuracy_Kiuic_2008, Accuracy_Kiuic_2011, Accuracy_Kiuic_2015, Accuracy_Kiuic_2015)
# #table <- t(table_accuracy)                            
# #write.csv(table, file="accuracy_table.csv")
# 
# 
# total = nrow(positive_positive) + nrow(false_positive) + nrow(false_negative) + nrow(negative_negative)
#     
# 
# confu_matrix <-data.frame(nrow(i_i), nrow(o_i), nrow(i_o), nrow(o_o))
# write.csv(confu_Kiuic_2015, file = "confu_Kiuic_2015.csv")
# 
#S= sum((wi_c^2*(nrow(tru_positive)/(nrow(tru_positive) + nrow(f_positive)))) 
#      + (wi_nc^2 * (nrow(tru_negative)/(nrow(tru_negative))+(nrow(f_negative)))))
#   
# pi = 0.06*(1034/1066) + 0.094*(129/415)
