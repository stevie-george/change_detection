# FCP site 
### 

site = 'FCP'
defo_perc = 50
list = c("25", "50", "75", "90", "99")
# Output the desired time periods
#par(mfrow=c(5,5), mar=c(1,1,1,1))
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1))
load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp585_MOHC_1899_2100_site_',site,'_',defo_perc,'.RData', sep = ""))    
states_all_FCP_site_FCP_met_585 <- states_all
load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp434_MOHC_1899_2100_site_',site,'_',defo_perc,'.RData', sep = ""))
states_all_FCP_site_FCP_met_434 <- states_all
load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp370_MOHC_1899_2100_site_',site,'_',defo_perc,'.RData', sep = ""))    
states_all_FCP_site_FCP_met_370 <- states_all
load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp245_MOHC_1899_2100_site_',site,'_',defo_perc,'.RData', sep = ""))    
states_all_FCP_site_FCP_met_245 <- states_all
load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp126_MOHC_1899_2100_site_',site,'_',defo_perc,'.RData', sep = ""))
states_all_FCP_site_FCP_met_126 <- states_all
getwd()
#defo_perc = 25

get_def <- function(defo_perc) {
    load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp585_MOHC_1899_2100_site_FCP_',defo_perc,'.RData', sep = ""))    
    states_all_FCP_site_FCP_met_585 <- states_all
    dat <- states_all_FCP_site_FCP_met_585$wood_gCm2[,1430:2424]
    for (i in 1:(length(dat))) {
    rate_585 = ((dat[1] - dat[994])/ length(dat)) * 100 }
    load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp434_MOHC_1899_2100_site_FCP_',defo_perc,'.RData', sep = ""))  
    states_all_FCP_site_FCP_met_434 <- states_all
    dat <-states_all_FCP_site_FCP_met_434$wood_gCm2[,1430:2424]
    for (i in 1:(length(dat))) {
    rate_434 = ((dat[1] - dat[994])/ length(dat)) * 100 }
    load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp370_MOHC_1899_2100_site_FCP_',defo_perc,'.RData', sep = ""))    
    states_all_FCP_site_FCP_met_370 <- states_all
    dat <- states_all_FCP_site_FCP_met_370$wood_gCm2[,1430:2424]
    for (i in 1:(length(dat))) {
      rate_370 = ((dat[1] - dat[994])/ length(dat)) * 100 }
    load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp245_MOHC_1899_2100_site_FCP_',defo_perc,'.RData', sep = ""))    
    states_all_FCP_site_FCP_met_245 <- states_all
    dat <- states_all_FCP_site_FCP_met_245$wood_gCm2[,1430:2424]
    for (i in 1:(length(dat))) {
      rate_245 = ((dat[1] - dat[994])/ length(dat)) * 100 }
    load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp126_MOHC_1899_2100_site_FCP_',defo_perc,'.RData', sep = ""))
    states_all_FCP_site_FCP_met_126 <- states_all
    dat <- states_all_FCP_site_FCP_met_126$wood_gCm2[,1430:2424]
    for (i in 1:(length(dat))) {
      rate_126 = ((dat[1] - dat[994])/ length(dat)) * 100 }
  rates <- data.frame(rate_126, rate_245, rate_370, rate_434, rate_585)
write.csv(rates, paste("rates", defo_perc,".csv", sep = ""))
}


for (i in 1:(length(list))) {
  get_def(list[i])
}

rates_25 <- read.csv('./rates25.csv')
rates_50 <- read.csv('./rates50.csv')
rates_75 <- read.csv('./rates75.csv')
rates_90 <- read.csv('./rates90.csv')
rates_99 <- read.csv('./rates99.csv')


rates_all <- rbind(rates_25, rates_50, rates_75, rates_90, rates_99)

par(mfrow=c(1,1))
plot(rates_all$rate_585, type = "l", col = "red", xlab = "% AGB removal",ylab ="Growth rate",  xaxt="n", main = paste(site, "growth rates"))
legend("topleft",bg='transparent', c("scenario 126","scenario 245", " scenario 370", "scenario 434", "scenario 585"),cex=1,col=c("blue","green", "yellow","orange", "red"), pch = 15, box.lty = 0)
axis(1, at=1:5, labels=list, cex.axis=1)
lines(rates_all$rate_434, type = "l", col = "orange", xlab = "% AGB removal", ylab ="Growth rate",  xaxt="n")
#axis(1, at=1:5, labels=list, cex.axis=1)
lines(rates_all$rate_370, type = "l", col = "yellow", xlab = "% AGB removal", ylab ="Growth rate",  xaxt="n")
#axis(1, at=1:5, labels=list, cex.axis=1)
lines(rates_all$rate_245, type = "l", col = "green", xlab = "% AGB removal", ylab ="Growth rate",  xaxt="n")
#axis(1, at=1:5, labels=list, cex.axis=1)
lines(rates_all$rate_126, type = "l", col = "blue", xlab = "% AGB removal",ylab ="Growth rate",  xaxt="n")
#axis(1, at=1:5, labels=list, cex.axis=1)



# load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Kiuic_climate_change_ssp434_MOHC_1899_2100_site_Kiuic_',defo_perc,'.RData', sep = ""))  
# states_all_Kiuic_site_Kiuic_met_434 <- states_all
# dat <-states_all_Kiuic_site_Kiuic_met_434$wood_gCm2
# for (i in 1:(length(dat))) {
#   rate_434 = (dat[i] - dat[i-1])/ dat[i-1] * 100 }
# load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Kiuic_climate_change_ssp370_MOHC_1899_2100_site_Kiuic_',defo_perc,'.RData', sep = ""))    
# states_all_Kiuic_site_Kiuic_met_370 <- states_all
# dat <- states_all_Kiuic_site_Kiuic_met_370$wood_gCm2
# for (i in 1:(length(dat))) {
#   rate_370 = (dat[i] - dat[i-1])/ dat[i-1] * 100 }
# load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Kiuic_climate_change_ssp245_MOHC_1899_2100_site_Kiuic_',defo_perc,'.RData', sep = ""))    
# states_all_Kiuic_site_Kiuic_met_245 <- states_all
# dat <- states_all_Kiuic_site_Kiuic_met_245$wood_gCm2
# for (i in 1:(length(dat))) {
#   rate_245 = (dat[i] - dat[i-1])/ dat[i-1] * 100 }
# load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Kiuic_climate_change_ssp126_MOHC_1899_2100_site_Kiuic_',defo_perc,'.RData', sep = ""))
# states_all_Kiuic_site_Kiuic_met_126 <- states_all
# dat <- states_all_Kiuic_site_Kiuic_met_126$wood_gCm2
# for (i in 1:(length(dat))) {
#   rate_126 = (dat[i] - dat[i-1])/ dat[i-1] * 100 }
# 
# 
# 
# get_rate <- function(dat) {
#   for (i in 1:(length(dat))) {
#     rate = (dat[i] - dat[i-1])/ dat[i-1] * 100
#   }
# }
# 
# 
# 
# 
# 
