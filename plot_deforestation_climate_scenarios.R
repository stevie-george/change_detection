###
## Make a list of deforestation percent
## 
par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1))
list = c("25", "50", "75", "90", "99")


## load data for specific scenario
load('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Site_run_tests_no_def/FCP_climate_change_ssp585_MOHC_1899_2100_test2_FCP.RData')    
states_all_FCP_585 <- states_all

## plot no deforestation as base
plot(apply(states_all_FCP_585$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,25000), lwd=3, type="l", col="red", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5, lty = "dashed")
##
# Begin looping through model experiment combinations
# Experiment option
for (i in 1:(length(list))) {
  load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp585_MOHC_1899_2100_site_FCP_',list[i],'.RData', sep = ""))    
  states_all_FCP_site_Kiuic_met_585 <- states_all
  # Load first experiment information to extract needed timing information
  
  lines(apply(states_all_FCP_site_Kiuic_met_585$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,29000), lwd=3, type="l", col="red", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5)
  
} # plot

### 434

load('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Site_run_tests_no_def/FCP_climate_change_ssp434_MOHC_1899_2100_test2_FCP.RData')    
states_all_FCP_434 <- states_all
plot(apply(states_all_FCP_434$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,25000), lwd=3, type="l", col="orange", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5, lty = "dashed")
##
# Begin looping through model experiment combinations
# Experiment option
for (i in 1:(length(list))) {
  load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp434_MOHC_1899_2100_site_FCP_',list[i],'.RData', sep = ""))    
  states_all_FCP_site_Kiuic_met_434 <- states_all
  # Load first experiment information to extract needed timing information
  
  lines(apply(states_all_FCP_site_Kiuic_met_434$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,29000), lwd=3, type="l", col="orange", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5)
  
} # export



list = c("25", "50", "75", "90", "99")
# filter for the selected site
#avail_scenarios = avail_scenarios[grepl(site, avail_scenarios)]
load('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Site_run_tests_no_def/FCP_climate_change_ssp370_MOHC_1899_2100_test2_FCP.RData')    
states_all_FCP_370 <- states_all
plot(apply(states_all_FCP_370$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,25000), lwd=3, type="l", col="yellow", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5, lty = "dashed")
##
# Begin looping through model experiment combinations
# Experiment option
for (i in 1:(length(list))) {
  load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp370_MOHC_1899_2100_site_FCP_',list[i],'.RData', sep = ""))    
  states_all_FCP_site_Kiuic_met_370 <- states_all
  # Load first experiment information to extract needed timing information
  
  lines(apply(states_all_FCP_site_Kiuic_met_370$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,29000), lwd=3, type="l", col="yellow", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5)
  
} # export



list = c("25", "50", "75", "90", "99")
# filter for the selected site
#avail_scenarios = avail_scenarios[grepl(site, avail_scenarios)]
load('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Site_run_tests_no_def/FCP_climate_change_ssp245_MOHC_1899_2100_test2_FCP.RData')    
states_all_FCP_245 <- states_all
plot(apply(states_all_FCP_245$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,25000), lwd=3, type="l", col="green", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5, lty = "dashed")
##
# Begin looping through model experiment combinations
# Experiment option
for (i in 1:(length(list))) {
  load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp245_MOHC_1899_2100_site_FCP_',list[i],'.RData', sep = ""))    
  states_all_FCP_site_Kiuic_met_245 <- states_all
  # Load first experiment information to extract needed timing information
  
  lines(apply(states_all_FCP_site_Kiuic_met_245$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,29000), lwd=3, type="l", col="green", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5)
  
} # export

list = c("25", "50", "75", "90", "99")
# filter for the selected site
#avail_scenarios = avail_scenarios[grepl(site, avail_scenarios)]
load('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/Site_run_tests_no_def/FCP_climate_change_ssp126_MOHC_1899_2100_test2_FCP.RData')    
states_all_FCP_126 <- states_all
plot(apply(states_all_FCP_126$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,25000), lwd=3, type="l", col="blue", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5, lty = "dashed")
##
# Begin looping through model experiment combinations
# Experiment option
for (i in 1:(length(list))) {
  load(paste('/home/stevie/CARDAMOM/regrowth/climate_change_scenario_running/outputs/FCP_climate_change_ssp126_MOHC_1899_2100_site_FCP_',list[i],'.RData', sep = ""))    
  states_all_FCP_site_Kiuic_met_126 <- states_all
  # Load first experiment information to extract needed timing information
  
  lines(apply(states_all_FCP_site_Kiuic_met_126$wood_gCm2,2,quantile,prob=c(0.50)), ylim=c(0,29000), lwd=3, type="l", col="blue", ylab="Wood (gC/m2)", xlab="Time (Months)", cex=1.5, cex.axis=1.5, cex.lab=1.5)
  
} # export

