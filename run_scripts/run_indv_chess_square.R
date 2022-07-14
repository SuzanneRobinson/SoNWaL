library(tidyr)
library(tidyverse)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(dplyr)
library(future)
library(ncdf4)
library(ggplot2)
library(httr)
library(furrr)
library(viridis)
library(SoNWaL)
library(tibble)
library(miscTools)
library(parallel)
library(sf)
library(rgdal)
library(raster)

timeStep="weekly"


# define location of data sets
clm_df_reg<-readRDS("data/regionalClmDat.RDS")

# list of soil data locations, _fao files are the map files from https://esdac.jrc.ec.europa.eu/content/maps-indicators-soil-hydraulic-properties-europe
# soildataNZero2 is the data from Astley derived from BGS data
soil_dat<-c(
  "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\wp_fao.tif",
  "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\fc_fao.tif",
  "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\ks_fao_octop.tif",
  "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\ths_fao_octop.tif",
  "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\spatialData\\soildataNZero2.RDS")

simDatLoc<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\misc_data\\AH_Harwood_CHESS_dat\\clm_rcp45_01_wide.RDS"

#simDat<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\misc_data\\AH_Harwood_CHESS_dat\\clm_hist_wide.RDS")

#paramsFile<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\reg_cals\\weekly_24_T.RDS")
paramsFile<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\final_results\\sitka\\HD\\weekly\\mcmcOut\\weekly_3_final_T.RDS")

#exampParams<-read.csv("exampParams.csv")
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")


# read in the mcmc results file to sample from the posterior
#paramsFile<-conc_chains(paramsFile,  clm_df_reg=clm_df_reg)
paramsFile<-mergeChains(paramsFile$chain)
paramsFile<-as.data.frame(paramsFile)
names(paramsFile)<-nm
codM<-colMedians(paramsFile)
#codM<-as.data.frame(t(miscTools::colMedians(as.data.frame(paramsFile) )))
#codM<-as.data.frame(t(codM)) 
#param_scaler<-createPriors_sitka(sitka=getParms(E_S1=1,E_S2 = 1))[[2]]
#codM<-codM*param_scaler
scape<-T

baseParms<-getParms(weather=as.data.frame(clm_df_reg),
                    waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                    timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
)




baseParms[nm]<-codM[nm]

baseParms<-getIndvClm(scape, "Harwood", soil_dat, simDatLoc, baseParms)
#weather<-baseParms$weather
#weather$Year<-weather$Year-30
#weather<-filter(weather,Year<1961)
#baseParms$weather<-rbind(weather,baseParms$weather)


  
output<-do.call(SoNWaL,baseParms)
summary(output$GPP)
plot(output$LAI)

ff<-output%>%group_by(Year, Month) %>% summarise (NPP_gC_m2_day=mean(NPP*7.14
                                                                      ), kl=mean(kl), littFall_gC_m2_day = mean(difLitter*7.14, na.rm = T)) %>% filter(Year>2014 & Year < 2019)

plot(ff$littFall_gC_m2_day)

saveRDS(ff, "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\misc_data\\NPP_sims_for_jon\\NPP_harwood_origMeteo.RDS")




paramsFile<-as.data.frame(paramsFile)
paramsFile$wiltPoint=simDat$wp
paramsFile$fieldCap=simDat$fc
paramsFile$satPoint=simDat$sp
paramsFile$K_s=simDat$soil_cond
paramsFile$startC = simDat$SOC
paramsFile$V_nr=simDat$soil_depth
paramsFile$maxRootDepth=simDat$soil_depth
paramsFile$SWpower0=simDat$no
paramsFile$SWconst0=simDat$co

results<-plotResultsNewMonthly(output,ShortTS=T,out=outSample,numSamps = 25)
ggarrange(results[[1]],results[[8]],results[[2]],results[[3]],results[[6]],results[[4]],nrow=3,ncol=2)
ggarrange(results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]],nrow=3,ncol=2)

