

######################
## START THE SCRIPT ##
######################

## Load necessary packages
library(fr3PGDN,quietly=TRUE)
install.packages("tidyverse")
install.packages("lubridate")
install.packages("coda")
install.packages("BayesianTools")
install.packages("miscTools")
install.packages("ggpubr")
## Years of data to use for calibration
startYear = 2015
endYear = 2018
#install.packages("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\fr3PGDN_2.0.tar.gz", repos = NULL, type="source")
timeStep<-"monthly"

## Met data
clm_df_full<-data.frame(getClimDat(timeStep))


clm_df_full<-clm_df_full%>%filter(Year<2019)
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  flxdata_daily <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
  }else
{
  data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
}
###########################
## Initialise Parameters ##
###########################
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                wiltPoint = 0.1, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
                fieldCap =0.24,#Field capacity
                satPoint = 0.3, #field saturation point
                K_s=0.1, #Soil conductivity
                shared_area=3, #shared area of rooting and non-rooting zone
                V_nr=4, #Volume of non-rooting zone
                maxRootDepth=2.5,
                sigma_zR =0.3, #area/depth explored by 1kg of root biomass
                SWC_nr=10, #SWC of non-rooting zone at time 0
                E_S1 =4.5, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                K_drain=0.16,
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################
#names of parameters to fit
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")

#out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_2_T.RDS")
codM<-as.data.frame(mergeChains(out$chain))
names(codM)<-nm
codM<-tail(as.data.frame(codM),1)

#priorSamp<-priorVals$sampler(35000)
#MCMCtrace(getSample(out,coda = T,thin=10,start=10000),wd="C:\\Users\\aaron.morris", post_zm=F,iter=5000,priors = priorSamp)

sitka<-getParms(waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
#sitka$E_S1<-2
#sitka$weather<-clm_df_pine
sitka[nm]<-codM[nm]
#sitka$SWpower0<-round(sitka$SWpower0)
#sitka$E_S2<-2
#sitka$E_S1<-2
#sitka$weather[sitka$weather$Year==2015,"Rain"][1]<-filter(clm_df_full,Year==2014)$Rain[1]

output<-do.call(fr3PGDN,sitka)
tail(output$GPP)
ff<-filter(output,Year>2014)
plot(ff$fSW)
plot(ff$volSWC_rz)
results<-plotResults(output,ShortTS=T,out=out)
#results
ggarrange(results[[1]],results[[2]],results[[3]],results[[5]],results[[4]],results[[6]])

ggarrange(results[[9]],results[[10]],results[[11]],results[[12]],results[[13]],results[[14]])

ggarrange(results[[1]],results[[2]],results[[3]],results[[5]],results[[4]],results[[6]],results[[9]],results[[10]],results[[11]],results[[12]],results[[13]],results[[14]])
fileNm<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\detailedOut\\1xVPD_03CV_01CVswc_flexPowConst\\monthly\\"
saveRDS(out,file=paste0(fileNm,timeStep,"_",iters,".RDS"))
ggsave(paste0(fileNm,timeStep,"_",iters,".png"),width = 400,
       height = 300,
       units = "mm")

#cv is ratio of sd to mean, cv * value of meanNEE
#sdmin_NEE               <- cv_EC * abs(data_NEEmean_value) #gives you the standard deviation for the averaged data points (assuming cv 20%)
#data_NEEmean_sd <- sapply( 1:length(sdmin_NEE), function(i) max(sdmin_NEE[i],0.5) )


##aggregate data by year for total annual ranfall
#annualPrecip <- weather%>%group_by(Year)%>%summarise(annual_precip=sum(Rain))
#
#
#hazprecip <- quantile(annualPrecip$annual_precip,0.1)
#
##Risk function
#calc_risk <- function(strtyr,endyr,df,hazval){
#  fldf    <- df %>% filter(Year >= strtyr & Year <= endyr)
#  lowyrs  <- fldf$Year[(fldf$annual_precip <= hazval)]
#  highyrs <- fldf$Year[!(fldf$Year %in% lowyrs)]
#  vuln    <- mean(output[weather$Year %in% highyrs,"GPP"]) - mean(output[weather$Year %in% lowyrs,"GPP"])
#  haz     <- length(lowyrs)/((endyr - strtyr)+1)
#  return(tibble("startYr"= strtyr,
#                "endYr"=endyr,
#                "vulnerability"= vuln,
#                "hazard"  = haz,
#                "risk" = vuln*haz))
#}
#
##' # Expected loss of GPP 
#inpt   <- tibble(strtyr=c(1985,1995,2005,2015),endyr=c(1988,1998,2008,2018))
#riskdf <- purrr::map2_df(inpt$strtyr,inpt$endyr,calc_risk,df=annualPrecip,hazval=hazprecip)
#knitr::kable(riskdf, digits=4, align = c(rep("c", 5)))
#