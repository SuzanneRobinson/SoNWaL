

######################
## START THE SCRIPT ##
######################

## Load necessary packages
library(fr3PGDN)
library("tidyverse")
library("lubridate")
library("coda")
library("BayesianTools")
library("miscTools")
library("ggpubr")
## Years of data to use for calibration
startYear = 2015
endYear = 2018
#install.packages("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\fr3PGDN_2.0.tar.gz", repos = NULL, type="source")
timeStep<-"weekly"

#Directory where climate data is stored (default is "data" in SonWal folder)
climDir<-("Data/")

## read in and format climate data
clm_df_full<-data.frame(getClimDatX(timeStep,climDir))%>%
  filter(Year<2019)


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
                #K_s=0.1, #Soil conductivity
                shared_area=3, #shared area of rooting and non-rooting zone
                V_nr=4, #Volume of non-rooting zone
                maxRootDepth=2.5,
                sigma_zR =0.3, #area/depth explored by 1kg of root biomass
                SWC_nr=10, #SWC of non-rooting zone at time 0
                E_S1 =4.5, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                #K_drain=0.16,
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################
#names of parameters to fit

#exampParams<-read.csv("exampParams.csv")
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC","kF")



out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_8_T.RDS")
#clm1<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\spatialChunk_28_22.RDS")

#out<-getSample(out,start=12000,thin=5,numSamples=500)
#codM<-out$chain[[2]][c(1:5000),]
codM<-as.data.frame(out$chain[[3]])
codM<-mergeChains(out$chain)

codM<-miscTools::colMedians(as.data.frame(codM))
codM<-tail(as.data.frame(codM),1)
names(codM)<-nm


#priorSamp<-priorVals$sampler(35000)
#MCMCtrace(getSample(out,coda = T,thin=1,start=100),wd="C:\\Users\\aaron.morris", post_zm=F,iter=10000,priors = priorSamp)

sitka<-getParms(waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
#sitka$E_S1<-2
#sitka$weather<-clm_df_pine
param_scaler<-createPriors_sitka(sitka=getParms(E_S1=1,E_S2 = 1))[[2]]

sitka[nm]<-codM[nm]*param_scaler

output<-do.call(fr3PGDN,sitka)
ff<-filter(output,Year>2014&Year<2100)
plot(ff$Rs*1.6)
plot(output$Nav)
plot(output$LAI)


plot(output$N)
plot(ff$EvapTransp)
tail(output$GPP)
#ff<-filter(output,Year>2014)
plot(ff$fSW)
plot(ff$volSWC_rz)
plot(ff$EvapTransp)
plot(output$totC)
plot(ff$GPP*7.14,col="red")


#Plot results
nmc = nrow(out$chain[[1]])
outSample   <- as.matrix(getSample(out,start=nmc/1.1,thin=1,numSamples = 100))
outSample<-as.data.frame(outSample%*%diag(param_scaler))
names(outSample)<-nm
results<-plotResultsNewMonthly(output,ShortTS=T,out=outSample,numSamps = 100)
ggarrange(results[[1]],results[[8]],results[[2]],results[[3]],results[[6]],results[[4]],nrow=3,ncol=2)
ggarrange(results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]],nrow=3,ncol=2)

diagPlots(out,flxdata_daily = flxdata_daily,nm=nm,param_scaler)


ggarrange(results[[1]],results[[2]],results[[3]],results[[5]],results[[4]],results[[6]],results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]])
fileNm<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\"
saveRDS(out,file=paste0(fileNm,timeStep,"_",iters,"_.RDS"))
ggsave(paste0(fileNm,timeStep,"_",iters,".png"),width = 400,
       height = 300,
       units = "mm")

