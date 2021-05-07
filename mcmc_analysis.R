##Analysis of MCMC
library(fr3PGDN,quietly=TRUE)
library(tidyverse,quietly=TRUE)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
library(MCMCvis)

## Met data
timeStep="monthly"
clm_df_full<-getClimDat(timeStep)
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
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
                E_S1 =1, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                K_drain=0.16,
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
)
#######################################################

nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")
##Set priors
priors<-createPriors_sitka(sitka=sitka)
pMaxima<-priors[[1]]
pMinima<-priors[[2]]
pMaxima[[30]]<-0.1
Uprior <- createPrior(lower = pMinima, upper = pMaxima)

out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_1_T.RDS")

codM<-getSample(out, start = 100, coda = TRUE, thin = 10)
priorSamp<-Uprior$sampler(35000)
MCMCtrace(codM,wd="C:\\Users\\aaron.morris", post_zm=F,iter=5000)
codM<-as.data.frame(mergeChains(out$chain))
names(codM)<-nm
codM<-colMedians(as.data.frame(codM))
sitka[nm]<-codM[nm]
sitka$waterBalanceSubMods<-T
sitkaBICcomp(data,startY=2015,endY=2018,pNum=36,sitka)



output<-do.call(fr3PGDN,sitka)
tail(output$GPP)
results<-plotResults(output,ShortTS=F)
ggarrange(results[[2]],results[[6]],results[[5]],results[[12]])


pine[nm]<-codM[nm]

pineBICcomp(data,startY=1996,endY=2014,pNum=36,pine)
