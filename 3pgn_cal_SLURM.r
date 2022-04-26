## Load necessary packages
library(SoNWaL)
library(BayesianTools)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)

##choose timestep size
timeStep<-"weekly"

#get climate data
clm_df_full<-getClimDat(timeStep)
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  flxdata <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  flxdata_daily <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  }else
{
  flxdata <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  flxdata_daily <- read.csv("/home/users/aaronm7/3pgData/harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  }


#get parameter values 
sitka<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
flxdata<-if(timeStep=="monthly")flxdata else flxdata_daily

observed<-observedVals(timeStep = timeStep,data=flxdata)[[1]]
dev<-observedVals(timeStep = timeStep,data=flxdata)[[2]]

nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")
##Set priors
priorVals<-createPriors_sitka(sitka=sitka)
startYear = 2015
endYear = 2018

observed<-if(sitka$waterBalanceSubMods==FALSE&&timeStep=="monthly")observed[-c(295:342)] else observed
dev<-if(sitka$waterBalanceSubMods==FALSE&&timeStep=="monthly")dev[-c(295:342)] else dev
observed<-if(sitka$waterBalanceSubMods==FALSE&&timeStep=="weekly")observed[-c(1279:1490)] else observed
dev<-if(sitka$waterBalanceSubMods==FALSE&&timeStep=="weekly")dev[-c(1279:1490)] else dev

if(sitka$waterBalanceSubMods==F&&timeStep=="monthly") likelihoodFunc<-NLL_noHYD
if(sitka$waterBalanceSubMods==F&&timeStep!="monthly") likelihoodFunc<-NLL_noHYDWeekly
if(sitka$waterBalanceSubMods==T&&timeStep!="monthly") likelihoodFunc<-NLL_weekly
if(sitka$waterBalanceSubMods==T&&timeStep=="monthly") likelihoodFunc<-NLL




for (i in c(1:15)){
  iters=100000

#Initiate bayesian setup
#BS3PGDN <- createBayesianSetup(likelihood = likelihoodFunc, prior = Uprior, names = nm, parallel = 8, catchDuplicates = F )
settings = list(
  iterations = iters,
  startValue = 7, # internal chain number - dont use these chains for convergence testing 
  nrChains = 1, # Number of chains
  pSnooker = 0.5,
  burnin = round(iters/100*10), #10% burnin
  parallel = T,
  message = TRUE)

#out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_3_T.RDS")

if(file.exists(paste0("/home/users/aaronm7/",timeStep,"_",chainNum,"_",args[4],".RDS"))==TRUE){
  print("previous file exists, restarting chain")
  out<-readRDS(paste0("/home/users/aaronm7/",timeStep,"_",chainNum,"_",args[4],".RDS"))
}

BS3PGDN <- createBayesianSetup(likelihood = likelihoodFunc, prior = priorVals, names = nm, parallel = 8, catchDuplicates = F )

out<-if(file.exists(paste0("/home/users/aaronm7/",timeStep,"_",chainNum,"_",args[4],".RDS"))==TRUE) runMCMC(bayesianSetup =out, sampler = "DEzs", settings = settings) 
else runMCMC(bayesianSetup =BS3PGDN, sampler = "DEzs", settings = settings)
#Save output
saveRDS(out,file=paste0(timeStep,"_",chainNum,"_",args[4],".RDS"))

stopParallel(BS3PGDN)
rm(BS3PGDN)


}
  
  


#
#priorSDfunc<-function(Uprior){
#  prVals<-Uprior$sampler(10000)
#  prVals[,7]<-exp(prVals[,7])
#  prVals[,8]<-exp(prVals[,8])
#  prValsX<-as.data.frame(round(HPDinterval(as.mcmc(prVals)),3))
#  prValsX$med<-colMedians(prVals)
#  prValsX$names<-nm
#}
#


