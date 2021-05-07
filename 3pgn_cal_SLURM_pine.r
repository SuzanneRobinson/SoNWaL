library(ProfoundData)
library(tidyverse)
library(lubridate)
library(RNetCDF)
library(zoo)
library(fr3PGDN)
library(BayesianTools)
library(dplyr)
library(coda)
library(miscTools)

setwd("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PROFOUND_data\\ProfoundData")
setDB(db_name ="ProfoundData.sqlite")
timeStep="weekly"

clm_df_full<-getClimDat(timeStep)
sitka<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
sites    <- c("peitz","soro" ,"hyytiala","solling_beech","solling_spruce")
fName=paste0("outx_",stringr::str_sub(Sys.time(), 0, -10),stringr::str_sub(Sys.time(), 15, -4),stringr::str_sub(Sys.time(), 18, -1),".RDS")

##hyytiala
site = 3

clm_df_pine<-getClmPine(timeStep)

#get parameters
pine<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
pine$weather<-clm_df_pine

#designate management cycle
presc<-data.frame(cycle=c(1,1,1),t=c(15,25,35),pNr=c(0.4,0.3,0.075),thinWl=c(0.4,0.3,0.075),
                  thinWsbr=c(0.4,0.3,0.075),thinWr=c(0.4,0.3,0.075),t.nsprouts=c(1,1,1))
pine$presc<-presc
###Get observed flux data
fluxDat <- getData(site=sites[site],dataset="FLUX")
observedPine<-observedValsPine(timeStep = timeStep,fluxDat = fluxDat)[[1]]
devPine<-observedValsPine(timeStep = timeStep,fluxDat = fluxDat)[[2]]

nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")


##Set priors
priors<-createPriors_pine(pine)
startYear = 1996
endYear = 2014

likelihoodFunc<-if(timeStep=="monthly") pineLL else pineLL_weekly

iters=250000
#Initiate bayesian setup
BS3PGDN <- createBayesianSetup(likelihood = likelihoodFunc, prior = Uprior, names = nm, parallel = 7, catchDuplicates = F )
settings = list(
  iterations = iters,
  startValue = 7, # internal chain number - dont use these chains for convergence testing 
  nrChains = 1, # Number of chains
  pSnooker = 0.5,
  burnin = round(iters/100*10), #10% burnin
  parallel = T,
  message = TRUE)

#run calibration with all parameters and priors based on initial hydro model runs
out2 <- runMCMC(bayesianSetup = BS3PGDN, sampler = "DEzs", settings = settings)


#Save output
saveRDS(out2,file=paste0(timeStep,"_pine_",fName))

