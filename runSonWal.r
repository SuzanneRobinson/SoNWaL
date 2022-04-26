

## Load necessary packages
library(SoNWaL)
library(tidyverse)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
library(matrixStats)
library(future)
library(furrr)
library(parallel)
## Years of data to use for calibration
startYear = 2015
endYear = 2018

#Time step to run SonWal with
timeStep<-"weekly"

#Directory where climate data is stored (default is "data" in SonWal folder)
climDir<-("data\\")

## read in and format climate data
clm_df_full<-data.frame(getClimDatX("weekly",climDir))%>%
  filter(Year<2019)

## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
  flxdata_daily <- read.csv("data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
#################################
## load in default parameters  ##
#################################
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################

#Names of fitted parameters  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
        "SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  
  
##update with some example calibrated parameters (parameters are sample from full MCMC chain of Harwood fitting)
exampParams<-as.data.frame(readRDS("data//exampParams.RDS"))
exParms<-miscTools::colMedians(as.data.frame(exampParams))
names(exParms)<-nm
sitka[nm]<-exParms[nm]
#run SoNWal
output<-do.call(SoNWaL,sitka)

##plotting QUICK PLOT - DOES NOT INCLUDE UNCERTAINTY! BUT QUICK :) -grouping aggregates the data, can be either "week" or "month"
quickPlot(flxdata_daily,output,grouping="month")

#full plots - much slower but gives credible intervals and uncertainty of observed data - num samps is how many samples from posterior to use (>=500 ideal but 50-100 will give a pretty solid output for a quick checking)
results<-plotResultsNewMonthly(output,ShortTS=T,out=exampParams,numSamps = 50)
ggarrange(results[[1]],results[[2]],results[[8]],results[[3]],results[[5]],results[[4]],ncol=2,nrow=3)
ggarrange(results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]],ncol=2,nrow=3)
