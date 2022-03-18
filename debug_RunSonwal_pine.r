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
library(ggpubr)

##read in the arguments listed at the command line
#args=(commandArgs(TRUE))
args<-c("weekly","weekly_1_","1","T")

print(args)
timeStep=args[1]
print(timeStep)
chainNum=args[3]
chainID=args[2]

setDB(db_name ="C:/Users/aaron.morris/OneDrive - Forest Research/Documents/Projects/PRAFOR/models/PROFOUND_data/ProfoundData/ProfoundData.sqlite")
flxdata_daily <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))

clm_df_full<-getClimDat(timeStep)
sitka<-getParms(waterBalanceSubMods=ifelse(args[4]=="T",TRUE,FALSE),timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
sites    <- c("peitz","soro" ,"hyytiala","solling_beech","solling_spruce")

##hyytiala
site = 3
clm_df_pine<-getClmPine(timeStep)

nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")


#Plot results Pine
outP<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_pine_2_T.RDS")
pine<-getParmsPine(waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
pine$weather<-clm_df_pine
#Set management
presc<-data.frame(cycle=c(1,1,1),t=c(15,25,35),pNr=c(0.4,0.3,0.075),thinWl=c(0.4,0.3,0.075),
                  thinWsbr=c(0.4,0.3,0.075),thinWr=c(0.4,0.3,0.075),t.nsprouts=c(1,1,1))
pine$presc<-presc
resultsPine2<-plotResultsPine(outP)

pineDiags<-diagPlotsPine(outP=outP,nm = nm )


