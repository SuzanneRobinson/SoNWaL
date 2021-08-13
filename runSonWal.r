

## Load necessary packages
library(fr3PGDN)
library(tidyverse)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
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

  
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn")
  
##update with some example calibrated parameters
exampParams<-read.csv("data//exampParams.csv")
sitka[nm]<-exampParams[nm]

output<-do.call(fr3PGDN,sitka)%>%
  filter(Year>2014)

quickPlotWeekly(flxdata_daily,output)


##quick function for plotting weekly simulated data against observed, does not include uncertainty as this requires mcmc posteriors
quickPlotWeekly<-function(flxdata_daily,output){
  
  output<-filter(output,Year>2014)
  
  #conversion factor to gCcm-2 for weekly data
  cf=7.142857
  
  flxdata_daily$week<-week(flxdata_daily$timestamp)
  flxdata_weekly<-flxdata_daily%>%
    filter(year>2014)%>%
    group_by(year,week)%>%
    summarise(gpp=mean(gpp),npp=mean(npp),nee=mean(nee),reco=mean(reco),rs=mean(rs),et=sum(et),swc=mean(swc),timestamp=median(timestamp))%>%
    add_column(simGPP=output$GPP*cf,simNPP=output$NPP*cf,simNEE=output$NEE*cf,simRS=output$Rs*cf,simET=output$EvapTransp,simSWC=output$volSWC_rz)

  gpp1<-ggplot()+theme_bw()+
    geom_line(data=flxdata_weekly,aes(x=timestamp,y=simGPP),colour="purple",size=1)+
    geom_point(data=flxdata_weekly,aes(x=timestamp, y=gpp),colour="black",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("GPP [gC"," ",cm^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  gpp2<-ggplot()+theme_bw()+
    geom_line(data=flxdata_weekly,aes(x=timestamp,y=simNPP),colour="purple",size=1)+
    geom_point(data=flxdata_weekly,aes(x=timestamp, y=npp),colour="black",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("NPP [gC"," ",cm^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  gpp3<-ggplot()+theme_bw()+
    geom_line(data=flxdata_weekly,aes(x=timestamp,y=simNEE),colour="purple",size=1)+
    geom_point(data=flxdata_weekly,aes(x=timestamp, y=nee),colour="black",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("NEE [gC"," ",cm^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  gpp4<-ggplot()+theme_bw()+
    geom_line(data=flxdata_weekly,aes(x=timestamp,y=simET),colour="purple",size=1)+
    geom_point(data=flxdata_weekly,aes(x=timestamp, y=et),colour="black",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("GPP [gC"," ",cm^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  gpp5<-ggplot()+theme_bw()+
    geom_line(data=flxdata_weekly,aes(x=timestamp,y=simSWC),colour="purple",size=1)+
    geom_point(data=flxdata_weekly,aes(x=timestamp, y=swc),colour="black",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("swc %",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  gpp6<-ggplot()+theme_bw()+
    geom_line(data=flxdata_weekly,aes(x=timestamp,y=simRS),colour="purple",size=1)+
    geom_point(data=flxdata_weekly,aes(x=timestamp, y=rs),colour="black",size=2)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("swc %",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  ggarrange(gpp1,gpp2,gpp3,gpp4,gpp5,gpp6)
  
}


