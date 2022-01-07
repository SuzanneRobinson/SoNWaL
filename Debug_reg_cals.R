## Load necessary packages
library(fr3PGDN,quietly=TRUE)
library(BayesianTools,quietly=TRUE)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)
library(Rmpi)
library(lubridate)
library(doMPI)
library(parallel)
library(matrixStats)
#read in arguments from batch file
args=(commandArgs(TRUE))
args<-c("weekly","weekly_1_","1","T")

print(args)
timeStep=args[1]
chainID=args[2]
chainNum=args[3]
#print arguments to log for easy reference
print(timeStep)
print(chainID)
print("sitka")

#Time step to run SonWal with
timeStep<-"weekly"

#Directory where climate data is stored (default is "data" in SonWal folder)
climDir<-("Data/")

clm_df_reg<-readRDS("Data/regionalClmDat.RDS")
clm_df_reg$RH<-calcRH(Tmean=clm_df_reg$Tmean,Tref=273.16,p=clm_df_reg$psurf_pa,q=clm_df_reg$specHumid)
clm_df_reg$VPD <-((( (0.61078 * exp(17.269 * (clm_df_reg$Tmean)/(237.3 + clm_df_reg$Tmean)))*(1-(clm_df_reg$RH)/100))))
clm_df_reg<-data.frame(site=clm_df_reg$siteName, Year=clm_df_reg$year,week=week(clm_df_reg$date),Month=clm_df_reg$month,Tmax=clm_df_reg$Tmax,Tmin=clm_df_reg$Tmin,
                       Tmean=clm_df_reg$Tmean,Rain=clm_df_reg$precip_mm,SolarRad=clm_df_reg$solarRad_MJ,FrostDays=ifelse(clm_df_reg$Tmin<=0,1,0),MonthIrrig=0,VPD=clm_df_reg$VPD,RH=clm_df_reg$RH,rainDays=ifelse(clm_df_reg$precip_mm>=1,1,0),
                       wp=clm_df_reg$wp,fc=clm_df_reg$fc,wp=clm_df_reg$wp,soilDepth=clm_df_reg$soilDepth,soilCond=clm_df_reg$soilCond,soilTex=clm_df_reg$soilTex)

clm_df_reg<-clm_df_reg%>%group_by(site,Year,week)%>%summarise(Year=median(Year),week=median(week),Month=median(Month),
                                                              Tmax=max(Tmax),Tmin=min(Tmin),Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad),
                                                              FrostDays=sum(FrostDays),MonthIrrig=0,VPD=mean(VPD),RH=mean(RH),rainDays=sum(rainDays),
                                                              wp=median(wp),fc=median(fc),wp=median(wp),soilDepth=median(soilDepth),soilCond=median(soilCond),soilTex=median(soilTex))

## read in and format climate data
#clm_df_full<-data.frame(getClimDatX(timeStep,climDir))%>%
#  filter(Year<2019)

bMarkDat<-read.csv("Data/bMarkCalibrationData.csv")
bMarkDat2<-bMarkDat[,c(2,5:12)]
bMarkDat2<-unique(bMarkDat2)
bMarkDat<-bMarkDat[,c(2,13:17)]

clm_df_reg<-merge(clm_df_reg,bMarkDat2,by.x="site",by.y="SiteIdentification",all=T)
clm_df_reg<-merge(clm_df_reg,bMarkDat,by.x=c("site","Year"),by.y=c("SiteIdentification","Year"),all=T)

clm_df_full<-clm_df_reg


priorVals<-createPriors_sitka_rg(sitka)

  
  runModReg<-function(g,dg=T){
    
   plotRes<-function(gx,sY,eY,dg){
        sitka[.GlobalEnv$nm]<-gx[1:length(nm)]
        sitka$weather<-filter(sitka$weather,site==gx$site)
        
        if(sitka$weather$soilDepth[1]==1){ 
          sitka$V_nr<- 0.25
          sitka$maxRootDepth<- 0.25
        }
        
        if(sitka$weather$soilDepth[1]==2){
          sitka$V_nr<- 0.5
          sitka$maxRootDepth<- 0.5
        }
        if(sitka$weather$soilDepth[1]==3){
          sitka$V_nr<- 0.8
          sitka$maxRootDepth<- 0.8
        }
        if(sitka$weather$soilDepth[1]==4){
          sitka$V_nr<- 1
          sitka$maxRootDepth<- 1
        }
        if(sitka$weather$soilDepth[1]==5){
          sitka$V_nr<- 1.5
          sitka$maxRootDepth<- 1.5
        }
        
        if(sitka$weather$soilTex[1]!=0){
          sitka$wiltPoint<-sitka$weather$wp[1]
          sitka$fieldCap<-sitka$weather$fc[1]
          sitka$satPoint<-sitka$weather$sp[1]
          sitka$K_s<-sitka$weather$soilCond[1]
          sitka$K_drain<-sitka$weather$soilCond[1]
        }
        if(sitka$weather$soilTex[1]==1) {
          sitka$E_S1<-0.05
          sitka$E_S2<-0.3
          sitka$SWpower0<-8
          sitka$SWconst0<-0.65
        }
        
        if(sitka$weather$soilTex[1]==0) {
          sitka$E_S1<-0.05
          sitka$E_S2<-0.3
          sitka$SWpower0<-8
          sitka$SWconst0<-0.65
        }
        
        if(is.na(sitka$weather$soilTex[1])==T)  {
          sitka$E_S1<-0.05
          sitka$E_S2<-0.3
          sitka$SWpower0<-8
          sitka$SWconst0<-0.65
        }
        
        if(sitka$weather$soilTex[1]==4)  {
          sitka$E_S1<-0.1
          sitka$E_S2<-0.3
          sitka$SWpower0<-5
          sitka$SWconst0<-0.5
        }
        
        
        if(sitka$weather$soilTex[1]==9) {
          sitka$E_S1<-0.3
          sitka$E_S2<-0.6
          sitka$SWpower0<-5
          sitka$SWconst0<-0.5
        }
        
        output<-   do.call(fr3PGDN,sitka)
        if(dg==T) return(output%>%filter(Year>=sY&Year<=eY)%>%group_by(Year)%>%summarise(dg=mean(dg))) else    
        return(output%>%filter(Year>=sY&Year<=eY)%>%group_by(Year,Month)%>%summarise(dg=mean(GPP)))

   }
   
   observed<-sitka$weather%>%filter(site==g$site[1])%>%filter(is.na(mean_dbh_cm)==F)%>%group_by(Year)%>%summarise(dbh=median(mean_dbh_cm),dbhSD=median(dbhSD_cm))
   observed$dbhSD<-ifelse(observed$dbhSD==0,0.0001,observed$dbhSD)
   sY=min(observed$Year)
   eY=max(observed$Year)
   
   gx<-split(g,seq(nrow(g)))
   ff<-mapply(plotRes,gx,MoreArgs = list(sY=sY,eY=eY,dg=dg),SIMPLIFY = F)
   
   getIntv<-function(paramName,modLst){
     simDat =   do.call(cbind,lapply(modLst, "[",  paramName))
     res<-as.data.frame(rowQuantiles(as.matrix(simDat), probs = c(0.11, 0.5, 0.89)))
     return(res)
   }
 
        #ff<-lapply(ff, function(x) filter(x, Year >=sY & Year <=eY))
        
        paramName<-("dg")
        intvsS<-mapply(getIntv,paramName,MoreArgs = list(modLst=ff))
        
   if(dg==T){
   observed$high  <- intvsS[,1]$`89%` + observed$dbhSD
   observed$low  <- intvsS[,1]$`11%` - observed$dbhSD
   observed$med  <- intvsS[,1]$`50%`# - 2 * 0.3
   
  res<- ggplot(data=observed)+
    geom_point(aes(y=dbh,x=Year))+
     geom_line(aes(y=med,x=Year))+
   geom_ribbon(aes(ymax=high,ymin=low,x=Year),alpha=0.3)
   } else {
     
     simul<-data.frame(Date=as.Date(paste0(ff$`1`$Year,"-",ff$`1`$Month,"-01")),Year=ff$`1`$Year,med=intvsS[,1]$`50%`*7.14,high=intvsS[,1]$`89%`*7.14,low=intvsS[,1]$`11%`*7.14)
     res<-ggplot(data=simul)+
       geom_line(aes(y=med,x=Date))+
       geom_ribbon(aes(ymax=high,ymin=low,x=Date),alpha=0.3)+
       ylab("GPP")
   }
        return(res)
  }
  
  
  out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_2_T.RDS")
  nmc<-nrow(out$chain[[1]])
  outSample   <- as.data.frame(getSample(out,start=round(nmc/1.1),thin=1,numSamples = 10))
  sitka<-getParms(
    waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
  sitka$weather<-clm_df_reg
  siteLst<-(unique(sitka$weather$site))
  outSampleX<-outSample[rep(seq_len(nrow(outSample)), length(siteLst)), ]
  outSampleX$site<-rep(siteLst,each=nrow(outSample))
  
  out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_3_final_T.RDS")
  codM<-mergeChains(out$chain)
  codM<-miscTools::colMedians(as.data.frame(codM))
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  names(codM)<-nm
  #get parameter values and observed data
  sitka<-getParms(
    waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
  sitka$weather<-clm_df_reg
  sitka[nm]<-codM[nm]
  
  nm<-c("sigma_zR","shared_area",
        "Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","startN","startC")
 # nm<-c("sigma_zR","shared_area",
 #       "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
 #       "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
 #       "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","Qa","Qb","MaxIntcptn","k","startN","startC")
  

  
  g <- split(outSampleX,outSampleX$site)

  
    g2<-mapply( runModReg,g[c(1:6)],MoreArgs = list(dg=T),SIMPLIFY = F)

