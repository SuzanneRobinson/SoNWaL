###plotting functions for FR3pgn
library(rlang)

runModel<- function(p){
  pine[.GlobalEnv$nm]<-p[nm]
  pull(do.call(SoNWaL,pine)%>%
    filter(Year>=1996)%>%
     group_by(Year,Month) %>%
   dplyr::summarise(GPP=mean(!!sym(prm))))
}

#'plotResultsPine Function to plot the model against the finnish data.
#'@return pine plots
#'@export
plotResultsPine <- function(outP,ShortTS=F){
  
  library(matrixStats)
  nmc = nrow(outP$chain[[1]])
  outSample   <- getSample(outP,start=nmc/1.2,thin=1)
  numSamples = 100# min(1000, nrow(outSample))
  codM<-miscTools::colMedians(as.data.frame(outSample))
  pine[.GlobalEnv$nm]<-codM[.GlobalEnv$nm]
  df<-do.call(SoNWaL,pine)
  
  dt=12
  df <- df[c(2:nrow(df)),]
  if(nrow(df)>600) df$week<-c(1:52) else df$week<-1
  df <- df %>% dplyr::group_by(Year)%>%mutate(cumGPP = cumsum(GPP),
                                              cumNPP = cumsum(NPP),
                                              timestamp = as.POSIXct(paste(sprintf("%02d",Year),sprintf("%02d",Month),sprintf("%02d",1),sep="-"),tz="GMT")) 
  df2<-df%>%filter(Year>=1996)
  
  flxdata<-flxdata_daily%>%
    mutate(week=week(as.Date(yday, origin = paste0(year,"-01-01"))))%>%
    mutate(month=month(as.Date(yday, origin = paste0(year,"-01-01"))))
  
  #get stand data
  standDat <- getData(site=sites[site],dataset="STAND")%>%
    filter(species_id=="pisy")
  #get soil data
  soildata <-  getData( sites[site], dataset="SOIL")
  
  df2<- (df2%>% 
           group_by(Year,Month) %>%
           dplyr::summarise(GPP=mean(GPP),EvapTransp=mean(EvapTransp),Reco=mean(Reco),Rs=mean(Rs),NPP=mean(NPP),
                            volSWC_rz=mean(volSWC_rz),NEE=mean(NEE),timestamp=median(timestamp),LAI=mean(LAI),t.proj=median(t.proj)))
  
  modif<-ifelse(nrow(df)<600,1.66,7.142857)
  nDays<-ifelse(nrow(df)<600,30,7)
  
  runModel<- function(p){
    pine<-getParmsPine(waterBalanceSubMods=T, timeStp = 52)
    presc<-data.frame(cycle=c(1,1,1),t=c(15,25,35),pNr=c(0.4,0.3,0.075),thinWl=c(0.4,0.3,0.075),
                      thinWsbr=c(0.4,0.3,0.075),thinWr=c(0.4,0.3,0.075),t.nsprouts=c(1,1,1))
    pine[.GlobalEnv$nm]<-p[nm]
    pine$presc<-presc
    res<- do.call(SoNWaL,pine)%>%
      filter(Year>=1996)%>%
      group_by(Year,Month)%>%
      dplyr::summarise(GPP=mean(GPP),NEE=mean(NEE),volSWC_rz=mean(volSWC_rz),EvapTransp=mean(EvapTransp)/nDays,Reco=mean(Reco),Rs=mean(Rs),N=mean(N),LAI=mean(LAI),dg=mean(dg),totC=mean(totC),totN=mean(totN),NPP=mean(NPP),alphaAn=mean(Reco/Rs))
    
    return(res)
  }
  
  getIntv<-function(paramName,modLst){
    
    simDat =   do.call(cbind,lapply(modLst, "[",  paramName))
    res<-as.data.frame(rowQuantiles(as.matrix(simDat), probs = c(0.11, 0.5, 0.89)))
    return(res)
  }
  
  
  outSample<-getSample(outP,start=round(nmc/1.2),numSamples = 50)
  outSample<-as.data.frame(outSample)
  outSample <- split(outSample, seq(nrow(outSample)))
  outRes<-lapply(outSample, runModel)
  
  
  paramName<-list("GPP","NEE","volSWC_rz","EvapTransp","Reco","Rs","N","LAI","dg","totC","totN","NPP","alphaAn")
  intvsS<-mapply(getIntv,paramName,MoreArgs = list(modLst=outRes))
  
  fluxDat <- getData(site=sites[site],dataset="FLUX")
  SWCData<-getSWC()
  fluxDat$week<-week(fluxDat$date)
  fluxDat$week[fluxDat$week==53]<-52
  SWCData$week[SWCData$week==53]<-52
  
  #mean LAI - could not find exact dates so using yearly average
  LAI<-aggregate(standDat$lai~standDat$year,FUN=mean)
  names(LAI)<-c("year","LAI")
  
  #tree density
  treeDen<-aggregate(standDat$density_treeha~standDat$year,FUN=mean)
  names(treeDen)<-c("year","treeDen")
  
  #DBH
  dbh<-aggregate(standDat$dbhBA_cm~standDat$year,FUN=mean)
  names(dbh)<-c("year","dbh")
  
  #carbon and nitrogen measured at differing soil depths - currently taking average
  massExt<-function(soildata,varMass){
    sProfM<- ((soildata$lowerDepth_cm - soildata$upperDepth_cm))*0.01
    t_ha<-10000*sProfM*soildata$density_gcm3*varMass
    return(t_ha)
  }
  totC<-sum(massExt(soildata,soildata$c_percent),na.rm=TRUE)
  totN<-sum(massExt(soildata,soildata$n_percent),na.rm=TRUE)
  
  
  
  
  convrtUnit=(12.011 * 24 * 60 * 60)/1000000 # convert to grams per day of C
  
  obsDat<-fluxDat%>%
    group_by(year,mo)%>%
    summarise(GPP=mean(gppDtCutRef_umolCO2m2s1*convrtUnit),NEE=mean(neeCutRef_umolCO2m2s1*convrtUnit),reco=mean(recoDtCutRef_umolCO2m2s1*convrtUnit),
              month=median(mo))%>%
    mutate(timestamp=as.POSIXct(paste0(year,"-",month,"-01")))
  
  obsDatSWC<-obsDat%>%filter(year>=1997)
  
  obsDatSWC$swc<-pull(SWCData%>%
                        filter(year>=1996)%>%
                        group_by(year,month)%>%
                        summarise(swc=mean(swc))%>%
                        select(swc))
  
  
  coefVar<-0.2
  
  predPos  <- intvsS[,1]$`89%`*modif + 2  * sapply( 1:length(obsDat$GPP), function(i) max( coefVar* abs(obsDat$GPP[i]),0.05))
  predNeg  <- intvsS[,1]$`11%`*modif - 2 * sapply( 1:length(obsDat$GPP), function(i) max( coefVar* abs(obsDat$GPP[i]),0.05))
  predm  <- intvsS[,1]$`50%`*modif# - 2 * 0.3
  

  
  gpp<-ggplot()+theme_bw()+
    geom_line(data=obsDat,aes(x=timestamp,y=predm),colour="black",size=1)+
    geom_point(data=obsDat,aes(x=timestamp,y=GPP),colour="red",size=2)+
      geom_ribbon(aes(x=obsDat$timestamp,ymin=predNeg,ymax=predPos),alpha=0.5,fill="orange")+
    scale_x_datetime(limits=c(as.POSIXct("1996-01-01",tz="GMT"),as.POSIXct("2015-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("GPP [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  
  predPosNEE  <- intvsS[,2]$`89%`*modif + 2  * sapply( 1:length(obsDat$NEE), function(i) max( coefVar* abs(obsDat$NEE[i]),0.05))
  predNegNEE  <- intvsS[,2]$`11%`*modif - 2 * sapply( 1:length(obsDat$NEE), function(i) max( coefVar* abs(obsDat$NEE[i]),0.05))
  predmNEE  <- intvsS[,2]$`50%`*modif# - 2 * 0.3

  nee<-ggplot()+theme_bw()+
    geom_line(data=obsDat,aes(x=timestamp,y=predmNEE),colour="black",size=1)+
    geom_point(data=obsDat,aes(x=timestamp,y=NEE),colour="red",size=2)+
      geom_ribbon(aes(x=obsDat$timestamp,ymin=predNegNEE,ymax=predPosNEE),alpha=0.5,fill="orange")+
    scale_x_datetime(limits=c(as.POSIXct("1996-01-01",tz="GMT"),as.POSIXct("2015-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("NEE [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  predPosSWC  <- intvsS[,3]$`89%`[13:228] + 2  * sapply( 1:length(obsDatSWC$swc), function(i) max( coefVar* abs(obsDatSWC$swc[i]),0.05))
  predNegSWC  <- intvsS[,3]$`11%`[13:228] - 2 * sapply( 1:length(obsDatSWC$swc), function(i) max( coefVar* abs(obsDatSWC$swc[i]),0.05))
  predmSWC  <- intvsS[,3]$`50%`[13:228]# - 2 * 0.3

  swc<-ggplot()+theme_bw()+
    geom_line(data=obsDatSWC,aes(x=timestamp,y=predmSWC),colour="black",size=1)+
    geom_point(data=obsDatSWC,aes(x=timestamp,y=swc),colour="red",size=2)+
    geom_ribbon(aes(x=obsDatSWC$timestamp,ymin=predNegSWC,ymax=predPosSWC),alpha=0.5,fill="orange")+
    scale_x_datetime(limits=c(as.POSIXct("1996-01-01",tz="GMT"),as.POSIXct("2015-01-01",tz="GMT")))+    
    labs(x="Year",y=("SWC"))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  coefVar2<-0.1
  LAIdat<-data.frame(LAI=intvsS[,8]$`50%`,lower=intvsS[,8]$`11%`,upper=intvsS[,8]$`89%`,year=obsDat$year)%>%group_by(year)%>%summarise(LAI=mean(LAI),upper=mean(upper),lower=mean(lower))%>%filter(year<2012)
 LAIdat$obs<-LAI$LAI[-1]
   LAIdat$upper  <- LAIdat$upper+ 2  * sapply( 1:length(LAIdat$obs), function(i) max( coefVar2* abs(LAIdat$obs[i]),0.05))
  LAIdat$lower  <- LAIdat$lower - 2 * sapply( 1:length(LAIdat$obs), function(i) max( coefVar2* abs(LAIdat$obs[i]),0.05))
  LAIdat$LAI  <- LAIdat$LAI# - 2 * 0.3
  
  
  
  LAI<-ggplot()+theme_bw()+
    geom_line(data=LAIdat,aes(x=year,y=LAI),colour="black",size=1)+
    geom_point(data=LAIdat,aes(x=year,y=obs),colour="red",size=2)+
    geom_ribbon(data=LAIdat,aes(x=year,ymin=lower,ymax=upper),alpha=0.5,fill="orange")+
    labs(x="Year",y=("LAI"))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  Ndat<-data.frame(N=intvsS[,7]$`50%`,lower=intvsS[,7]$`11%`,upper=intvsS[,7]$`89%`,year=obsDat$year)%>%group_by(year)%>%summarise(N=mean(N),upper=mean(upper),lower=mean(lower))%>%filter(year<2012)
  Ndat$obs<-treeDen$treeDen[-1]
  Ndat$upper  <- Ndat$upper+ 2  * sapply( 1:length(Ndat$obs), function(i) max( coefVar2* abs(Ndat$obs[i]),0.05))
  Ndat$lower  <- Ndat$lower - 2 * sapply( 1:length(Ndat$obs), function(i) max( coefVar2* abs(Ndat$obs[i]),0.05))
  Ndat$N  <- Ndat$N# - 2 * 0.3
  
  Nplot<-ggplot()+theme_bw()+
    geom_line(data=Ndat,aes(x=year,y=N),colour="black",size=1)+
    geom_point(data=Ndat,aes(x=year,y=obs),colour="red",size=2)+
    geom_ribbon(data=Ndat,aes(x=year,ymin=lower,ymax=upper),alpha=0.5,fill="orange")+
    labs(x="Year",y=("N"))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  dbhdat<-data.frame(N=intvsS[,9]$`50%`,lower=intvsS[,9]$`11%`,upper=intvsS[,9]$`89%`,year=obsDat$year)%>%group_by(year)%>%summarise(N=mean(N),upper=mean(upper),lower=mean(lower))%>%filter(year<2012)
  dbhdat$obs<-dbh$dbh[-1]
  dbhdat$upper  <- dbhdat$upper+ 2  * sapply( 1:length(dbhdat$obs), function(i) max( coefVar2* abs(dbhdat$obs[i]),0.05))
  dbhdat$lower  <- dbhdat$lower - 2 * sapply( 1:length(dbhdat$obs), function(i) max( coefVar2* abs(dbhdat$obs[i]),0.05))
  dbhdat$N  <- dbhdat$N# - 2 * 0.3
  
  dbhplot<-ggplot()+theme_bw()+
    geom_line(data=dbhdat,aes(x=year,y=N),colour="black",size=1)+
    geom_point(data=dbhdat,aes(x=year,y=obs),colour="red",size=2)+
    geom_ribbon(data=dbhdat,aes(x=year,ymin=lower,ymax=upper),alpha=0.5,fill="orange")+
    labs(x="Year",y=("DBH (cm)"))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  

  
  
  gpp1<-ggpubr::ggarrange(gpp,nee,swc,LAI,dbhplot,Nplot)
  
  return(list(gpp1))
}


############################################################################################



#'plotResultsNewMonthly Function to plot the model against the Harwood data.
#'@return pine plots
#'@export
plotResultsNewMonthly <- function(df,outSample,ShortTS=F,numSamps=500){
  library(matrixStats)
  library(future)
  outSample<-sample_n(outSample,numSamps)
  codM<-miscTools::colMedians(as.data.frame(outSample))
  sitka[.GlobalEnv$nm]<-codM[.GlobalEnv$nm]
  df<-do.call(SoNWaL,sitka)
  
  df <- df[c(2:nrow(df)),]
  if(nrow(df)>600) df$week<-c(1:53) else df$week<-1
  df <- df%>%mutate(timestamp = as.POSIXct(paste(sprintf("%02d",Year),sprintf("%02d",Month),sprintf("%02d",1),sep="-"),tz="GMT")) 
  
  flxdata<-flxdata_daily%>%
    mutate(week=week(as.Date(yday, origin = paste0(year,"-01-01"))))%>%
    mutate(month=month(as.Date(yday, origin = paste0(year,"-01-01"))))
  
  dataX<- flxdata%>% 
    group_by(year,month) %>%
    dplyr::summarise(gppOb=mean(gpp),gppsum=sum(gpp),nppOb=mean(npp),nppsum=sum(npp),etOb=mean(et),recoOb=mean(reco),rsOb=mean(rs),
                     swcOb=mean(swc),neeOb=mean(nee))%>%mutate(cumGppObs=cumsum(gppsum),cumNppObs=cumsum(nppsum))
  df2<- (df%>% 
           group_by(Year,Month) %>%
           dplyr::summarise(GPPsum=sum(GPP*7),NPPsum=sum(NPP*7),GPP=mean(GPP),EvapTransp=mean(EvapTransp),Reco=mean(Reco),Rs=mean(Rs),NPP=mean(NPP),
                            volSWC_rz=mean(volSWC_rz),NEE=mean(NEE),timestamp=median(timestamp),LAI=mean(LAI),t.proj=median(t.proj)))
  df3<-df2%>%filter(Year>=1971)
  df2<-df2%>%filter(Year>=2015)
  
  modif<-ifelse(nrow(df)<600,1.66,7.142857)
  nDays<-ifelse(nrow(df)<600,30,7)
  
  dataX<-dataX %>% right_join(df2, by=c("year"="Year","month"="Month"))
  
  dataX$simGpp<-df2$GPP*modif
  dataX$simCGpp<-dplyr::pull(df2%>%dplyr::mutate(gppC=cumsum(GPPsum*modif))%>%dplyr::select(gppC))
  dataX$simCNpp<-dplyr::pull(df2%>%dplyr::mutate(nppC=cumsum(NPPsum*modif))%>%dplyr::select(nppC))
  
  dataX$simReco<-df2$Reco*modif
  dataX$simNEE<-df2$NEE*modif
  dataX$simEtransp<-df2$EvapTransp
  dataX$simswc<-df2$volSWC_rz
  dataX$simRs<-df2$Rs*modif
  dataX$timestamp<-df2$timestamp
  dataX$month<-month(dataX$timestamp)  
  

  
  
  runModel<- function(p){
    require("dplyr")
    sitka[nm]<-p
      res<- do.call(SoNWaL,sitka)%>%
                   group_by(Year,Month)%>%
        dplyr::summarise(GPPsum=sum(GPP*nDays),NPPsum=sum(NPP*nDays),GPP=mean(GPP),NEE=mean(NEE),volSWC_rz=mean(volSWC_rz),EvapTransp=mean(EvapTransp)/nDays,Reco=mean(Reco),Rs=mean(Rs),N=mean(N),LAI=mean(LAI),dg=mean(dg),totC=mean(totC),totN=mean(totN),NPP=mean(NPP),alphaAn=mean(Reco/Rs))
   
    return(res)
  }
  
  getIntv<-function(paramName,modLst){
    simDat =   do.call(cbind,lapply(modLst, "[",  paramName))
   res<-as.data.frame(rowQuantiles(as.matrix(simDat), probs = c(0.11, 0.5, 0.89)))
  return(res)
   }
  

  outSample <- split(outSample, seq(nrow(outSample)))
  outResX<- lapply(outSample, runModel) 
  outResX<-lapply(outResX, function(x) filter(x, Year <= 2018))
 # outResX<-lapply(outResX, function(x) filter(x, Year >= 2015))
  outRes<-lapply(outResX, function(x) filter(x, Year >= 2015))

  
paramName<-list("GPP","NEE","volSWC_rz","EvapTransp","Reco","Rs","N","LAI","dg","totC","totN","NPP","alphaAn","GPPsum","NPPsum")
intvsS<-mapply(getIntv,paramName,MoreArgs = list(modLst=outRes))

  data<-flxdata_daily%>%
    mutate(grp=month(as.Date(flxdata_daily$yday, origin = paste0(flxdata_daily$year,"-01-01"))))
  sdMin<-data%>% group_by(year,grp) %>%
    dplyr::summarise(sdgpp=mean(gpp),sdnpp=mean(npp),sdnee=mean(nee),sdreco=mean(reco),
                     sdrs=mean(rs),sdet=mean(et),sdswc=mean(swc,na.rm=T))
  
  
  coefVar<-0.2
  
  predPos  <- intvsS[,1]$`89%`*modif + 2  * sapply( 1:length(dataX$gppOb), function(i) max( coefVar* abs(dataX$gppOb[i]),0.05))
  predNeg  <- intvsS[,1]$`11%`*modif - 2 * sapply( 1:length(dataX$gppOb), function(i) max( coefVar* abs(dataX$gppOb[i]),0.05))
  predm  <- intvsS[,1]$`50%`*modif# - 2 * 0.3
  
  predNPPPos  <- intvsS[,12]$`89%`*modif + 2  * sapply( 1:length(sdMin$sdnpp), function(i) max( coefVar* abs(sdMin$sdnpp[i]),0.05))
  predNPPNeg  <- intvsS[,12]$`11%`*modif - 2 * sapply( 1:length(sdMin$sdnpp), function(i) max( coefVar* abs(sdMin$sdnpp[i]),0.05))
  predmNPP  <- intvsS[,12]$`50%`*modif# - 2 * 0.3
  
  predPosSum  <- intvsS[,14]$`89%`*modif + 2  * sapply( 1:length(dataX$gppsum), function(i) max( coefVar* abs(dataX$gppsum[i]),0.05))
  predNegSum  <- intvsS[,14]$`11%`*modif - 2 * sapply( 1:length(dataX$gppsum), function(i) max( coefVar* abs(dataX$gppsum[i]),0.05))
  predmSum  <- intvsS[,14]$`50%`*modif# - 2 * 0.3
  
  predNPPPosSum  <- intvsS[,15]$`89%`*modif + 2  * sapply( 1:length(dataX$nppsum), function(i) max( coefVar* abs(dataX$nppsum[i]),0.05))
  predNPPNegSum  <- intvsS[,15]$`11%`*modif - 2 * sapply( 1:length(dataX$nppsum), function(i) max( coefVar* abs(dataX$nppsum[i]),0.05))
  predmNPPSum  <- intvsS[,15]$`50%`*modif# - 2 * 0.3
  
  
  
  
  dataX$simCGppPos<-predPosSum
  dataX$simCGppNeg<-predNegSum
  dataX$simCGpp<-predmSum
  
  dataX$simCGppNeg<-dplyr::pull(dataX%>%dplyr::mutate(GPPCneg=cumsum(simCGppNeg))%>%dplyr::select(GPPCneg))
  dataX$simCGppPos<-dplyr::pull(dataX%>%dplyr::mutate(GPPCpos=cumsum(simCGppPos))%>%dplyr::select(GPPCpos))
  dataX$simCGpp<-dplyr::pull(dataX%>%dplyr::mutate(GPPCpos=cumsum(simCGpp))%>%dplyr::select(GPPCpos))
  
  dataX$simCNppPos<-predNPPPosSum
  dataX$simCNppNeg<-predNPPNegSum
  dataX$predmNPPSum<-predmNPPSum
  
  dataX$simCNppNeg<-dplyr::pull(dataX%>%dplyr::mutate(NPPCneg=cumsum(simCNppNeg))%>%dplyr::select(NPPCneg))
  dataX$simCNppPos<-dplyr::pull(dataX%>%dplyr::mutate(NPPCpos=cumsum(simCNppPos))%>%dplyr::select(NPPCpos))
  dataX$simCNpp<-dplyr::pull(dataX%>%dplyr::mutate(NPPC=cumsum(predmNPPSum))%>%dplyr::select(NPPC))
  
  
  
  
  # dataX2<-dataX%>% 
  #   group_by(year,month) %>%
  #   dplyr::summarise(simGpp=mean(simGpp),timestamp=median(timestamp),gpp=mean(gpp))
  #
  gpp1<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predm),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=timestamp, y=gppOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNeg,ymax=predPos,x=dataX$timestamp),fill="orange",alpha=0.3)+
    ## geom_ribbon(data=data,aes(x=timestamp,ymin=gpp-gpp.sd,ymax=gp+gpp.sd),alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("GPP [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  gppC<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=simCGpp,group=year),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=timestamp, y=cumGppObs),colour="black",size=2)+
    geom_ribbon(data=dataX,aes(ymin=simCGppNeg,ymax=simCGppPos,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Cumulative GPP [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  nppC<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=simCNpp,group=year),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=timestamp, y=cumNppObs),colour="black",size=2)+
    geom_ribbon(data=dataX,aes(ymin=simCNppNeg,ymax=simCNppPos,x=dataX$timestamp),fill="orange",alpha=0.3)+
        scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Cumlative NPP [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  predPosEtransp  <- intvsS[,4]$`89%`  + 2  * sapply( 1:length(sdMin$sdet), function(i) max( coefVar* abs(sdMin$sdet[i]),0.01))
  predNegEtransp  <- intvsS[,4]$`11%` - 2  * sapply( 1:length(sdMin$sdet), function(i) max( coefVar* abs(sdMin$sdet[i]),0.01))
  predmEtransp  <-   intvsS[,4]$`50%` # - 2 * sd(dataX$gpp)
  
  etrans<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmEtransp),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=etOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegEtransp,ymax=predPosEtransp,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste(E[t]," [mm]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  predPosSWC  <- intvsS[,3]$`89%`  + 2  * sapply( 1:length(sdMin$sdswc), function(i) max( 0.05* abs(sdMin$sdswc[i]),0.01))
  predNegSWC  <- intvsS[,3]$`11%` - 2   * sapply( 1:length(sdMin$sdswc), function(i) max( 0.05* abs(sdMin$sdswc[i]),0.01))
  predmSWC  <-   intvsS[,3]$`50%` # - 2 2 * sd(dataX$gpp)
#  newSWC<-clm_df_full%>% suzie- debug attempt swapped in to read in flux data
#    filter(Year>=2015)%>%
#    group_by(Year,Month)%>%
#    summarise(swc=mean(SWC))
  newSWC <- flxdata%>%
    filter(year>=2015)%>%
    group_by(year,month)%>%
    summarise(swc=mean(swc))
  
  swcPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmSWC),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=newSWC$swc),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegSWC,ymax=predPosSWC,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y="SWC")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  predPosNEE  <-intvsS[,2]$`89%`*modif  + 2 * sapply( 1:length(sdMin$sdnee), function(i) max( coefVar* abs(sdMin$sdnee[i]),0.05))
  predNegNEE  <-intvsS[,2]$`11%`*modif - 2  * sapply( 1:length(sdMin$sdnee), function(i) max( coefVar* abs(sdMin$sdnee[i]),0.05))
  predmNEE  <- intvsS[,2]$`50%`*modif # - 2  sd(dataX$gpp)
  
  
  
  
  NEEPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmNEE),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=neeOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegNEE,ymax=predPosNEE,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("NEE [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  predPosreco  <- intvsS[,5]$`89%`*modif  + 2 * sapply( 1:length(sdMin$sdreco), function(i) max( coefVar* abs(sdMin$sdreco[i]),0.1))
  predNegreco  <- intvsS[,5]$`11%`*modif - 2  * sapply( 1:length(sdMin$sdreco), function(i) max( coefVar* abs(sdMin$sdreco[i]),0.1))
  predmreco  <-   intvsS[,5]$`50%`*modif # - 2  * sd(dataX$gpp)
  
  recoPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmreco),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=recoOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegreco,ymax=predPosreco,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Reco [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  
  predPosRs  <- intvsS[,6]$`89%`*modif  + 2 * sapply( 1:length(sdMin$sdrs), function(i) max( coefVar* abs(sdMin$sdrs[i]),0.3))
  predNegRs  <- intvsS[,6]$`11%`*modif - 2  * sapply( 1:length(sdMin$sdrs), function(i) max( coefVar* abs(sdMin$sdrs[i]),0.3))
  predmRs  <-   intvsS[,6]$`50%`*modif # - 2  * sd(dataX$gpp)
  
  rsPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmRs),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=rsOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegRs,ymax=predPosRs,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Rs [gC"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  paramName<-list("GPP","NEE","volSWC_rz","EvapTransp","Reco","Rs","N","LAI","dg","totC","totN","NPP","alphaAn","GPPsum","NPPsum")
  intvsS2<-mapply(getIntv,paramName,MoreArgs = list(modLst=outResX))
  
  dfLeaf<-df2%>%group_by(Year,Month) %>%
    dplyr::summarise(LAI=mean(LAI),t.proj=median(t.proj))
# suzie- getting error: non-numeric argument to binary operator- issue with the dt
  dt=12
# suzie re-defined dt as 12 this allows the plotting function to run ** need to check why?
  leaf<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2015&Month==8)$Year+dplyr::filter(dfLeaf,Year==2015&Month==8)$Month/dt,dplyr::filter(dfLeaf,Year==2015&Month==8)$t.proj,5.7),
    c(dplyr::filter(dfLeaf,Year==2018&Month==8)$Year+dplyr::filter(dfLeaf,Year==2018&Month==8)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==8)$t.proj,5.56)
  )
  )
  names(leaf)<-c("year","age","lai")
  
  stemNo<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2018&Month==8)$Year+dplyr::filter(dfLeaf,Year==2018&Month==8)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==8)$t.proj,1348)
  )
  )
  names(stemNo)<-c("year","age","stem")
  
  DBH<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2018&Month==8)$Year+dplyr::filter(dfLeaf,Year==2018&Month==8)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==8)$t.proj,24.1)
  )
  )
  names(DBH)<-c("year","age","dbh")
  
  Wr<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2015&Month==8)$Year+dplyr::filter(dfLeaf,Year==2015&Month==8)$Month/dt,dplyr::filter(dfLeaf,Year==2015&Month==8)$t.proj,4.88)
  )
  )
  names(Wr)<-c("year","age","wr")
  
  difRoots<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2015&Month==8)$Year+dplyr::filter(dfLeaf,Year==2015&Month==8)$Month/dt,dplyr::filter(dfLeaf,Year==2015&Month==8)$t.proj,0.53)
  )
  )
  names(difRoots)<-c("year","age","difroots")
  
  totC<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2018&Month==7)$Year+dplyr::filter(dfLeaf,Year==2018&Month==7)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==7)$t.proj,214.76)
  )
  )
  names(totC)<-c("year","age","totC")
  
  totN<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2018&Month==7)$Year+dplyr::filter(dfLeaf,Year==2018&Month==7)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==7)$t.proj,7.15)
  )
  )
  names(totN)<-c("year","age","totN")
  
  alphaN<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2015&Month==1)$Year+dplyr::filter(dfLeaf,Year==2015&Month==1)$Month/dt,dplyr::filter(dfLeaf,Year==2015&Month==1)$t.proj,2.06),
    c(dplyr::filter(dfLeaf,Year==2016&Month==1)$Year+dplyr::filter(dfLeaf,Year==2016&Month==1)$Month/dt,dplyr::filter(dfLeaf,Year==2016&Month==1)$t.proj,2.62),
    c(dplyr::filter(dfLeaf,Year==2017&Month==1)$Year+dplyr::filter(dfLeaf,Year==2017&Month==1)$Month/dt,dplyr::filter(dfLeaf,Year==2017&Month==1)$t.proj,2.69),
    c(dplyr::filter(dfLeaf,Year==2018&Month==1)$Year+dplyr::filter(dfLeaf,Year==2018&Month==1)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==1)$t.proj,3.11)
    )
  )
  names(alphaN)<-c("year","age","alphaN")
  
  predPosN  <- intvsS2[,7]$`89%` # + 2 * 0.01
  predNegN   <- intvsS2[,7]$`11%` #- 2  * 0.01
  predmN   <-   intvsS2[,7]$`50%`
  
  predPosLAI  <- intvsS2[,8]$`89%` # + 2 * 0.01
  predNegLAI   <- intvsS2[,8]$`11%` #- 2  * 0.01
  predmLAI   <-   intvsS2[,8]$`50%`
  
  predPosDG  <- intvsS2[,9]$`89%` # + 2 * 0.01
  predNegDG    <- intvsS2[,9]$`11%` #- 2  * 0.01
  predmDG    <-   intvsS2[,9]$`50%`
  

  predPosTotC  <- intvsS2[,10]$`89%` # + 2 * 0.125
  predNegTotC    <- intvsS2[,10]$`11%` #- 2  * 0.125
  predmTotC    <-   intvsS2[,10]$`50%`
  
  predPosTotN  <- intvsS2[,11]$`89%`  #+ 2 * 0.125
  predNegTotN    <- intvsS2[,11]$`11%` #- 2  * 0.125
  predmTotN    <-   intvsS2[,11]$`50%`

  dataX$predmAlphaAn<-intvsS[,13]$`50%`
  dataX$predPosAlphaAn<-intvsS[,13]$`89%`
  dataX$predNegAlphaAn<-intvsS[,13]$`11%`
  alphaAnMean<-dataX%>%group_by(year)%>%summarise(predmAlphaAn=mean(predmAlphaAn),predPosAlphaAn=mean(predPosAlphaAn),predNegAlphaAn=mean(predNegAlphaAn))
  
  pLAI<-ggplot()+theme_bw()+
    geom_line(data=df3,aes(x=Year+Month/dt,y=predmLAI))+
    geom_point(data=leaf,aes(x=year,y=lai),shape=16,size=3,colour="red")+
    geom_pointrange(data=leaf,aes(x=year,y=lai,ymin=lai-lai*0.1, ymax=lai+lai*0.1),colour="red")+
    geom_ribbon(aes(ymin=predNegLAI,ymax=predPosLAI,x=df3$Year+df3$Month/dt),fill="orange",alpha=0.3)+
    labs(x="Year",y=expression(paste("LAI"," ","[",cm^-2," ",cm^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pStemNo<-ggplot()+theme_bw()+
    geom_line(data=df3,aes(x=Year+Month/dt,y=predmN))+
    geom_point(data=stemNo,aes(x=year,y=stem),shape=16,size=3,colour="red")+
    geom_pointrange(data=stemNo,aes(x=year,y=stem,ymin=stem-stem*0.1, ymax=stem+stem*0.1),colour="red")+
    geom_ribbon(aes(ymin=predNegN,ymax=predPosN,x=df3$Year+df3$Month/dt),fill="orange",alpha=0.3)+
    labs(x="Year",y=expression(paste("N"," ","[",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pDBH<-ggplot()+theme_bw()+
    geom_line(data=df3,aes(x=Year+Month/dt,y=predmDG))+
    geom_point(data=DBH,aes(x=year,y=dbh),shape=16,size=3,colour="red")+
    geom_pointrange(data=DBH,aes(x=year,y=dbh,ymin=dbh-dbh*0.1, ymax=dbh+dbh*0.1),colour="red")+
    geom_ribbon(aes(ymin=predNegDG,ymax=predPosDG,x=df3$Year+df3$Month/dt),fill="orange",alpha=0.3)+
    labs(x="Year",y="DBH [cm]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pWr<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=Wr))+
    geom_point(data=Wr,aes(x=year,y=wr),shape=16,size=3,colour="red")+
    labs(x="Year",y="Wr [tDM/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pdifRoots<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=difRoots))+
    geom_point(data=difRoots,aes(x=year,y=difroots),shape=16,size=3,colour="red")+
    labs(x="Year",y="difRoots [tDM/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))    
  
  ptotC<-ggplot()+theme_bw()+
    geom_line(data=df3,aes(x=Year+Month/dt,y=predmTotC))+
    geom_point(data=totC,aes(x=year,y=totC),shape=16,size=3,colour="red")+
    geom_pointrange(data=totC,aes(x=year,y=totC,ymin=totC-totC*0.5, ymax=totC+totC*0.5),colour="red")+
    geom_ribbon(aes(ymin=predNegTotC,ymax=predPosTotC,x=df3$Year+df3$Month/dt),fill="orange",alpha=0.3)+
    labs(x="Year",y="totC [t/ha]")+ 
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  ptotN<-ggplot()+theme_bw()+
    geom_line(data=df3,aes(x=Year+Month/dt,y=predmTotN))+
    geom_point(data=totN,aes(x=year,y=totN),shape=16,size=3,colour="red")+
    geom_pointrange(data=totN,aes(x=year,y=totN,ymin=totN-totN*0.5, ymax=totN+totN*0.5),colour="red")+
    geom_ribbon(aes(ymin=predNegTotN,ymax=predPosTotN,x=df3$Year+df3$Month/dt),fill="orange",alpha=0.3)+
    labs(x="Year",y="totN [t/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  ggAlphaAn<-ggplot()+theme_bw()+
    geom_line(data=alphaAnMean,aes(x=year,y=predmAlphaAn))+
    geom_point(data=alphaN,aes(x=year-0.083,y=alphaN),shape=16,size=3,colour="red")+
    geom_pointrange(data=alphaN,aes(x=year-0.083,y=alphaN,ymin=alphaN-alphaN*0.3, ymax=alphaN+alphaN*0.3),colour="red")+
    geom_ribbon(data=alphaAnMean,aes(ymin=predNegAlphaAn,ymax=predPosAlphaAn,x=year),fill="orange",alpha=0.3)+
    labs(x="Year",y="AlphaAnnual")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  return(list(gpp1,swcPlot,NEEPlot,etrans,recoPlot,rsPlot,gppC,nppC, pLAI,pStemNo,pDBH,pWr,ptotC,ptotN,ggAlphaAn))
}




