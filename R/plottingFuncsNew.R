###plotting functions for FR3pgn
library(rlang)
runModel<- function(p){
  pine[.GlobalEnv$nm]<-p
  pull(do.call(fr3PGDN,pine)%>%
    filter(Year>=1996)%>%
     group_by(Year,Month) %>%
   dplyr::summarise(GPP=mean(!!sym(prm))))
}

## Function to plot the model against the Harwood data.
plotResultsPine <- function(df,ShortTS=F){
  
  ##if plyr is loaded before dplyr can cause problems with groupby
  dt=12
  

  df$week<-c(1:53)
  
  df<-df[-1,]
  convrtUnit=(12.011 * 24 * 60 * 60)/1000000 # convert to grams per day of C
  fluxDat <- getData(site=sites[site],dataset="FLUX")
  SWCData<-getSWC()
  fluxDat$week<-week(fluxDat$date)
  fluxDat$week[fluxDat$week==53]<-52
  SWCData$week[SWCData$week==53]<-52
  
  obsDat<-fluxDat%>%
    group_by(year,mo)%>%
    summarise(GPP=mean(gppDtCutRef_umolCO2m2s1*convrtUnit),NEE=mean(neeCutRef_umolCO2m2s1*convrtUnit),reco=mean(recoDtCutRef_umolCO2m2s1*convrtUnit),
             month=median(mo))%>%
    mutate(timestamp=as.POSIXct(paste0(year,"-",month,"-01")))
  
    obsDat$simGPP<-pull(df%>%
    filter(Year>=1996)%>%
    group_by(Year,Month)%>%
    summarise(gpp=mean(GPP*100/7))%>%
    select(gpp))

    obsDat$simGPP<-pull(df%>%
                          filter(Year>=1996)%>%
                          group_by(Year,Month)%>%
                          summarise(gpp=mean(GPP*100/7))%>%
                          select(gpp))
    
    obsDat$simNEE<-pull(df%>%
                          filter(Year>=1996)%>%
                          group_by(Year,Month)%>%
                          summarise(nee=mean(NEE*100/7))%>%
                          select(nee))
    
    obsDat$simSWC<-pull(df%>%
                          filter(Year>=1996)%>%
                          group_by(Year,Month)%>%
                          summarise(swc=mean(volSWC_rz))%>%
                          select(swc))
    
    
    nmc = nrow(out$chain[[1]])
    outSample   <- getSample(out,start=nmc/2)
    numSamples = min(1000, nrow(outSample))
    
    runMltMod<-function(prm){
      print(prm)
      prm<<-prm
    pred     <- getPredictiveIntervals(parMatrix=outSample, model=runModel, numSamples=numSamples, quantiles=c(0.025,0.5,0.975))
    return(pred)
    }
    

    plan(multisession, workers = 4,.cleanup = TRUE)
    intvs<-map(c("GPP","NEE","volSWC_rz"),runMltMod)
    
    
    
    
    
    predPosGPP  <- intvs[[1]]$posteriorPredictiveCredibleInterval[3,]*100/7 #+ .5  * sd(obsDat$GPP)
    predNegGPP  <- intvs[[1]]$posteriorPredictiveCredibleInterval[1,]*100/7 #- .5 * sd(obsDat$GPP)

  gpp<-ggplot()+theme_bw()+
    geom_line(data=obsDat,aes(x=timestamp,y=simGPP),colour="black",size=1)+
    geom_point(data=obsDat,aes(x=timestamp,y=GPP),colour="red",size=2)+
   #  geom_ribbon(aes(x=obsDat$timestamp,ymin=predNegGPP,ymax=predPosGPP),alpha=0.5,fill="blue")+
    scale_x_datetime(limits=c(as.POSIXct("1996-01-01",tz="GMT"),as.POSIXct("2015-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("GPP [gDM"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  predPosNEE  <- intvs[[2]]$posteriorPredictiveCredibleInterval[3,]*100/7# + 2  * sd(obsDat$NEE)
  predNegNEE  <- intvs[[2]]$posteriorPredictiveCredibleInterval[1,]*100/7# - 2 * sd(obsDat$NEE)
  
  nee<-ggplot()+theme_bw()+
    geom_line(data=obsDat,aes(x=timestamp,y=simNEE),colour="black",size=1)+
    geom_point(data=obsDat,aes(x=timestamp,y=NEE),colour="red",size=2)+
   # geom_ribbon(aes(x=obsDat$timestamp,ymin=predNegNEE,ymax=predPosNEE),alpha=0.5,fill="blue")+
    scale_x_datetime(limits=c(as.POSIXct("1996-01-01",tz="GMT"),as.POSIXct("2015-01-01",tz="GMT")))+       
    labs(x="Year",y=expression(paste("NEE [gDM"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  obsDatSWC<-obsDat%>%filter(year>=1997)
  
  obsDatSWC$swc<-pull(SWCData%>%
                     filter(year>=1996)%>%
                     group_by(year,month)%>%
                     summarise(swc=mean(swc))%>%
                     select(swc))
  
  predPosSWC  <- intvs[[3]]$posteriorPredictiveCredibleInterval[3,c(13:228)] #+ 2  * sd(obsDat$NEE)
  predNegSWC  <- intvs[[3]]$posteriorPredictiveCredibleInterval[1,c(13:228)] #- 2 * sd(obsDat$NEE)
  predSWC  <- intvs[[3]]$posteriorPredictiveCredibleInterval[2,c(13:228)] #- 2 * sd(obsDat$NEE)
  
  swc<-ggplot()+theme_bw()+
    geom_line(data=obsDatSWC,aes(x=timestamp,y=predSWC),colour="black",size=1)+
    geom_point(data=obsDatSWC,aes(x=timestamp,y=swc),colour="red",size=2)+
   # geom_ribbon(aes(x=obsDatSWC$timestamp,ymin=predNegSWC,ymax=predPosSWC),alpha=0.5,fill="blue")+
    scale_x_datetime(limits=c(as.POSIXct("1996-01-01",tz="GMT"),as.POSIXct("2015-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("swc",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  gpp1<-ggarrange(gpp,nee,swc)
  
  return(list(gpp1))
}


## Function to find the yield class of the stand based on stand dominant height (assuming it is the same as top height)
YC.finder <- function(HT,AGE) 
{
  YC.site = ((HT/(1-exp(-0.033329*AGE))^1.821054)-14.856317)/1.425397
  if(YC.site>24) 24
  else if (YC.site<4) 4
  else unlist(sub("\\]","",unlist(strsplit(as.character(cut(YC.site,breaks=c(6,8,10,12,14,16,18,20,22,24),right=T)),split="\\,"))[2]))
}

############################################################################################





## Function to plot the model against the Harwood data.
plotResults <- function(df,out,ShortTS=F){
  
  ##if plyr is loaded before dplyr can cause problems with groupby
  dt=12
df <- df[c(2:nrow(df)),]
 if(nrow(df)>600) df$week<-c(1:53) else df$week<-1
  df <- df %>% dplyr::group_by(Year)%>%mutate(cumGPP = cumsum(GPP),
                                              cumNPP = cumsum(NPP),
                                              timestamp = as.POSIXct(paste(sprintf("%02d",Year),sprintf("%02d",Month),sprintf("%02d",1),sep="-"),tz="GMT")) 
  df2<-df%>%filter(Year>=2015)
  
  flxdata<-flxdata_daily%>%
    mutate(week=week(as.Date(yday, origin = paste0(year,"-01-01"))))%>%
    mutate(month=month(as.Date(yday, origin = paste0(year,"-01-01"))))
  
  dataX<- flxdata%>% 
             group_by(year,month) %>%
             dplyr::summarise(gppOb=mean(gpp),nppOb=mean(npp),etOb=sum(et),recoOb=mean(reco),rsOb=mean(rs),
                              swcOb=mean(swc),neeOb=mean(nee))%>%mutate(cumGppObs=cumsum(gppOb),cumNppObs=cumsum(nppOb))
  df2<- (df2%>% 
             group_by(Year,Month) %>%
             dplyr::summarise(GPP=mean(GPP),Etransp=sum(Etransp),Reco=mean(Reco),Rs=mean(Rs),NPP=mean(NPP),
                              volSWC_rz=mean(volSWC_rz),NEE=mean(NEE),timestamp=median(timestamp),LAI=mean(LAI),t.proj=median(t.proj)))
  
  modif<-ifelse(nrow(df)<600,100/30,100/7)
  
  dataX<-dataX %>% right_join(df2, by=c("year"="Year","month"="Month"))
  
  dataX$simGpp<-df2$GPP*modif
  dataX$simCGpp<-pull(df2%>%mutate(gppC=cumsum(GPP*modif))%>%select(gppC))
  dataX$simCNpp<-pull(df2%>%mutate(nppC=cumsum(NPP*modif))%>%select(nppC))
  
  dataX$simReco<-df2$Reco*modif
  dataX$simNEE<-df2$NEE*modif
  dataX$simEtransp<-df2$Etransp
  dataX$simswc<-df2$volSWC_rz
  dataX$simRs<-df2$Rs*modif
  dataX$timestamp<-df2$timestamp
  dataX$month<-month(dataX$timestamp)  
  
  
  
 
 nmc = nrow(out$chain[[1]])
 outSample   <- getSample(out,start=nmc/2)
 numSamples = 50# min(1000, nrow(outSample))
 
 runModel<- function(p){
   sitka[.GlobalEnv$nm]<-p
   if(prm!="Etransp"){
  res<- pull(do.call(fr3PGDN,sitka)%>%
          filter(Year>=2015)%>%
          group_by(Year,Month) %>%
          dplyr::summarise(mean=mean(!!sym(prm)))%>%
          select(mean)
        )
   } else {
     res<- pull(do.call(fr3PGDN,sitka)%>%
                  filter(Year>=2015)%>%
                  group_by(Year,Month) %>%
                  dplyr::summarise(sum=sum(!!sym(prm)))%>%
                  select(sum))
   }
   return(res)
 }
 
 runMltMod <- function(prm=prm){
   print(prm)
   prm<<-prm
   pred     <- getPredictiveIntervals(parMatrix=outSample, model=runModel, numSamples=numSamples, quantiles=c(0.025,0.5,0.975))
   return(pred)
 }
 
 intvsS<-map(c("GPP","NEE","volSWC_rz","Etransp","Reco","Rs"),runMltMod)
 
 data<-flxdata_daily%>%
   mutate(grp=month(as.Date(flxdata_daily$yday, origin = paste0(flxdata_daily$year,"-01-01"))))
 sdMin<-data%>% group_by(year,grp) %>%
   dplyr::summarise(sdgpp=mean(gpp),sdnpp=mean(npp),sdnee=mean(nee),sdreco=mean(reco),
                    sdrs=mean(rs),sdet=sum(et),sdswc=mean(swc))
 
 
 coefVar<-0.1

predPos  <- intvsS[[1]]$posteriorPredictiveCredibleInterval[3,]*modif + 2  * sapply( 1:length(sdMin$sdgpp), function(i) max( 0.03* abs(sdMin$sdgpp[i]),0.05))
predNeg  <- intvsS[[1]]$posteriorPredictiveCredibleInterval[1,]*modif - 2 * sapply( 1:length(sdMin$sdgpp), function(i) max( 0.03* abs(sdMin$sdgpp[i]),0.05))
predm  <- intvsS[[1]]$posteriorPredictiveCredibleInterval[2,]*modif# - 2 * 0.3


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
    labs(x="Year",y=expression(paste("GPP [gDM"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  

gppC<-ggplot()+theme_bw()+
  geom_line(data=dataX,aes(x=timestamp,y=simCGpp,group=year),colour="purple",size=1)+
  geom_point(data=dataX,aes(x=timestamp, y=cumGppObs),colour="black",size=2)+
  ## geom_ribbon(data=data,aes(x=timestamp,ymin=gpp-gpp.sd,ymax=gp+gpp.sd),alpha=0.3)+
  scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
  labs(x="Year",y=expression(paste("Cumulative GPP [gDM"," ",m^-2,"]",sep="")))+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=14))


nppC<-ggplot()+theme_bw()+
  geom_line(data=dataX,aes(x=timestamp,y=simCNpp,group=year),colour="purple",size=1)+
  geom_point(data=dataX,aes(x=timestamp, y=cumNppObs),colour="black",size=2)+
  ## geom_ribbon(data=data,aes(x=timestamp,ymin=gpp-gpp.sd,ymax=gp+gpp.sd),alpha=0.3)+
  scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
  labs(x="Year",y=expression(paste("Cumlative NPP [gDM"," ",m^-2,"]",sep="")))+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=14))

predPosEtransp  <- intvsS[[4]]$posteriorPredictiveCredibleInterval[3,] + 2  * sapply( 1:length(sdMin$sdet), function(i) max( 0.2* abs(sdMin$sdet[i]),0.1))
predNegEtransp  <- intvsS[[4]]$posteriorPredictiveCredibleInterval[1,] - 2 * sapply( 1:length(sdMin$sdet), function(i) max( 0.2* abs(sdMin$sdet[i]),0.1))
predmEtransp  <- intvsS[[4]]$posteriorPredictiveCredibleInterval[2,]# - 2 * sd(dataX$gpp)
  
  etrans<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmEtransp),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=etOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegEtransp,ymax=predPosEtransp,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste(E[t]," [mm]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  predPosSWC  <- intvsS[[3]]$posteriorPredictiveCredibleInterval[3,] + 2  * sapply( 1:length(sdMin$sdswc), function(i) max( coefVar* abs(sdMin$sdswc[i]),0.01))
  predNegSWC  <- intvsS[[3]]$posteriorPredictiveCredibleInterval[1,] - 2 * sapply( 1:length(sdMin$sdswc), function(i) max( coefVar* abs(sdMin$sdswc[i]),0.01))
  predmSWC  <- intvsS[[3]]$posteriorPredictiveCredibleInterval[2,]# - 2 * sd(dataX$gpp)
  newSWC<-clm_df_full%>%
    filter(Year>=2015)%>%
    group_by(Year,Month)%>%
    summarise(swc=mean(SWC))
  
  
  swcPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmSWC),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=newSWC$swc),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegSWC,ymax=predPosSWC,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste(E[t]," [mm]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  predPosNEE  <- intvsS[[2]]$posteriorPredictiveCredibleInterval[3,]*modif + 2  * sapply( 1:length(sdMin$sdnee), function(i) max( 0.03* abs(sdMin$sdnee[i]),0.05))
  predNegNEE  <- intvsS[[2]]$posteriorPredictiveCredibleInterval[1,]*modif - 2 * sapply( 1:length(sdMin$sdnee), function(i) max( 0.03* abs(sdMin$sdnee[i]),0.05))
  predmNEE  <- intvsS[[2]]$posteriorPredictiveCredibleInterval[2,]*modif# - 2 * sd(dataX$gpp)
  
  NEEPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmNEE),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=neeOb),colour="black",size=2)+
   geom_ribbon(aes(ymin=predNegNEE,ymax=predPosNEE,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("NEE [gDM"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  predPosreco  <- intvsS[[5]]$posteriorPredictiveCredibleInterval[3,]*modif + 2  * sapply( 1:length(sdMin$sdreco), function(i) max( coefVar* abs(sdMin$sdreco[i]),0.1))
  predNegreco  <- intvsS[[5]]$posteriorPredictiveCredibleInterval[1,]*modif - 2 * sapply( 1:length(sdMin$sdreco), function(i) max( coefVar* abs(sdMin$sdreco[i]),0.1))
  predmreco  <- intvsS[[5]]$posteriorPredictiveCredibleInterval[2,]*modif# - 2 * sd(dataX$gpp)
  
  recoPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmreco),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=recoOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegreco,ymax=predPosreco,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Reco [gDM"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  
  predPosRs  <- intvsS[[6]]$posteriorPredictiveCredibleInterval[3,]*modif + 2  * sapply( 1:length(sdMin$sdrs), function(i) max( coefVar* abs(sdMin$sdrs[i]),0.3))
  predNegRs  <- intvsS[[6]]$posteriorPredictiveCredibleInterval[1,]*modif - 2 * sapply( 1:length(sdMin$sdrs), function(i) max( coefVar* abs(sdMin$sdrs[i]),0.3))
  predmRs  <- intvsS[[6]]$posteriorPredictiveCredibleInterval[2,]*modif# - 2 * sd(dataX$gpp)
  
  rsPlot<-ggplot()+theme_bw()+
    geom_line(data=dataX,aes(x=timestamp,y=predmRs),colour="purple",size=1)+
    geom_point(data=dataX,aes(x=df2$timestamp, y=rsOb),colour="black",size=2)+
    geom_ribbon(aes(ymin=predNegRs,ymax=predPosRs,x=dataX$timestamp),fill="orange",alpha=0.3)+
    scale_x_datetime(limits=c(as.POSIXct("2015-01-01",tz="GMT"),as.POSIXct("2019-01-01",tz="GMT")))+    
    labs(x="Year",y=expression(paste("Rs [gDM"," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  
  
  dfLeaf<-df2%>%group_by(Year,Month) %>%
    dplyr::summarise(LAI=mean(LAI),t.proj=median(t.proj))
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
    c(dplyr::filter(dfLeaf,Year==2018&Month==7)$Year+dplyr::filter(dfLeaf,Year==2018&Month==7)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==7)$t.proj,429.52)
  )
  )
  names(totC)<-c("year","age","totC")
  
  totN<-data.frame(rbind(
    c(dplyr::filter(dfLeaf,Year==2018&Month==7)$Year+dplyr::filter(dfLeaf,Year==2018&Month==7)$Month/dt,dplyr::filter(dfLeaf,Year==2018&Month==7)$t.proj,14.30)
  )
  )
  names(totN)<-c("year","age","totN")
  
  pLAI<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=LAI))+
    geom_point(data=leaf,aes(x=year,y=lai),shape=16,size=3,colour="red")+
    labs(x="Year",y=expression(paste("L"," ","[",m^-2," ",m^-2,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pStemNo<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=N))+
    geom_point(data=stemNo,aes(x=year,y=stem),shape=16,size=3,colour="red")+
    labs(x="Year",y=expression(paste("N"," ","[",ha^-1,"]",sep="")))+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  pDBH<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=dg))+
    geom_point(data=DBH,aes(x=year,y=dbh),shape=16,size=3,colour="red")+
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
    geom_line(data=df,aes(x=Year+Month/dt,y=totC))+
    geom_point(data=totC,aes(x=year,y=totC),shape=16,size=3,colour="red")+
    labs(x="Year",y="totC [t/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  ptotN<-ggplot()+theme_bw()+
    geom_line(data=df,aes(x=Year+Month/dt,y=totN))+
    geom_point(data=totN,aes(x=year,y=totN),shape=16,size=3,colour="red")+
    labs(x="Year",y="totN [t/ha]")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=14))
  
  return(list(gpp1,swcPlot,NEEPlot,etrans,recoPlot,rsPlot,gppC,nppC, pLAI,pStemNo,pDBH,pWr,ptotC,ptotN))
}

## Function to plot the key output pools and nutritional modifier of the model
plotModel<-function(output.df){
  
  cls<-iprior::gg_colour_hue(2)
  
  flux<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=GPP,colour=cls[1]))+
    geom_line(aes(x=t.proj,y=NPP,colour=cls[2]))+scale_colour_manual(values=cls,name="Fluxes",labels=c("NPP","GPP"))+theme(panel.grid=element_blank())+labs(x="Year",y="Flux")
  
  
  lai<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=LAI,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("LAI"))+theme(panel.grid=element_blank())+labs(x="Year",y="Flux")
  
  
  g<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=G,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("G"))+theme(panel.grid=element_blank())+labs(x="Year",y="G")
  
  
  dg<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=dg,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("DG"))+theme(panel.grid=element_blank())+labs(x="Year",y="Dg")
  
  
  hdom<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=hdom,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("HDOM"))+theme(panel.grid=element_blank())+labs(x="Year",y="H")
  
  
  rad<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=RAD,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("RAD"))+theme(panel.grid=element_blank())+labs(x="Year",y="Rg")
  
  fN<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=fN,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("fNn"))+theme(panel.grid=element_blank())+labs(x="Year",y="fNn")
  
  nav<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=Nav,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("Nav"))+theme(panel.grid=element_blank())+labs(x="Year",y="Nav")
  
  
  totc<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=totC,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("totC"))+theme(panel.grid=element_blank())+labs(x="Year",y="TotC")
  
  
  totn<-ggplot(data=output.df)+theme_bw()+geom_line(aes(x=t.proj,y=totN,colour=cls[1]))+
    scale_colour_manual(values=cls,name="Fluxes",labels=c("totN"))+theme(panel.grid=element_blank())+labs(x="Year",y="TotN")
  
  
  return(list(flux,lai,rad,g,dg,hdom,fN,nav,totc,totn))
}