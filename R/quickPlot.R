
#'quickPlot Function to plot the model against harwood.
#'@param flxdata_daily harwood flux data
#'@param output model output
#'@grouping whether to plot monthly or weekly outputs
#'@return pine plots
#'@export
quickPlot<-function(flxdata_daily,output,grouping="month"){
  
  if(grouping=="week") output$grouping<-c(1:53) else output$grouping<-output$Month
  output<-filter(output,Year>2014)%>%
    group_by(Year,grouping)%>%
    summarise(GPP=mean(GPP),NPP=mean(NPP),NEE=mean(NEE),Reco=mean(Reco),Rs=mean(Rs),EvapTransp=sum(EvapTransp),volSWC_rz=mean(volSWC_rz))
  
  #conversion factor to gCcm-2 for weekly data
  cf=7.142857
  
  if(grouping=="week") flxdata_daily$grouping<-week(flxdata_daily$timestamp) else flxdata_daily$grouping<-month(flxdata_daily$timestamp)
  
  flxdata_weekly<-flxdata_daily%>%
    filter(year>2014)%>%
    group_by(year,grouping)%>%
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
  
  
  
  ggpubr::ggarrange(gpp1,gpp2,gpp3,gpp4,gpp5, ncol=2, nrow=3)
  
}



diagPlots<-function(out,flxdata_daily,nm,param_scaler=NULL){
  lm_eqn <- function(m){
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  
  
  flxdata_daily$month<-month(flxdata_daily$timestamp)
  cf<-7.14
  flxdata_weekly<-flxdata_daily%>%
    filter(year>2014)%>%
    group_by(year,month)%>%
    summarise(gpp=mean(gpp),npp=mean(npp),nee=mean(nee),reco=mean(reco),rs=mean(rs),et=mean(et),swc=mean(swc),timestamp=median(timestamp))
  
  
  codM<-mergeChains(out$chain)
  codM<-miscTools::colMedians(as.data.frame(codM))
  names(codM)<-nm
  
  sitka<-getParms(waterBalanceSubMods=T, timeStp =  52)
  sitka[nm]<-if(is.null(param_scaler)==T) codM[nm] else codM[nm]*param_scaler 
  
  output<-do.call(SoNWaL,sitka)
  ff<-output%>%
    filter(Year>2014&Year<2100)%>%
    group_by(Year,Month)%>%
    summarise(GPP=mean(GPP),NEE=mean(NEE),EvapTransp=mean(EvapTransp),volSWC_rz=mean(volSWC_rz))
  
  flxdata_weekly$simGPP<-ff$GPP*7.14
  flxdata_weekly$simET<-ff$EvapTransp/7
  flxdata_weekly$simSWC<-ff$volSWC_rz
  flxdata_weekly$simNEE<-ff$NEE*7.14
  
  g1<-ggplot(data=flxdata_weekly,aes(x=as.numeric(simGPP), y=as.numeric(gpp)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(0,15)+
    ylim(0,15)+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    ggtitle(expression(paste("GPP [gC"," ",cm^-2,day^-1,"]",sep="")))+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  g2<-ggplot(data=flxdata_weekly,aes(x=as.numeric(simET), y=as.numeric(et)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(0,3.5)+
    ylim(0,3.5)+
    ggtitle(expression(paste("Evapotranspiration [mm"," ",day^-1,"]",sep="")))+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  g3<-ggplot(data=flxdata_weekly,aes(x=as.numeric(simSWC), y=as.numeric(swc)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(0.22,0.32)+
    ylim(0.22,0.32)+
    ggtitle(expression(paste("SWC",sep="")))+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  
  g4<-ggplot(data=flxdata_weekly,aes(x=as.numeric(simNEE), y=as.numeric(nee)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(-6,2)+
    ylim(-6,2)+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    ggtitle(expression(paste("NEE [gC"," ",cm^-2,day^-1,"]",sep="")))+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  mod1<-lm(simGPP~gpp,data=flxdata_weekly)
  mod2<-lm(simET~et,data=flxdata_weekly)
  mod3<-lm(simSWC~swc,data=flxdata_weekly)
  mod4<-lm(simNEE~nee,data=flxdata_weekly)
  
  
  g1 <- g1 + geom_text(x = 10.8, y = 1, label = lm_eqn(mod1), parse = TRUE,size=3)
  g2 <- g2 + geom_text(x = 2.475, y = 0.25, label = lm_eqn(mod2), parse = TRUE,size=3)
  g3 <- g3 + geom_text(x = 0.288, y = 0.225, label = lm_eqn(mod3), parse = TRUE,size=3)
  g4 <- g4 + geom_text(x = -0.25, y = -5.5, label = lm_eqn(mod4), parse = TRUE,size=3)
  
  return(list(g1,g2,g3,g4))
  
}




diagPlotsPine<-function(outP,nm){
  setDB(db_name ="C:/Users/aaron.morris/OneDrive - Forest Research/Documents/Projects/PRAFOR/models/PROFOUND_data/ProfoundData/ProfoundData.sqlite")
  
  lm_eqn <- function(m){
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }
  
  
  fluxDat <- getData(site=sites[site],dataset="FLUX")
  fluxDat$timestamp<-as.POSIXct(fluxDat$date)
  fluxDat$month<-month(fluxDat$timestamp)
  swcDat<-getSWC()
  swcDat$timestamp<-as.POSIXct(swcDat$date)
  
  convrtUnit=(12.011 * 24 * 60 * 60)/1000000 # convert to grams per day of C
  
  flxdata<-fluxDat%>%
    filter(year>1995)%>%
    group_by(year,month)%>%
    summarise(gpp=mean(gppDtCutRef_umolCO2m2s1*convrtUnit),nee=mean(neeCutRef_umolCO2m2s1*convrtUnit),timestamp=median(timestamp))
  
  swcDat<-swcDat%>%
    filter(year>1995)%>%
    group_by(year,month)%>%
    summarise(swc=mean(swc),timestamp=median(timestamp))
  
  flxdata<-merge(flxdata,swcDat,by.x=c("year","month"),by.y = c("year","month"),all=T)
  flxdata<-flxdata[with(flxdata, order(year, month)), ]
  
  codM<-mergeChains(outP$chain)
  codM<-miscTools::colMedians(as.data.frame(codM))
  names(codM)<-nm
  
  pine<-getParmsPine(waterBalanceSubMods=T, timeStp =  52)
  pine[nm]<-codM[nm]
  
  output<-do.call(SoNWaL,pine)
  ff<-output%>%
    filter(Year>1995&Year<2100)%>%
    group_by(Year,Month)%>%
    summarise(GPP=mean(GPP),NEE=mean(NEE),EvapTransp=mean(EvapTransp),volSWC_rz=mean(volSWC_rz))
  
  flxdata$simGPP<-ff$GPP*7.14
  flxdata$simET<-ff$EvapTransp/7
  flxdata$simSWC<-ff$volSWC_rz
  flxdata$simNEE<-ff$NEE*7.14
  
  flxdata<-filter(flxdata,year!=2007)
  
  g1<-ggplot(data=flxdata,aes(x=as.numeric(simGPP), y=as.numeric(gpp)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(0,11)+
    ylim(0,11)+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    ggtitle(expression(paste("GPP [gC"," ",cm^-2,day^-1,"]",sep="")))+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  
  
  g3<-ggplot(data=flxdata,aes(x=as.numeric(simSWC), y=as.numeric(swc)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(0.18,0.45)+
    ylim(0.18,0.45)+
    ggtitle(expression(paste("SWC"," [",day^-1,"]",sep="")))+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  
  g4<-ggplot(data=flxdata,aes(x=as.numeric(simNEE), y=as.numeric(nee)))+
    theme_bw()+
    geom_point(colour="red",size=2)+
    geom_smooth(method='lm',col="black")+
    ylab("Observed")+
    xlab("Simulated")+
    xlim(-5.5,1)+
    ylim(-5.5,1)+
    geom_abline(intercept =0 , slope = 1, lty="dashed")+
    ggtitle(expression(paste("NEE [gC"," ",cm^-2,day^-1,"]",sep="")))+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10))
  
  mod1<-lm(simGPP~gpp,data=flxdata)
  mod3<-lm(simSWC~swc,data=flxdata)
  mod4<-lm(simNEE~nee,data=flxdata)
  
  
  g1 <- g1 + geom_text(x = 9, y = 1, label = lm_eqn(mod1), parse = TRUE,size=3)
  g3 <- g3 + geom_text(x = 0.4, y = 0.2, label = lm_eqn(mod3), parse = TRUE,size=3)
  g4 <- g4 + geom_text(x = -0.25, y = -5.3, label = lm_eqn(mod4), parse = TRUE,size=3)
  
  ggpubr::ggarrange(g1,g3,g4)
  
}


