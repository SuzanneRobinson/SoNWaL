
observedVals<-function(timeStep,data,sY=2015,eY=2018){

  data<-filter(data,year>=sY&year<=eY)
  data<-data[-1,]
  if(timeStep =="weekly"){
    data<-data%>%
      mutate(grp=week(as.Date(data$yday, origin = paste0(data$year,"-01-01"))))
  }
  
  if(timeStep =="monthly"){
    data<-data%>%
      mutate(grp=month(as.Date(data$yday, origin = paste0(data$year,"-01-01"))))

  }
  



  #add missing swc values, can cause problems otherwise
  mod1<-lm(swc~npp+nee+reco,data=data)

  #remove some negative values
  data$et[data$et<0]<-0
  
  data<-data %>% 
    mutate(swc = coalesce(swc,predict(mod1,newdata = data)))

  sdMin<-data%>% group_by(year,grp) %>%
    dplyr::summarise(sdgpp=mean(gpp),sdnpp=mean(npp),sdnee=mean(nee),sdreco=mean(reco),sdrs=mean(rs),sdet=mean(et),sdswc=mean(swc))
  
  sdAlphaAnn<-data%>% group_by(year) %>%
    dplyr::summarise(sdAlphaAnn=mean(reco/rs))
  
  
    observed <- c(sdAlphaAnn$sdAlphaAnn,
                  pull(data%>% 
                         group_by(year,grp) %>%
                         dplyr::summarise(gpp=mean(gpp))%>%
                         select(gpp)),                ## GPP
                #  pull(data%>%
                #         group_by(year,grp) %>%
                #         dplyr::summarise(npp=mean(npp))%>%
                #         select(npp)),                ## NPP
                  pull(data%>%
                         group_by(year,grp) %>%
                         dplyr::summarise(nee=mean(nee))%>%
                         select(nee)),                ## NEE
                #  pull(data%>%
                #         group_by(year,grp) %>%
                #         dplyr::summarise(reco=mean(reco))%>%
                #         select(reco)),               ## Reco
                #  pull(data%>%
                #         group_by(year,grp) %>%
                #         dplyr::summarise(rs=mean(rs))%>%
                #         select(rs)),                 ## Rs
                  pull(data%>%
                         group_by(year,grp) %>%
                         dplyr::summarise(et=mean(et))%>%
                         select(et)),                 ## Etransp
                  #  data$gs[2:nrow(data)],   ## CanCond
                  5.7,5.56,                ## LAI
                  1348,                    ## N - fairly well known
                  24.1,                    ## dg
                  #  4.88,                    ## Wr
                  # 0.53,                    ## difRoots
                  (214.76),                    ## totC, see jarvis_total_soil.ods
                  (7.15),                     ## totN, 40 C:N ratio
                  pull(data%>%
                         group_by(year,grp) %>%
                         dplyr::summarise(swc=mean(swc))%>%
                         select(swc))                ## SWC

    )
    
    #observed<-data.frame(observed)
   # observed$lab<-c(rep("alphaAn",4),rep("gpp",212),rep("nee",212),rep("et",212),"LAI","LAI","N","dg","totC","totN",rep("swc",212))
      
    coefVar1=0.1
    coefVar2=0.3
    coefVar3=0.2
    
    dev <- c(sapply( 1:length(sdAlphaAnn$sdAlphaAnn), function(i) max( coefVar2* abs(sdAlphaAnn$sdAlphaAnn[i]),0.01) ),
             sapply( 1:length(sdMin$sdgpp), function(i) max( coefVar3* abs(sdMin$sdgpp[i]),0.01) ),
             # sapply( 1:length(sdMin$sdnpp), function(i) max( 0.05* abs(sdMin$sdnpp[i]),0.05) ),
             sapply( 1:length(sdMin$sdnee), function(i) max( coefVar3* abs(sdMin$sdnee[i]),0.01) ),
             # sapply( 1:length(sdMin$sdreco), function(i) max( coefVar1* abs(sdMin$sdreco[i]),0.1) ),
             # sapply( 1:length(sdMin$sdrs), function(i) max( coefVar2* abs(sdMin$sdrs[i]),0.01) ),
             sapply( 1:length(sdMin$sdet), function(i) max( coefVar2* abs(sdMin$sdet[i]),0.01) ),
             # rep(0.5,(nrow(dplyr::filter(data,year>=startYear&year<=endYear))-1)),
             5.7*coefVar1,5.56*coefVar1,
             1348*coefVar1,
             24.1*coefVar1,
             #  2,
             #  1,
             214.76*0.5,
             7.15*0.5,
             sapply( 1:length(sdMin$sdswc), function(i) max( coefVar3* abs(sdMin$sdswc[i]),0.01) )
             
    )
 
 
 return(list(observed,dev))

}

observedValsPine<-function(timeStep,fluxDat,SWCData){
  convrtUnit=(12.011 * 24 * 60 * 60)/1000000 # convert to grams per day of C
  if(timeStep=="monthly"){
GPP<-aggregate(fluxDat$gppDtCutRef_umolCO2m2s1*convrtUnit~fluxDat$mo+fluxDat$year,FUN=mean)
names(GPP)<-c("month","year","GPP")
GPP<-filter(GPP,year!=2007)
NEE<-aggregate(fluxDat$neeCutRef_umolCO2m2s1*convrtUnit~fluxDat$mo+fluxDat$year,FUN=mean)
names(NEE)<-c("month","year","NEE")
NEE<-filter(NEE,year!=2007)
reco<-aggregate(fluxDat$recoDtCutRef_umolCO2m2s1*convrtUnit~fluxDat$mo+fluxDat$year,FUN=mean)
names(reco)<-c("month","year","reco")
reco<-filter(reco,year!=2007)
swc<-aggregate(SWCData$swc~SWCData$month +SWCData$year,FUN=mean)
names(swc)<-c("month","year","swc")
swc<-filter(swc,year!=2007)
  }else{
    fluxDat$week<-week(fluxDat$date)
    fluxDat$week[fluxDat$week==53]<-52
    SWCData$week[SWCData$week==53]<-52
    GPP<-aggregate(fluxDat$gppDtCutRef_umolCO2m2s1*convrtUnit~fluxDat$week+fluxDat$year,FUN=mean)
    names(GPP)<-c("week","year","GPP")
    GPP<-filter(GPP,year!=2007)
    NEE<-aggregate(fluxDat$neeCutRef_umolCO2m2s1*convrtUnit~fluxDat$week+fluxDat$year,FUN=mean)
    names(NEE)<-c("week","year","NEE")
    NEE<-filter(NEE,year!=2007)
    reco<-aggregate(fluxDat$recoDtCutRef_umolCO2m2s1*convrtUnit~fluxDat$week+fluxDat$year,FUN=mean)
    names(reco)<-c("week","year","reco")
    reco<-filter(reco,year!=2007)
    swc<-aggregate(SWCData$swc~SWCData$week +SWCData$year,FUN=mean)
    names(swc)<-c("month","year","swc")
    swc<-filter(swc,year!=2007)
  }

#get stand data
standDat <- getData(site=sites[site],dataset="STAND")%>%
  filter(species_id=="pisy")
#get soil data
soildata <-  getData( sites[site], dataset="SOIL")

#convert soil mass data and variable percentage to tons per hectare
massExt<-function(soildata,varMass){
  sProfM<- ((soildata$lowerDepth_cm - soildata$upperDepth_cm))*0.01
  t_ha<-10000*sProfM*soildata$density_gcm3*varMass
  return(t_ha)
}

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
totC<-sum(massExt(soildata,soildata$c_percent),na.rm=TRUE)
totN<-sum(massExt(soildata,soildata$n_percent),na.rm=TRUE)

observedPine <- c(GPP$GPP,                ## GPP - monthly avg
                  NEE$NEE,                ## NPP - monthly avg
                  reco$reco,              ## NEE - monthly avg
                  LAI[,2],                    ##LAI -yearly average 
                  dbh[,2],                    ## DBH
                  totC,                   ## totC, 1995-1996
                  totN,                    ## totN, 1995-1996
                  treeDen$treeDen, #tree density
                  swc$swc
)


startYear<-1996
endYear<-2014
coefVar=0.2
devPine <- c(             sapply( 1:length(GPP$GPP), function(i) max( coefVar* abs(GPP$GPP[i]),0.001) ),
                        sapply( 1:length(NEE$NEE), function(i) max( coefVar* abs(NEE$NEE[i]),0.001) ),
                        sapply( 1:length(reco$reco), function(i) max( coefVar* abs(reco$reco[i]),0.001) ),
                        sapply( 1:length(LAI$LAI), function(i) max( 0.1* abs(LAI$LAI[i]),0.001) ),
                        sapply( 1:length(dbh$dbh), function(i) max( 0.1* abs(dbh$dbh[i]),0.001) ),
                        20,
             5,
             sapply( 1:length(treeDen$treeDen), function(i) max( 0.1* abs(treeDen$treeDen[i]),0.001) ),
             sapply( 1:length(swc$swc), function(i) max( coefVar* abs(swc$swc[i]),0.001) ))


return(list(observedPine,devPine))

}


getSWC<-function(){
  #getSWC
  ll <- getData(dataset="SOILTS",site="hyytiala")
  ll$date <- lubridate::as_datetime(ll$date)
  ll$yday <- floor_date(ll$date, "day")
  ld <- ll %>% group_by(yday) %>% summarize(swcd1= mean(swcFMDS1_degC),swcd2= mean(swcFMDS2_degC),swcd3= mean(swcFMDS3_degC),swcd4= mean(swcFMDS4_degC),swcd5= mean(swcFMDS5_degC))
  ## soil levels weighting
  ls <- getData(dataset="SOIL",site="hyytiala")
  layerThickcm <- ls$lowerDepth_cm-ls$upperDepth_cm
  layerThickcm
  soilThickcm = sum(layerThickcm,na.rm=TRUE)
  soilThickcm
  fracLayer = layerThickcm/soilThickcm
  fracLayer
  #swc = ld$swcd2 * fracLayer[2] +
  swc <- ld %>% mutate(swc = ((ld$swcd2 * fracLayer[2]) +
                                (ld$swcd3 * fracLayer[3]) +
                                (ld$swcd4 * fracLayer[4]) +
                                (ld$swcd5 * fracLayer[5]))/100.)
  
  swc2 <- swc[367:length(swc$yday),]
  SWCData <- swc2 %>% mutate(year=year(yday), doy=yday(yday),date=as.Date(yday(yday), origin = paste0(year(yday),"-01-01"))) %>% select(date,year,doy,swc) %>%
    mutate(month=month(date),week=week(date))

  return(SWCData)
  
}