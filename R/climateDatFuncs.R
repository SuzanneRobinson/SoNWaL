
#' getClmPine get climate data for finnish scots pine using PREDICT database
#' @param timeStep time step to get data for ("weekly", "monthly" or "daily")
#' @return climate data
#' @export
#' 
getClmPine<-function(timeStep="monthly"){
  
  ## 1996-2012
  NinetySixtoEnd <- getData(site=sites[site],dataset="CLIMATE_LOCAL")
  ## 2001-2005
  twentyZeroOnetoZeroFive <- filter(NinetySixtoEnd, year %in% seq(2001,2005))
  ## 2001-2010
  twentyZeroOnetoTen      <- filter(NinetySixtoEnd, year %in% seq(2001,2010))
  ## 1996-2010
  nineteenNinetySixtoTen  <- filter(NinetySixtoEnd, year %in% seq(1996,2010))
  ## combine to give twenty yrs
  twentyYrs                     <- rbind(twentyZeroOnetoZeroFive,nineteenNinetySixtoTen)
  ## 1981 to 1995
  nineteenEightyOnetoNinetyFive <- rbind(twentyZeroOnetoTen, twentyZeroOnetoZeroFive)
  ## 1961 to 2012
  
  weathdata                     <- rbind(twentyYrs, nineteenEightyOnetoNinetyFive, NinetySixtoEnd)
  extendDates                   <- seq(as.Date("1961/1/1"), as.Date(NinetySixtoEnd$date[length(NinetySixtoEnd$date)]), "days")
  extendYrs                     <- as.integer(format(extendDates,format="%Y"))
  weathdata                     <- weathdata %>% mutate(date = extendDates, year = extendYrs)
  
  weathdata$week<-week(weathdata$date)
  weathdata$week[weathdata$week==53]<-52
  
  weathdata$FrostDays<-ifelse(weathdata$tmin_degC<=0,1,0)
  
  if(timeStep=="monthly"){
    clm_df_pine<-aggregate(weathdata$tmean_degC~weathdata$mo+weathdata$year,FUN=mean)
    clm_df_pine<-data.frame(Year=clm_df_pine$`weathdata$year`,Month=clm_df_pine$`weathdata$mo`,Tmean=clm_df_pine$`weathdata$tmean_degC`)
    
    clm_df_pine$Tmax         <- aggregate(weathdata$tmax_degC~weathdata$mo+weathdata$year,FUN=max)[,3]
    clm_df_pine$Tmin         <- aggregate(weathdata$tmin_degC~weathdata$mo+weathdata$year,FUN=min)[,3]
    clm_df_pine$Rain      <- aggregate(weathdata$p_mm~weathdata$mo+weathdata$year,FUN=sum)[,3]
    clm_df_pine$SolarRad        <- aggregate(weathdata$rad_Jcm2*((100*100)/1000000)~weathdata$mo+weathdata$year,FUN=mean)[,3] # convert to MJ/m*m/day
    clm_df_pine$FrostDays <- aggregate(weathdata$FrostDays~weathdata$mo+weathdata$year,FUN=sum)[,3]
    clm_df_pine$MonthIrrig<-0
    clm_df_pine$RH<-aggregate(weathdata$relhum_percent~weathdata$mo+weathdata$year,FUN=mean)[,3]
    clm_df_pine$date<-as.Date(paste0(clm_df_pine$Year,"-",clm_df_pine$Month,"-01"))
    
    
  }
  
  if(timeStep=="weekly"){
    clm_df_pine<-aggregate(weathdata$tmean_degC~weathdata$week+weathdata$year,FUN=mean)
    clm_df_pine<-data.frame(Year=clm_df_pine$`weathdata$year`,Week=clm_df_pine$`weathdata$week`,Tmean=clm_df_pine$`weathdata$tmean_degC`)
    
    clm_df_pine$Tmax         <- aggregate(weathdata$tmax_degC~weathdata$week+weathdata$year,FUN=max)[,3]
    clm_df_pine$Tmin         <- aggregate(weathdata$tmin_degC~weathdata$week+weathdata$year,FUN=min)[,3]
    clm_df_pine$Rain      <- aggregate(weathdata$p_mm~weathdata$week+weathdata$year,FUN=sum)[,3]
    clm_df_pine$SolarRad        <- aggregate(weathdata$rad_Jcm2*((100*100)/1000000)~weathdata$week+weathdata$year,FUN=mean)[,3] # convert to MJ/m*m/day
    clm_df_pine$FrostDays<- aggregate(weathdata$FrostDays~weathdata$week+weathdata$year,FUN=sum)[,3]
    clm_df_pine$MonthIrrig<-0
    clm_df_pine$RH<-aggregate(weathdata$relhum_percent~weathdata$week+weathdata$year,FUN=mean)[,3]
    clm_df_pine$date<-as.Date(paste(clm_df_pine$Year, clm_df_pine$Week, 1, sep="-"), "%Y-%U-%u")
    
    clm_df_pine$Month<-month(clm_df_pine$date)
    clm_df_pine<-clm_df_pine[,c(1,11,3:10,2)]
    
  }
  
  if(timeStep=="daily"){
    clm_df_pine<-aggregate(weathdata$tmean_degC~weathdata$week+weathdata$year,FUN=mean)
    clm_df_pine<-data.frame(Year=clm_df_pine$`weathdata$year`,Month=clm_df_pine$`weathdata$week`,Tmean=clm_df_pine$`weathdata$tmean_degC`)
    
    clm_df_pine$Tmax         <- aggregate(weathdata$tmax_degC~weathdata$week+weathdata$year,FUN=max)[,3]
    clm_df_pine$Tmin         <- aggregate(weathdata$tmin_degC~weathdata$week+weathdata$year,FUN=min)[,3]
    clm_df_pine$Rain      <- aggregate(weathdata$p_mm~weathdata$week+weathdata$year,FUN=sum)[,3]
    clm_df_pine$SolarRad        <- aggregate(weathdata$rad_Jcm2*((100*100)/1000000)~weathdata$week+weathdata$year,FUN=mean)[,3] # convert to MJ/m*m/day
    clm_df_pine$FrostDays<- aggregate(weathdata$FrostDays~weathdata$week+weathdata$year,FUN=sum)[,3]
    clm_df_pine$MonthIrrig<-0
    
  }
  
 
  clm_df_pine <- PredictWeatherVariables(weather = clm_df_pine)
  
  clm_df_pine$Month<-month(clm_df_pine$date)
  
 
  
  clm_df_pine<-dplyr::select(clm_df_pine,c("Year","Week","Month", "Tmax" , "Tmin" ,"Tmean",  "Rain", "SolarRad" ,"FrostDays", "MonthIrrig",      "VPD",    "RH"))
  
  return(clm_df_pine)
  
}


#' getClimDatX get climate data for Harwood - sitka spruce calibrations
#' @param timeStep time step to get data for ("weekly", "monthly" or "daily")
#' @param climDir directory where climate data is stored
#' @return climate data
#' @export
getClimDatX<-function(timeStep="monthly",climDir){
  library(dplyr)
  library(lubridate)
  

  clm_df_full<-read.csv(paste0(climDir,"clm_df_full.csv"))
  clm_df_daily<-read.csv(paste0(climDir,"weather_day_basfor.csv"))
  
  #Add date
  
  clm_df_daily<- rbind(do.call("rbind", replicate(22, (clm_df_daily[1:730,]), simplify = FALSE)),clm_df_daily[732:1461,])
  
  clm_df_daily$Year<-rep(1973:2018,each=365)
  
  clm_df_full$date<-as.Date(paste(clm_df_full$Year,"-",clm_df_full$Month,"-01",sep=""))
  clm_df_full$week<-week(clm_df_full$date)
  clm_df_daily$Date<-as.Date(clm_df_daily$DOY, origin = paste0(clm_df_daily$Year,"-01-01"))
  clm_df_daily$week<-week(clm_df_daily$Date)
  clm_df_daily$month<-month(clm_df_daily$Date)
  clm_df_daily[which(clm_df_daily$DOY==365&clm_df_daily$week==1),"week"]<-52
  
  modDat<-filter(clm_df_daily,Year>2014)
  predDat<-filter(clm_df_daily,Year<=2014)
  
  
  
  mod1<-lm(SolarRad~Tmax+Tmean+Rain,data=modDat[-c(1275,2103,615),])
  
  # predDat$SolarRad<-predict(mod1,newdata = predDat)
  
  
  clm_df_daily<-rbind(predDat,modDat)
  
  clm_df_daily$SolarRad[clm_df_daily$SolarRad<0]<-0
  
  
  clm_df_daily <- PredictWeatherVariables(weather = clm_df_daily)
  clm_df_daily$FrostHours<-ifelse(clm_df_daily$Tmin<=0,1,0)
  clm_df_daily$rainDays<-ifelse(clm_df_daily$Rain>0,1,0)
  
  clm_df_weekly<-clm_df_daily%>%
    group_by(Year,week)%>%
    summarise(Year=median(Year),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=sum(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100),rainDays=sum(rainDays))
  
  
  
  clm_df_full<-clm_df_daily%>%
    group_by(Year,month)%>%
    summarise(Year=median(Year),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=sum(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100),rainDays=sum(rainDays))
  
  
  clm_df_daily<-clm_df_daily%>%
    group_by(Year,DOY)%>%
    summarise(Year=median(Year),week=median(week),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=sum(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100),rainDays=sum(rainDays))
  
  
  if(timeStep=="monthly") return (clm_df_full)
  if(timeStep=="weekly") return (clm_df_weekly)
  if(timeStep=="daily") return (clm_df_daily)
  #if(timeStep=="pseudo") return (clm_df_daily)
  
}
