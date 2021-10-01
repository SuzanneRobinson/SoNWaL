##################get climate data for calibration################
###needs updating and cleaning when we finalise data, currently just aggregates some daily data
getClimDat2019<-function(timeStep="monthly"){
  library(dplyr)
  library(lubridate)
  
if(Sys.info()[1]=="Windows"){
  clm_df_full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\clm_df_full.csv")
  clm_df_daily<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\weather_day_basfor.csv")
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
}else
{
  clm_df_full<-read.csv("/home/users/aaronm7/3pgData/clm_df_full.csv")
  clm_df_daily<-read.csv("/home/users/aaronm7/3pgData/weather_day.csv")
  data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
  
}
  

  
#  dups<- clm_df_daily[duplicated(clm_df_daily$Date)==T,]
#  for (i in c(1:nrow(dups))){
#dupY<-filter(clm_df_daily,SolarRad==dups[i,7]) 
#medDup<-median(dupY$DOY)  
#print(paste0("medDup= ",medDup))
#print(paste0("dups =", dups[i,2]))
##if(dups[i,2]==medDup) print("yes") else print("no")
#
#}
  
  
#Add date
clm_df_full$date<-as.Date(paste(clm_df_full$Year,"-",clm_df_full$Month,"-01",sep=""))
clm_df_full$week<-week(clm_df_full$date)
clm_df_daily$Date<-as.Date(clm_df_daily$DOY, origin = paste0(clm_df_daily$Year,"-01-01"))

clm_df_daily$dup<-duplicated(clm_df_daily$Date)
clm_df_daily<-clm_df_daily[duplicated(clm_df_daily$Date)==F,]

clm_df_daily$week<-week(clm_df_daily$Date)
clm_df_daily$month<-month(clm_df_daily$Date)
clm_df_daily[which(clm_df_daily$DOY==365&clm_df_daily$week==1),"week"]<-52


modDat<-filter(clm_df_daily,Year>2014)
predDat<-filter(clm_df_daily,Year<=2006)
missingDat<-clm_df_daily%>%filter(Year<2015)
repYear<-clm_df_daily%>%filter(Year==2019)
fillDat<-rbind(do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(4,repYear ,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1),
               do.call("rbind",replicate(2,repYear ,simplify=F)),
                              do.call("rbind",replicate(4, repYear,simplify=F)),
               clm_df_daily%>%filter(Year==2019&DOY==1))

missingDat<-missingDat[!duplicated(missingDat$Date),]

fillDat$Year<-missingDat$Year

 ##use last 5 years to fill in climate gaps?
 #clm_df_dailyYavg<-clm_df_daily%>%filter(Year<2007&Year>2002)%>%
 #  group_by(DOY)%>%
 #  summarise(Tmax=median(Tmax),Tmean=median(Tmean),Tmin=median(Tmin),Rain=median(Rain))
 #
 #clm_df_u<-filter(clm_df_daily,Year>=2007&Year<2015)%>%left_join(clm_df_dailyYavg,by=c("DOY"="DOY"),keep=F)%>%
 #  mutate(Year=Year,Tmax=Tmax.y,Tmin=Tmin.y,Tmean=Tmean.y,Rain=Rain.y)%>%
 #  select(Year,DOY,Tmax,Tmin,Tmean,Rain,SolarRad,VPD,FrostHours,DayIrrig,RH,SWC,Date,week,month)
 
 
 mod1<-lm(SolarRad~Tmax+Tmean+Rain,data=modDat[-c(1275,2103,615),])

 predDat$SolarRad<-predict(mod1,newdata = predDat)

 
 clm_df_daily<-rbind(fillDat,modDat)

clm_df_daily$SolarRad[clm_df_daily$SolarRad<0]<-0


clm_df_daily <- PredictWeatherVariables(weather = clm_df_daily)


clm_df_weekly<-clm_df_daily%>%
  group_by(Year,week)%>%
  summarise(Year=median(Year),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
            Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
          ,FrostDays=mean(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))
  


clm_df_full<-clm_df_daily%>%
  group_by(Year,month)%>%
  summarise(Year=median(Year),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
            Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
            ,FrostDays=mean(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))


clm_df_daily<-clm_df_daily%>%
  group_by(Year,DOY)%>%
  summarise(Year=median(Year),week=median(week),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
            Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
            ,FrostDays=mean(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))


if(timeStep=="monthly") return (clm_df_full)
if(timeStep=="weekly") return (clm_df_weekly)
if(timeStep=="daily") return (clm_df_daily)
#if(timeStep=="pseudo") return (clm_df_daily)

}
################################################################################


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
  
  if(timeStep=="monthly"){
    clm_df_pine<-aggregate(weathdata$tmean_degC~weathdata$mo+weathdata$year,FUN=mean)
    clm_df_pine<-data.frame(Year=clm_df_pine$`weathdata$year`,Month=clm_df_pine$`weathdata$mo`,Tmean=clm_df_pine$`weathdata$tmean_degC`)
    
    clm_df_pine$Tmax         <- aggregate(weathdata$tmax_degC~weathdata$mo+weathdata$year,FUN=max)[,3]
    clm_df_pine$Tmin         <- aggregate(weathdata$tmin_degC~weathdata$mo+weathdata$year,FUN=min)[,3]
    clm_df_pine$Rain      <- aggregate(weathdata$p_mm~weathdata$mo+weathdata$year,FUN=sum)[,3]
    clm_df_pine$SolarRad        <- aggregate(weathdata$rad_Jcm2*((100*100)/1000000)~weathdata$mo+weathdata$year,FUN=mean)[,3] # convert to MJ/m*m/day
    clm_df_pine$FrostDays<-0
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
    clm_df_pine$FrostDays<-0
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
    clm_df_pine$FrostDays<-0
    clm_df_pine$MonthIrrig<-0
    
  }
  
 
  clm_df_pine <- PredictWeatherVariables(weather = clm_df_pine)
  
  clm_df_pine$Month<-month(clm_df_pine$date)
  clm_df_pine<-select(clm_df_pine,c("Year","Week","Month", "Tmax" , "Tmin" ,"Tmean",  "Rain", "SolarRad" ,"FrostDays", "MonthIrrig",      "VPD",    "RH"))
  
  return(clm_df_pine)
  
}


##################get climate data for calibration################
###needs updating and cleaning when we finalise data, currently just aggregates some daily data
getClimDat<-function(timeStep="monthly"){
  library(dplyr)
  library(lubridate)
  
  if(Sys.info()[1]=="Windows"){
    clm_df_full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\clm_df_full.csv")
    clm_df_daily<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\weather_day_basfor.csv")
    data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
    
  }else
  {
    clm_df_full<-read.csv("/home/users/aaronm7/3pgData/clm_df_full.csv")
    clm_df_daily<-read.csv("/home/users/aaronm7/3pgData/weather_day.csv")
    data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
    
    
  }
  
  #Add date
  clm_df_full$date<-as.Date(paste(clm_df_full$Year,"-",clm_df_full$Month,"-01",sep=""))
  clm_df_full$week<-week(clm_df_full$date)
  clm_df_daily$Date<-as.Date(clm_df_daily$DOY, origin = paste0(clm_df_daily$Year,"-01-01"))
  clm_df_daily$week<-week(clm_df_daily$Date)
  clm_df_daily$month<-month(clm_df_daily$Date)
  clm_df_daily[which(clm_df_daily$DOY==365&clm_df_daily$week==1),"week"]<-52
  
  
  modDat<-filter(clm_df_daily,Year>2014)
  predDat<-filter(clm_df_daily,Year<=2006)
  missingDat<-clm_df_daily%>%filter(Year>=2007&Year<2015)
  repYear<-clm_df_daily%>%filter(Year>1997&Year<2006)
  #fillDat<-rbind(do.call("rbind",replicate(4,repYear1 ,simplify=F)),clm_df_daily%>%filter(Year==2019&DOY==1),
  #                             do.call("rbind",replicate(4, repYear,simplify=F)),clm_df_daily%>%filter(Year==2019&DOY==1))
  
  repYear<-repYear[!duplicated(repYear$Date),]
  repYear$Year<-missingDat$Year
  fillDat<-repYear
  ##use last 5 years to fill in climate gaps?
  #clm_df_dailyYavg<-clm_df_daily%>%filter(Year<2007&Year>2002)%>%
  #  group_by(DOY)%>%
  #  summarise(Tmax=median(Tmax),Tmean=median(Tmean),Tmin=median(Tmin),Rain=median(Rain))
  #
  #clm_df_u<-filter(clm_df_daily,Year>=2007&Year<2015)%>%left_join(clm_df_dailyYavg,by=c("DOY"="DOY"),keep=F)%>%
  #  mutate(Year=Year,Tmax=Tmax.y,Tmin=Tmin.y,Tmean=Tmean.y,Rain=Rain.y)%>%
  #  select(Year,DOY,Tmax,Tmin,Tmean,Rain,SolarRad,VPD,FrostHours,DayIrrig,RH,SWC,Date,week,month)
  
  
  mod1<-lm(SolarRad~Tmax+Tmean+Rain,data=modDat[-c(1275,2103,615),])
  
  predDat$SolarRad<-predict(mod1,newdata = predDat)
  
  
  clm_df_daily<-rbind(predDat,fillDat,modDat)
  
  clm_df_daily$SolarRad[clm_df_daily$SolarRad<0]<-0
  
  
  clm_df_daily <- PredictWeatherVariables(weather = clm_df_daily)
  
  
  clm_df_weekly<-clm_df_daily%>%
    group_by(Year,week)%>%
    summarise(Year=median(Year),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=mean(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))
  
  
  
  clm_df_full<-clm_df_daily%>%
    group_by(Year,month)%>%
    summarise(Year=median(Year),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=mean(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))
  
  
  clm_df_daily<-clm_df_daily%>%
    group_by(Year,DOY)%>%
    summarise(Year=median(Year),week=median(week),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=mean(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))
  
  
  if(timeStep=="monthly") return (clm_df_full)
  if(timeStep=="weekly") return (clm_df_weekly)
  if(timeStep=="daily") return (clm_df_daily)
  #if(timeStep=="pseudo") return (clm_df_daily)
  
}
################################################################################

##################get climate data for calibration################
###needs updating and cleaning when we finalise data, currently just aggregates some daily data
getClimDatX<-function(timeStep="monthly",climDir){
  library(dplyr)
  library(lubridate)
  
    clm_df_full<-read.csv(paste0(climDir,"clm_df_full.csv"))
    clm_df_daily<-read.csv(paste0(climDir,"weather_day_basfor.csv"))
    #clm_df_daily2<-read.csv(paste0(climDir,"weather_day.csv"))
   # clm_df_daily2[which(clm_df_daily2$Year>2014 & clm_df_daily2$Year <2019),]<-clm_df_daily[,-13]
    #clm_df_daily<-clm_df_daily2
    
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
  clm_df_daily$FrostHours<-0#ifelse(clm_df_daily$Tmean<=0,1,0)
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
              ,FrostDays=sum(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))
  
  
  clm_df_daily<-clm_df_daily%>%
    group_by(Year,DOY)%>%
    summarise(Year=median(Year),week=median(week),Month=median(month(Date)),Tmax=max(Tmax),Tmin=min(Tmin),
              Tmean=mean(Tmean),Rain=sum(Rain),SolarRad=mean(SolarRad)
              ,FrostDays=sum(FrostHours),MonthIrrig=mean(DayIrrig), VPD=mean(VPD),RH=mean(RH),SWC=mean(SWC/100))
  
  
  if(timeStep=="monthly") return (clm_df_full)
  if(timeStep=="weekly") return (clm_df_weekly)
  if(timeStep=="daily") return (clm_df_daily)
  #if(timeStep=="pseudo") return (clm_df_daily)
  
}
################################################################################
