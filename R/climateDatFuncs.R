##################get climate data for calibration################
###needs updating and cleaning when we finalise data, currently just aggregates some daily data
getClimDat<-function(timeStep="monthly"){
  library(dplyr)
  library(lubridate)
  
if(Sys.info()[1]=="Windows"){
  clm_df_full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\clm_df_full.csv")
  clm_df_daily<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\weather_day.csv")
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

clm_df_dailyX<-clm_df_daily[!duplicated(clm_df_daily[c(1,12)]),]

weeklyRfall<-aggregate(clm_df_daily$Rain~clm_df_daily$week+clm_df_daily$Year,FUN=sum)
names(weeklyRfall)<-c("week","Year","Rain")

weeklyRfall$SolarRad<-aggregate(clm_df_daily$SolarRad~clm_df_daily$week+clm_df_daily$Year,FUN=mean)[,3]
weeklyRfall$Tmax<-aggregate(clm_df_daily$Tmax~clm_df_daily$week+clm_df_daily$Year,FUN=max)[,3]
weeklyRfall$Tmin<-aggregate(clm_df_daily$Tmin~clm_df_daily$week+clm_df_daily$Year,FUN=min)[,3]
weeklyRfall$Tmean<-aggregate(clm_df_daily$Tmean~clm_df_daily$week+clm_df_daily$Year,FUN=mean)[,3]

weeklyRfall$Month<-month(clm_df_dailyX$Date)
weeklyRfall$FrostDays<-(aggregate(clm_df_daily$FrostHours~clm_df_daily$week+clm_df_daily$Year,FUN=sum)[,3])/24
weeklyRfall$MonthIrrig<-0

#Split into weekly data
clm_df_fullX<-NULL
clm_df<-clm_df_full
for(i in c(1:nrow(clm_df_full))){
  reps<-ifelse(clm_df_full$Month[i]!=12,clm_df_full[i+1,"week"]-clm_df_full[i,"week"],4)
  clm_df[i,]$Rain<-clm_df[i,]$Rain/reps
  clm_df[i,]$MonthIrrig<-clm_df[i,]$MonthIrrig/reps
  
  clm_df_fullX<-rbind(clm_df_fullX,do.call("rbind",(replicate(reps, clm_df[i,], simplify = FALSE))))
}
clm_df_fullX$week<-c(1:51)

weeklyRfall<-weeklyRfall[,c(names(clm_df_fullX[,c(1:9)]))]


if(timeStep=="monthly") return (clm_df_full)
if(timeStep=="weekly") return (weeklyRfall)

}
################################################################################

