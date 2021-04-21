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

setwd("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PROFOUND_data\\ProfoundData")
setDB(db_name ="ProfoundData.sqlite")
timeStep="weekly"

clm_df_full<-getClimDat(timeStep)
sitka<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
sites    <- c("peitz","soro" ,"hyytiala","solling_beech","solling_spruce")
fName=paste0("outx_",stringr::str_sub(Sys.time(), 0, -10),stringr::str_sub(Sys.time(), 15, -4),stringr::str_sub(Sys.time(), 18, -1),".RDS")

##hyytiala
site = 3

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

return(clm_df_pine)

}
clm_df_pine<-getClmPine(timeStep)

#get parameters
pine<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
pine$weather<-clm_df_pine
  
#management
#15.5	4230	0	0	0
#20.5	2120	0.554071739	0.561434582	0.573399222
#25.5	1275	0.55064136	0.557545216	0.568770838

#1976 0.4 removed
#1991 0.3 removed
#1995 0.007201646 removed

presc<-data.frame(cycle=c(1,1,1),t=c(15,25,35),pNr=c(0.4,0.3,0.075),thinWl=c(0.4,0.3,0.075),
                  thinWsbr=c(0.4,0.3,0.075),thinWr=c(0.4,0.3,0.075),t.nsprouts=c(1,1,1))


pine$presc<-presc
###Get observed flux data
fluxDat <- getData(site=sites[site],dataset="FLUX")
GPP<-aggregate(fluxDat$gppDtCutRef_umolCO2m2s1~fluxDat$mo+fluxDat$year,FUN=mean)
names(GPP)<-c("month","year","GPP")
GPP<-filter(GPP,year!=2007)
NEE<-aggregate(fluxDat$neeCutRef_umolCO2m2s1~fluxDat$mo+fluxDat$year,FUN=mean)
names(NEE)<-c("month","year","NEE")
NEE<-filter(NEE,year!=2007)
reco<-aggregate(fluxDat$recoDtCutRef_umolCO2m2s1~fluxDat$mo+fluxDat$year,FUN=mean)
names(reco)<-c("month","year","reco")
reco<-filter(reco,year!=2007)


#get stand data
standDat <- getData(site=sites[site],dataset="STAND")
#get soil data
soildata <-  getData( sites[site], dataset="SOIL")

#get field capacity and wilting point for soil profile
#FC <- sum(soildata$fcapv_percent[2:5] * (soildata$lowerDepth_cm[2:5] - soildata$upperDepth_cm[2:5])/sum(soildata$lowerDepth_cm[2:5] - soildata$upperDepth_cm[2:5]))/100
#WP <- sum(soildata$wiltpv_percent[2:5] * (soildata$lowerDepth_cm[2:5] - soildata$upperDepth_cm[2:5])/sum(soildata$lowerDepth_cm[2:5] - soildata$upperDepth_cm[2:5]))/100

#convert soil mass data and variable percentage to tons per hectare
massExt<-function(soildata,varMass){
sProfM<- ((soildata$lowerDepth_cm - soildata$upperDepth_cm))*0.01
t_ha<-10000*sProfM*soildata$density_gcm3*varMass
return(t_ha)
}

#mean LAI - could not find exact dates so using yearly average
LAI<-aggregate(standDat$lai~standDat$year,FUN=mean)
names(LAI)<-c("year","LAI")

#DBH
dbh<-aggregate(standDat$dbhBA_cm~standDat$year,FUN=mean)
names(LAI)<-c("year","dbh")

#carbon and nitrogen measured at differing soil depths - currently taking average
totC<-sum(massExt(soildata,soildata$c_percent),na.rm=TRUE)
totN<-sum(massExt(soildata,soildata$n_percent),na.rm=TRUE)

observedPine <- c(GPP$GPP,                ## GPP - monthly avg
              NEE$NEE,                ## NPP - monthly avg
              reco$reco,              ## NEE - monthly avg
              LAI[,2],                    ##LAI -yearly average 
              dbh[,2],                    ## DBH
              totC,                   ## totC, 1995-1996
              totN                    ## totN, 1995-1996
)


startYear<-1996
endYear<-2014

devPine <- c(rep(.3,nrow(dplyr::filter(GPP,year>=startYear&year<=endYear&year!=2007))),
         rep(.3,nrow(dplyr::filter(NEE,year>=startYear&year<=endYear&year!=2007))),
         rep(.3,nrow(dplyr::filter(reco,year>=startYear&year<=endYear&year!=2007))),
         rep(.1,nrow(LAI)),
         rep(.3,nrow(dbh)),
         20,
         5
         
)


## Extract simulated data for use in likelihood function
#2007 is removed as a year due to some odd observed data values, basically there's no GPP or NEE etc.
sampleOutputPine<-function(df,sY,eY){
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Month,Year)%>%
           summarise(sum=sum(GPP,na.rm=TRUE))%>%
           select(sum)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Month,Year)%>%
           summarise(sum=sum(NEE,na.rm=TRUE))%>%
           select(sum)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Month,Year)%>%
           summarise(sum=sum(Reco,na.rm=TRUE))%>%
           select(sum)),
   
      pull(filter(df,Year>=1995&Year<=2011)%>%
         group_by(Year)%>%
         summarise(mean=mean(LAI,na.rm=TRUE))%>%
        select(mean)),
      
      pull(filter(df,Year>=1995&Year<=2011)%>%
        group_by(Year)%>%
        summarise(mean=mean(dg,na.rm=TRUE))%>%
        select(mean)),
       
      pull(filter(df,Year==1996)%>%
        group_by(Year)%>%
        summarise(mean=mean(totC,na.rm=TRUE))%>%
        select(mean)),
      
      pull(filter(df,Year==1996)%>%
        group_by(Year)%>%
        summarise(mean=mean(totN,na.rm=TRUE))%>%
        select(mean))
      

  )
  return(m)
}




nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")



##Set priors
priors<-createPriors_pine(pine)

pMaxima<-priors[[1]]
pMinima<-priors[[2]]
pMaxima[[30]]<-0.1
Uprior <- createPrior(lower = pMinima, upper = pMaxima)
## Set observed calibration data and uncertainty
startYear = 1996
endYear = 2014


Uprior <- createTruncatedNormalPrior(mean = priors[[3]], sd=(pMaxima-pMinima)*0.3,
                                     lower = pMinima, upper = pMaxima)



## Likelihood function
pineLL <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPine(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}


##management - not sure if needed just yet
#treedens <- na.trim.ts(standDat$density_treeha/10000)
#myrs     <- standDat$year[!is.na(standDat$density_treeha)]
#nyrs     <- length(myrs)
#
#treedens <- c(0.0875, 0.0755, 0.0684)
#myrs     <- c(1995, 2002, 2010)
#nyrs     <- length(myrs)



iters=250000
#Initiate bayesian setup
BS3PGDN <- createBayesianSetup(likelihood = pineLL, prior = Uprior, names = nm, parallel = 7, catchDuplicates = F )
settings = list(
  iterations = iters,
  ## Z = NULL,
  startValue = 7, # internal chain number - dont use these chains for convergence testing 
  nrChains = 1, # Number of chains
  pSnooker = 0.5,
  burnin = round(iters/100*10), #10% burnin
  ## thin = 1,
  ## f = 2.38,
  ## eps = 0,
  parallel = T,
  ## pGamma1 = 0.1,
  ## eps.mult = 0.2,
  ## eps.add = 0,
  ## consoleUpdates = 100,
  ## zUpdateFrequency = 1,
  ## currentChain = 3,
  ## blockUpdate  = list("none",
  ##                     k = NULL,
  ##                     h = NULL,
  ##                     pSel = NULL,
  ##                     pGroup = NULL,
  ##                     groupStart = 1000,
  ##                     groupIntervall = 1000),
  message = TRUE)

#run calibration with all parameters and priors based on initial hydro model runs
out2 <- runMCMC(bayesianSetup = BS3PGDN, sampler = "DEzs", settings = settings)


#Save output
saveRDS(out2,file=paste0(timeStep,"_pine_",fName))

