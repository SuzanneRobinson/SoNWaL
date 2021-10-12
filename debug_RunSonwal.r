

######################
## START THE SCRIPT ##
######################

## Load necessary packages
library(fr3PGDN)
library("tidyverse")
library("lubridate")
library("coda")
library("BayesianTools")
library("miscTools")
library("ggpubr")
## Years of data to use for calibration
startYear = 2015
endYear = 2018
#install.packages("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\fr3PGDN_2.0.tar.gz", repos = NULL, type="source")
timeStep<-"monthly"

#Directory where climate data is stored (default is "data" in SonWal folder)
climDir<-("Data/")


## read in and format climate data
clm_df_full<-data.frame(getClimDatX("monthly",climDir))%>%
  filter(Year<2019)


## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  flxdata_daily <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
  }else
{
  data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
}
###########################
## Initialise Parameters ##
###########################
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                wiltPoint = 0.1, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
                fieldCap =0.24,#Field capacity
                satPoint = 0.3, #field saturation point
                #K_s=0.1, #Soil conductivity
                shared_area=3, #shared area of rooting and non-rooting zone
                V_nr=4, #Volume of non-rooting zone
                maxRootDepth=2.5,
                sigma_zR =0.3, #area/depth explored by 1kg of root biomass
                SWC_nr=10, #SWC of non-rooting zone at time 0
                E_S1 =4.5, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                #K_drain=0.16,
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################
#names of parameters to fit

#exampParams<-read.csv("exampParams.csv")
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
sitka[nm]<-exampParams[nm]



out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_2_T.RDS")

#out<-getSample(out,start=12000,thin=5,numSamples=500)
#codM<-out$chain[[2]][c(1:5000),]
codM<-as.data.frame(out$chain[[4]])
codM<-mergeChains(out$chain)

codM<-miscTools::colMedians(as.data.frame(codM))
codM<-tail(as.data.frame(codM),1)
names(codM)<-nm


#priorSamp<-priorVals$sampler(35000)
#MCMCtrace(getSample(out,coda = T,thin=2,start=5000),wd="C:\\Users\\aaron.morris", post_zm=F,iter=10000,priors = priorSamp)

sitka<-getParms(waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
#sitka$E_S1<-2
#sitka$weather<-clm_df_pine
sitka[nm]<-codM[nm]
#sitka$SWpower0<-round(sitka$SWpower0)
#sitka$E_S2<-0.1
#sitka$E_S1<-2
#sitka$weather[sitka$weather$Year==2015,"Rain"][1]<-filter(clm_df_full,Year==2014)$Rain[1]


#presc<-data.frame(cycle=c(1,1),t=c(20,30),pNr=c(0.1,0.15),thinWl=c(.1,.1),
 #                thinWsbr=c(1,1),thinWr=c(.5,.5),t.nsprouts=c(1,.1))
#sitka$presc<-presc


#sitka$klmax<-0.0001
#for(i in c(1:length(codM))){
#  sitka<-getParms(     wiltPoint = 0.06037,fieldCap =0.28403,satPoint  =0.50171 ,   K_s=13.36049,  V_nr=3.596,
#                       sigma_zR =0.98164 ,  E_S1 =1.51894  , E_S2 =3.09456,shared_area =4.30381,maxRootDepth =2.00764,
#                       K_drain =447.7182 ,  pFS2=1.03867 ,  pFS20 =0.21525  ,   aS    =0.01731  ,nS =3.00313     ,pRx=0.59633,
#                       pRn =0.3305, gammaFx =0.02076, gammaF0 =0.00848 , tgammaF =95.99751, Rttover =0.18559 ,    mF =0.00104   ,  mR =0.52931   ,
#                       mS =0.50793 , SLA0 =5.2504 , SLA1 =8.1942 ,  tSLA =2.0036,  alpha  =0.04856  ,    Y  =0.52186, 
#                       m0=0.05019, MaxCond =0.01001 , LAIgcx =0.67487, CoeffCond =0.09998, BLcond =0.10841,
#                       Nf=0.12252,    Navm =0.00966 ,  Navx =28.0322 , klmax =0.01234,  krmax =0.00548 , komax  =0.00819 ,  
#                       hc  =0.13474  ,  qir  =249.3636,    qil =19.54624,      qh =22.92336,  
#                       qbc  =4.38353,    el =0.82495 ,   er=0.0298, SWconst0=0.30329 ,SWpower0 =2.00351   ,    Qa  =-94.98497 ,   Qb=0.85812,
#                       waterBalanceSubMods=ifelse(args[4]=="T",TRUE,FALSE), timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
#  sitka$weather<-clm_df_full
#sitka[nm[i]]<-codM[nm[i]]
#output<-do.call(fr3PGDN,sitka)
#plot(output$LAI,main=nm[i])
#}
#sitka$pFS20<-0.1
output<-do.call(fr3PGDN,sitka)
ff<-filter(output,Year>2014)
plot(output$LAI)




plot(output$N)
plot(ff$EvapTransp)
tail(output$GPP)
#ff<-filter(output,Year>2014)
plot(ff$fSW)
plot(ff$volSWC_rz)
plot(ff$EvapTransp)
plot(output$totN)
plot(ff$GPP*7.14,col="red")

#ff<-filter(output,t<=100)

sum(tail(ff$YlN,1),tail(ff$ON,1),tail(ff$YrN,1))
sum(tail(ff$YlC,1),tail(ff$OC,1),tail(ff$YrC,1))

mean((ff$YlN)/(ff$Nav))
mean((ff$ON)/(ff$Nav))
mean((ff$YrN)/(ff$Nav))


sum((ff$YlN[4000]),(ff$ON[4000]),(ff$YrN[4000]))
sum((ff$YlC[4000]),(ff$OC[4000]),(ff$YrC[4000]))


output<-do.call(fr3PGDN,sitka)%>%
  filter(Year>=2015)%>%
  group_by(Year,Month)%>%
  dplyr::summarise(mean=mean(GPP*7.14,na.rm=TRUE))
plot(output$mean)


#Plot results
nmc = nrow(out$chain[[1]])
outSample   <- getSample(out,start=nmc/1.1,thin=5)
results<-plotResultsNewMonthly(output,ShortTS=T,out=outSample,numSamps = 25)
ggarrange(results[[1]],results[[2]],results[[3]],results[[5]],results[[4]])
ggarrange(results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]])


ggarrange(results[[1]],results[[2]],results[[3]],results[[5]],results[[4]],results[[6]],results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]])
fileNm<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\"
saveRDS(out,file=paste0(fileNm,timeStep,"_",iters,"_.RDS"))
ggsave(paste0(fileNm,timeStep,"_",iters,".png"),width = 400,
       height = 300,
       units = "mm")



resultsPine<-plotResultsPine(outP)

#cv is ratio of sd to mean, cv * value of meanNEE
#sdmin_NEE               <- cv_EC * abs(data_NEEmean_value) #gives you the standard deviation for the averaged data points (assuming cv 20%)
#data_NEEmean_sd <- sapply( 1:length(sdmin_NEE), function(i) max(sdmin_NEE[i],0.5) )


codM<-as.data.frame(outP$chain[[2]])
codM<-mergeChains(outP$chain)
codM<-miscTools::colMedians(as.data.frame(codM))
codM<-tail(as.data.frame(codM),1)
names(codM)<-nm
pine<-getParmsPine(waterBalanceSubMods=T, timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)
pine[nm]<-codM[nm]
output<-do.call(fr3PGDN,pine)
ff<-filter(output,Year>1996)
plot(ff$GPP)


calc_risk <- function(strtyr,endyr,df,hazval){
  fldf    <- df %>% filter(Year >= strtyr & Year <= endyr)
  lowyrs  <- fldf$Year[(fldf$annual_precip <= hazval)]
  highyrs <- fldf$Year[!(fldf$Year %in% lowyrs)]
  vuln    <- mean(output[weather$Year %in% highyrs,"GPP"]) - mean(output[weather$Year %in% lowyrs,"GPP"])
  haz     <- length(lowyrs)/((endyr - strtyr)+1)
  return(tibble("startYr"= strtyr,
                "endYr"=endyr,
                "vulnerability"= vuln,
                "hazard"  = haz,
                "risk" = vuln*haz))
}


annualPrecip <- weather%>%group_by(Year)%>%summarise(annual_precip=sum(Rain))
hazprecip <- quantile(annualPrecip$annual_precip,0.1)
#' # Expected loss of GPP 
inpt   <- tibble(strtyr=c(1985,1995,2005,2015),endyr=c(1988,1998,2008,2018))
riskdf <- purrr::map2_df(inpt$strtyr,inpt$endyr,calc_risk,df=annualPrecip,hazval=hazprecip)
