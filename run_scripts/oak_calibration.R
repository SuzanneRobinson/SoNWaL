# oak calibration
library(SoNWaL)
library(BayesianTools)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)

#read in arguments from batch file
#args=(commandArgs(TRUE))
args<-c("weekly","weekly_1_","1","T")

print(args)
timeStep=args[1]
chainID=args[2]
chainNum=args[3]
#print arguments to log for easy reference
print(timeStep)
print(chainID)
print("oak")

#Time step to run SonWal with
timeStep<-"weekly"
args[1]<-timeStep


#ah_flux<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\AH_flux\\Flux data_GapFill\\Flux data_GapFill\\AliceHolt_partitioned_data_99_19")%>%
#  mutate(NEE = as.numeric(NEE), date = as.Date(DateTime), Rg = as.numeric(Rg),Tair = as.numeric(Tair),Tsoil = as.numeric(Tsoil),  rH = as.numeric(rH),
#         NEE_uStar_f = as.numeric(NEE_uStar_f), GPP_uStar_f = as.numeric(GPP_uStar_f), LE_uStar_f = as.numeric(LE_uStar_f)) %>%
#  mutate(NEE = ifelse(NEE < -500, NA, NEE), Month = month(date), Year = year(date), Week = week(date)) %>%
#  group_by(Year,Week) %>%
#  summarise(across(NEE:LE_uStar_f, ~ mean(.x, na.rm=T)), Month = first(Month), Date = first(date))
#
#convFac<-1800*12*1e-6


ah_flux<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\AH_flux\\AliceHolt_daily_EC_data.csv")[-1,]
ah_menst<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\AH_flux\\straits_mensuration_summary_2009_2020.csv")
# conversion factor 
co2_to_c<-12/44

# manip
ah_flux<-ah_flux%>%
  mutate(month = month(date), year = year(date))%>%
  group_by(year, month) %>%
  summarise(across(NEE:H, ~ mean(as.numeric(.x), na.rm=T))) %>%
  mutate(NEE = NEE * co2_to_c, GPP = GPP * co2_to_c, Reco = Reco * co2_to_c) %>%
  filter(year < 2018)


ah_menst <- ah_menst %>%
  filter(species == "OK") %>% 
  group_by(year) %>%
  summarise (dbh_cm = mean(as.numeric(dbh), na.rm = T), dbh_sd = sd(as.numeric(dbh), na.rm = T), height = mean(as.numeric(total.height), na.rm = T)) %>%
  filter(year<2019)

coefVar = 0.25
obs<- as.vector(c(ah_flux$GPP, ah_flux$NEE, ah_flux$Reco, ah_menst$dbh_cm))
dev<- c(coefVar * abs(obs[1:(length(obs)-2)]), ah_menst$dbh_sd)

climDir<-("Data/")

## read in and format climate data
clm_df_full<-data.frame(getClimDatX(timeStep,climDir))%>%
  filter(Year<2019)

###########################
## Initialise Parameters ##
###########################
oak<-getParmsOak(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
)
#######################################################


baseParms<-getIndvClm(scape, "AH", soil_dat, simDatLoc, oak)
weather<-baseParms$weather
weather$Year<-weather$Year-30
weather<-filter(weather,Year<1961)
baseParms$weather<-rbind(weather,baseParms$weather)


output<-do.call(SoNWaL,baseParms)
ff<-filter(output,Year>2014&Year<2019)
plot(ff$Rs*7.4)
plot(output$Nav)
plot(output$LAI)




sampleOutputOak<-function(df,sY,eY,swc=T){
  #convert to average grams per m2 per day depending on timestep of model
  modif<- if(nrow(df) < 3000) 1.6 else  7.142857
  nDays<- if(nrow(df) < 3000) 30 else  7
  
  df<- filter(df,Year>=sY&Year<=eY)
  
  df<-df%>%
    mutate(GPP=GPP*modif)%>%
    mutate(NEE=NEE*modif)%>%
    mutate(Reco=Reco*modif)
  
  dbh<-df %>% group_by(Year) %>%
    summarise(dbh = mean(dg)) %>%
    filter(Year %in% c(2009, 2011, 2019, 2020))
  
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       dbh$dbh
       )
  
  return(m)
}


## Likelihood function
LL_oak<- function(p){
  p<-p*.GlobalEnv$param_scaler
  baseParms[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
      output<-   do.call(SoNWaL,baseParms)
      modelled <- sampleOutputOak(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      mod_sim <- data.frame(modelled = modelled, obs = .GlobalEnv$obs, dev = .GlobalEnv$dev)
      mod_sim <- mod_sim[is.na(mod_sim$obs) == F,]
      NlogLik  <- ifelse(any(is.na(modelled) == T), -Inf, flogL(data=mod_sim$obs, sims=mod_sim$modelled, data_s=mod_sim$dev))

    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


priorVals<-createPriors_oak(oak=getParmsOak(E_S1=1,E_S2 = 1))[[1]]
param_scaler<-createPriors_oak(oak=getParmsOak(E_S1=1,E_S2 = 1))[[2]]


nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC","Q10")


likelihoodFunc<-LL_oak

#run in loop and write to file every 100k in case of errors or problems with JASMIN - allows for easy restarting of mcmc chain if something goes wrong
  iters=100000
  #Initiate bayesian setup
  settings = list(
    iterations = iters,
    startValue = 7, # internal chain number 
    nrChains = 1, # Number of chains
    burnin = round(iters/100*10), #10% burnin
    parallel = T,
    consoleUpdates = 10,
    message = TRUE)
  
  
  #on JASMIN I found you need to create bayesian setup even if re-starting a chain, I think as this initiates the cluster needed to run in parallel  
  BS3PGDN <- createBayesianSetup(likelihood = likelihoodFunc, prior = priorVals, names = nm, parallel = 8, catchDuplicates = F )

    out<- runMCMC(bayesianSetup =out, sampler = "DEzs", settings = settings)
  

    codM<-as.data.frame(out$chain[[1]])

    codM<-tail(as.data.frame(codM),1)
    names(codM)<-nm
    
    
    baseParms[nm]<-codM[nm]*param_scaler
    
    
    output<-do.call(SoNWaL,baseParms)
    ff<-filter(output,Year>=2013&Year<2019)
    plot(ff$Rs*7.4)
    plot(ff$Nav)
    plot(ff$GPP*7.14)
plot(output$dg)    
