
#############################################
#  _____             _ _   _     _ _        #
# |   __|___ ___ ___|_| |_|_|_ _|_| |_ _ _  #
# |__   | -_|   |_ -| |  _| | | | |  _| | | #
# |_____|___|_|_|___|_|_| |_|\_/|_|_| |_  | #
#                                     |___| #
############################################# 
#########Morris sensitivity analysis#########

library(fr3PGDN)
library(BayesianTools)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)
library(sensitivity)
library(viridis)

##choose timestep size ("monthly", "weekly" or "daily")
timeStep<-"weekly"

#Directory where climate data is stored (default is "data" in SonWal folder)
climDir<-("Data/")

## read in and format climate data
clm_df_full<-data.frame(getClimDatX(timeStep,climDir))%>%
  filter(Year<2019)


## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  flxdata_daily <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
}else
{
  data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
}

#get parameter values 
sitka<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)

#param names to test for sensitivity
nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")


#read in current MCMC chains
out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_2_T.RDS")
#get some previous run parameter estimates from chains
codM<-as.data.frame(mergeChains(out$chain))
names(codM)<-nm
codM<-colMedians(as.data.frame(codM))
sitka[nm]<-codM[nm]

##Set priors
priorVals<-createPriors_sitka(sitka=sitka)
pMaxima<-priorVals$upper
pMinima<-priorVals$lower
#pMaxima[[30]]<-0.1
#Uprior <- createPrior(lower = pMinima, upper = pMaxima)

#Select years for fitting to
startYear = 2015
endYear = 2018

observed<-observedVals(timeStep = "monthly",data=flxdata_daily,sY=startYear,eY=endYear)[[1]]
dev<-observedVals(timeStep =  "monthly",data=flxdata_daily,sY=startYear,eY=endYear)[[2]]

#Choose likelihood function depending on length of time-step (NLL_weekly works for both weekly or daily)
likelihoodFunc<-NLL


##setup likelihood and model running functions
morris_setup <- createBayesianSetup(
  likelihood = likelihoodFunc, 
  prior = priorVals, 
  names = nm)

set.seed(50)

#run morris sensitivity analysis
morrisOut <- morris(
  model = morris_setup$posterior$density,
  factors = nm, 
  r = 25, 
  design = list(type = "oat", levels = 25, grid.jump = 3), 
  binf = pMinima, 
  bsup = pMaxima, 
  scale = TRUE)

# summarise the sensitivity analysis
morrisOut_df <- data.frame(
  parameter = nm,
  mu.star = apply(abs(morrisOut$ee), 2, mean, na.rm = T),
  sigma = apply(morrisOut$ee, 2, sd, na.rm = T)
) %>%
  arrange( mu.star )


#plot results
wSitka<-morrisOut_df %>%
  gather(variable, value, -parameter) %>%
  #top_n(n=20, value) %>%
  ggplot(aes(reorder(parameter, value), value, fill = variable), color = NA)+
  geom_bar(position = position_dodge(), stat = 'identity') +
  scale_fill_viridis_d("", option = "D", labels = c('mu.star' = expression(mu * "* (high influence) "), 
                                      'sigma' = expression(sigma ~"(high interaction) "))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
    axis.title = element_blank(),
    legend.position = c(0.05 ,0.95),legend.justification = c(0.05,0.95)
  )

morrisOut_df$mu.star<-log10(1+morrisOut_df$mu.star)
morrisOut_df$sigma<-log10(1+morrisOut_df$sigma)
