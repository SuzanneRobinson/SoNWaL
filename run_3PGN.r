

######################
## START THE SCRIPT ##
######################

## Load necessary packages
library(fr3PGDN,quietly=TRUE)
library(tidyverse,quietly=TRUE)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
## Years of data to use for calibration
startYear = 2015
endYear = 2018
install.packages("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\fr3PGDN_2.0.tar.gz", repos = NULL, type="source")


## Met data
clm_df_full<-getClimDat("monthly")
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
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
                K_s=0.1, #Soil conductivity
                shared_area=3, #shared area of rooting and non-rooting zone
                V_nr=4, #Volume of non-rooting zone
                maxRootDepth=2.5,
                sigma_zR =0.3, #area/depth explored by 1kg of root biomass
                SWC_nr=10, #SWC of non-rooting zone at time 0
                E_S1 =1, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                K_drain=0.16,
                timeStp = 12 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################

## Plot the timeseries of model output vs data
#i<-c(26,12,15,13)
#26 bad
#12 bad
#15 bad
#13 bad
#sitka[nm[-i]]<-codM[nm[-i]]*3

output<-do.call(fr3PGDN,sitka)
tail(output$GPP)
results<-plotResults(output,ShortTS=F)
#results
ggarrange(results[[1]],results[[2]],results[[3]],results[[4]],results[[5]],results[[12]])

codM<-getSample(out, start = 10000, coda = TRUE, thin = 10)
#bayesplot::mcmc_trace(codM)
priorSamp<-Uprior$sampler(35000)
MCMCtrace(codM[-7],wd="C:\\Users\\aaron.morris", priors = priorSamp, post_zm=F,iter=5000)
#get some previous run parameter estimates#
codM<-as.data.frame(mergeChains(out$chain))
names(codM)<-nm
codM<-colMedians(as.data.frame(codM))
#codM<-data.table::transpose(data.frame(colMedians(codM)))
#names(codM)<-nm
sitka[nm]<-codM[nm]
#
### Run the 3PGN model using the sitka parameters
##.GlobalEnv$interRad<-0
#output<-do.call(fr3PGDN,sitka)
output<-do.call(fr3PGDN,sitka)
plot(output$ASW[c(541:nrow(output))]~output$t[c(541:nrow(output))],col="white")
lines(output$ASW[c(541:nrow(output))]~output$t[c(541:nrow(output))],col="red")
lines(output$SWC_nr[c(1:nrow(output))]~output$t[c(1:nrow(output))],col="blue")


plot(output$ASW[c(541:nrow(output))]~output$t[c(541:nrow(output))],col="white")
lines(output$ASW[c(541:nrow(output))]~output$t[c(541:nrow(output))],col="red")
plot(output$ASW[c((541-12*3):(nrow(output)-12*3))]~output$t[c(541:nrow(output))],col="blue")

  plot(.GlobalEnv$interRad[c(1:2347)]~output$t[c(1:2347)],col="white")
lines(.GlobalEnv$interRad[c(1:2347)]~output$t[c(1:2347)],col="red")
lines(.GlobalEnv$Rad[c(1:2347)]~output$t[c(1:2347)],col="blue")

#outVals<-as.data.frame(cbind(.GlobalEnv$Rad,.GlobalEnv$interRad))
#names(outVals)<-c("Rad","interRad")
#
output<-filter(output, Year == "2017")

g1<-ggplot(output[-1,],aes(y=soilRad,x=t))+
  geom_line(col="red",alpha=0.7)+
  ggtitle("Solar radiation reaching soil")+
  ylab(expression(Solar~radiation~varphi~(J~m^2~day^-1)))+
  xlab("Time (years)")+
  theme_bw()
#
#
g2<-ggplot(output[-1,],aes(y=totalRad,x=t))+
  geom_line(col="red",alpha=0.7)+
  ggtitle("Total solar radiation")+
  ylab(expression(Solar~radiation~varphi~(J~m^2~day^-1)))+
  xlab("Time (years)")+
  theme_bw()

g3<-ggplot(output[-1,],aes(y=LAI,x=t))+
  geom_line(col="dark green",alpha=0.7)+
  ggtitle("LAI")+
  ylab(expression(Leaf~area~index))+
  xlab("Time (years)")+
  theme_bw()

g4<-ggplot(output[-1,],aes(y=potentialEvap,x=t))+
  geom_line(col="purple",alpha=0.7)+
  ggtitle("Potential evaporation")+
  ylab(expression(Potential~evaporation~(mm~day^-1)))+
  xlab("Time (years)")+
  theme_bw()

g5<-ggplot(output[-c(1:15),],aes(y=volSWC_rz,x=t))+
  geom_line(col="blue",alpha=0.7)+
  ggtitle("Volumetric SWC of root zone")+
  ylab(expression(theta[rz]~(m^-1~m^-3)))+
  xlab("Time (years)")+
  theme_bw()


g6<-ggplot(output[-c(1:15),],aes(y=ASW,x=t))+
  geom_line(col="blue",alpha=0.7)+
  ggtitle("Available soil water")+
  ylab(expression(ASW~(m^-1~m^-3)))+
  xlab("Time (years)")+
  theme_bw()


#
ggarrange(g2,g1,g3,g4)
#
### Plot model outputs
#pOut <- plotModel(output)
out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_outx_2021-02-265347.RDS")

#
#YC.finder <- function(HT,AGE) 
#  
#{
#  
#  YC.site = ((HT/(1-exp(-0.033329*AGE))^1.821054)-14.856317)/1.425397
#  
#  if(YC.site>24) 24
#  
#  else if (YC.site<4) 4
#  
#  else unlist(sub("\\]","",unlist(strsplit(as.character(cut(YC.site,breaks=c(6,8,10,12,14,16,18,20,22,24),right=T)),split="\\,"))[2]))
#  
#}
#
weather<-clm_df_full
#aggregate data by year for total annual ranfall
annualPrecip <- weather%>%group_by(Year)%>%summarise(annual_precip=sum(Rain))


hazprecip <- quantile(annualPrecip$annual_precip,0.1)

#Risk function
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

#' # Expected loss of GPP 
inpt   <- tibble(strtyr=c(1985,1995,2005,2015),endyr=c(1988,1998,2008,2018))
riskdf <- purrr::map2_df(inpt$strtyr,inpt$endyr,calc_risk,df=annualPrecip,hazval=hazprecip)
knitr::kable(riskdf, digits=4, align = c(rep("c", 5)))
