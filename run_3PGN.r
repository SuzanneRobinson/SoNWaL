

######################
## START THE SCRIPT ##
######################

## Load necessary packages
library(fr3PGDN,quietly=TRUE)
library(tidyverse,quietly=TRUE)
library(lubridate)
## Years of data to use for calibration
startYear = 2015
endYear = 2018


## Met data
clm_df_full<-getClimDat("monthly")
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
#data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))

###########################
## Initialise Parameters ##
###########################
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                wiltPoint = 0.12, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
                fieldCap =0.26,#Field capacity
                satPoint = 0.4, #field saturation point
                K_s=0.2, #Soil conductivity
                shared_area=3, #shared area of rooting and non-rooting zone
                V_nr=4, #Volume of non-rooting zone
                maxRootDepth=1.5,
                sigma_zR =0.1, #area/depth explored by 1kg of root biomass
                SWC_nr=10, #SWC of non-rooting zone at time 0
                E_S1 =5, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                K_drain=0.1,
                timeStp = 12 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################

## Plot the timeseries of model output vs data
output<-do.call(fr3PGDN,sitka)
results<-plotResults(output,ShortTS=F)
results[1]
ggarrange(results[[1]],results[[2]],results[[3]],results[[4]],results[[5]],results[[6]])

codM<-getSample(out, start = 100000, coda = TRUE, thin = 10)
bayesplot::mcmc_trace(codM)
#get some previous run parameter estimates#
codM<-mergeChains(out$chain)
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
g1<-ggplot(output[-1,],aes(y=soilRad,x=t))+
  geom_line(col="red",alpha=0.7)+
  ggtitle("Solar radiation reaching soil")+
  ylab(expression(Solar~radiation~varphi~(J~m^2~day^-1)))+
  xlab("Time (years)")+
  theme_bw()
#
#
g2<-ggplot(output[-1,],aes(y=netRad,x=t))+
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

#
ggarrange(g2,g1,g3,g4)
#
### Plot model outputs
#pOut <- plotModel(output)
