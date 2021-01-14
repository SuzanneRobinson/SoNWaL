

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

###########################
## Initialise Parameters ##
###########################
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                wiltPoint = 0.1, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
                fieldCap =0.26,#Field capacity
                satPoint = 0.3, #field saturation point
                K_s=0.2, #Soil conductivity
                shared_area=3, #shared area of rooting and non-rooting zone
                V_nr=3, #Volume of non-rooting zone
                maxRootDepth=2.5,
                sigma_zR =0.6, #area/depth explored by 1kg of root biomass
                SWC_nr=10, #SWC of non-rooting zone at time 0
                E_S1 =0.1, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
                E_S2 =0.2, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
                MaxASW_state=50,
                K_drain=0.99,
                timeStp = 12 # time step, 52 for weekly, 12 for monthly and 365 for daily
                )
#######################################################

## Plot the timeseries of model output vs data
output<-do.call(fr3PGDN,sitka)
results<-plotResults(output,ShortTS=F)
results[1]

#get some previous run parameter estimates#
codM<-mergeChains(out$chain)
codM<-tail(as.data.frame(codM),1)
#codM<-data.table::transpose(data.frame(colMedians(codM)))
#names(codM)<-nm
sitka[nm]<-codM[nm]
#
### Run the 3PGN model using the sitka parameters
##.GlobalEnv$interRad<-0
#output<-do.call(fr3PGDN,sitka)
output<-do.call(fr3PGDN,sitka)
plot(output$SWC_nr[c(1:nrow(output))]~output$t[c(1:nrow(output))],col="white")
lines(output$SWC_rz[c(1:nrow(output))]~output$t[c(1:nrow(output))],col="red")
lines(output$SWC_nr[c(1:nrow(output))]~output$t[c(1:nrow(output))],col="blue")

  plot(.GlobalEnv$interRad[c(1:2347)]~output$t[c(1:2347)],col="white")
lines(.GlobalEnv$interRad[c(1:2347)]~output$t[c(1:2347)],col="red")
lines(.GlobalEnv$Rad[c(1:2347)]~output$t[c(1:2347)],col="blue")

#outVals<-as.data.frame(cbind(.GlobalEnv$Rad,.GlobalEnv$interRad))
#names(outVals)<-c("Rad","interRad")
#
#g1<-ggplot(outVals,aes(y=interRad,x=output$t[-1]))+
#  geom_line(col="red",alpha=0.7)+
#  ggtitle("Solar radiation reaching soil")+
#  ylab(expression(Solar~radiation~varphi))+
#  xlab("Time (years)")+
#  theme_bw()
#
#
#g2<-ggplot(outVals,aes(y=Rad,x=output$t[-1]))+
#  geom_line(col="red",alpha=0.7)+
#  ggtitle("Solar radiation")+
#  ylab(expression(Solar~radiation~varphi))+
#  xlab("Time (years)")+
#  theme_bw()
#
#ggarrange(g2,g1)
#
### Plot model outputs
#pOut <- plotModel(output)
