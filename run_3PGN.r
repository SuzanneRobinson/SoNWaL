

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
clm.df.full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\clm_df_full.csv")

## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))


#######################################################
## Parameters Xenakis 2020 (unpublished) ##
#######################################################
#not sure if monthly rates need to be modified to whatever timestep is being used, depends on how they are used in the model
#may be easier to adjust them by timestep within the model rather than at proposal to keep things cleaner?
getParms<-function(){
sitka<-list(weather=clm.df.full,
            ## ~~ Initial pools ~~ ##
            Wl = 0.01,
            WlDormant = 0,
            Wr = 0.01,
            Wsbr = 0.1,
            Wlitt = 0,
            YrC = 42.95,
            YlC = 85.90,
            OC = 300.66,
            YrN = 1.43,
            YlN = 2.86,
            ON = 10.01,
            Nav = 8.76,
            ## ~~ Site ~~ ##
            N = 2000,
            rotation = 1,
            cycle = 1,
            rm.sprouts = F,
            nyears = 35,
            initial.month = 1,
            latitude = 57.06,
            soilclass = -1,
            ASW = 165,
            MaxASW = 500,
            MinASW = 0,
            CO2 = 400,
            ## ~~ Parameters ~~ ##
            pFS2 = 1.4,
            pFS20 = 0.8,
            aS = 0.138,
            nS = 2.3,
            pRx = 0.45,
            pRn = 0.3,
            Tmin = -5,
            Topt = 15,
            Tmax = 35,
            kF = 1,
            SWconst0 = 0.7,
            SWpower0 = 9,
            m0 = 0,
            MaxAge = 400,
            nAge = 4,
            rAge = 0.95,
            gammaFx = 0.01888,
            gammaF0 = 0.001,
            tgammaF = 36,
            Rttover = 0.017,
            MaxCond = 0.02,
            LAIgcx = 3.33,
            BLcond = 0.2,
            wSx1000 = 500,
            thinPower = 1.5,
            mF = 0.5,
            mR = 0.3,
            mS = 0.2,
            SLA0 = 5,
            SLA1 = 3,
            tSLA = 3,
            k = 0.5,
            fullCanAge = 15,
            MaxIntcptn = 0.15,
            LAImaxIntcptn = 5,
            alpha = 0.06,
            Y = 0.47,
            poolFractn = 0,
            e20 = 2.2,
            rhoAir = 1.2,
            lambda = 2460000,
            VPDconv = 0.000622,
            fracBB0 = 0.15,
            fracBB1 = 0.15,
            tBB = 10,
            rhoMin = 0.55,
            rhoMax = 0.55,
            tRho = 5,
            Qa = -90,
            Qb = 0.8,
            gDM_mol = 24,
            molPAR_MJ = 2.3,
            CoeffCond = 0.05,
            fCalpha700 = 1.433,
            fCg700 = 0.451,
            fCalphax = 2.33333333333333,
            fCg0 = 1.75,
            MinCond = 0.015,
            klmax = 0.01,
            krmax = 0.00423943,
            komax = 0.00045886,
            hc = 0.2,
            qir = 334.290515,
            qil = 49.0841127,
            qh = 23.6348669,
            qbc = 2.21427684,
            el = 0.24636719,
            er = 0.56122150,
            Nf = 0.00684,
            Navm = 0.01,
            Navx = 10,
            leaf.grow = 0,
            leaf.fall = 0,
            Wl.s = 0.01,
            Wsbr.s = 0.1,
            Wr.s = 0.01,
            pWl.sprouts = 0.5,
            pWsbr.sprouts = 0.9,
            cod.pred = "3PG",
            cod.clim = "Month",
            ## ~~ Almedia et al. Parameters ~~ ##
            waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
            theta_wp = 0.1, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
            theta_fc =0.29,#Field capacity
            theta_sat= 1, #field saturation point
            K_s=0.1, #Soil conductivity
            shared_area=4, #shared area of rooting and non-rooting zone
            V_nr=3, #Volume of non-rooting zone
            maxRootDepth=2,
            sigma_zR =0.3, #area/depth explored by 1kg of root biomass
            SWC_nr=10, #SWC of non-rooting zone at time 0
            E_S1 =0.3, #Cumulitive evap threshold (kg^m-2) - sensitive to length of time-step, e.g. monthly time-step means wetting event only occurs at end of month
            E_S2 =0.3, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
            MaxASW_state=50,
            K_drain=0.1,
            timeStp = 12 # time step, 52 for weekly, 12 for monthly and 365 for daily
            )
}
#######################################################


#Key of terminology 
#NPP - net primary production
#GPP - gross primary production
#NEE - Net ecosystem exchange
#Rs - soil respiration
#Reco - ecosystem respiration

## Plot the timeseries of model output vs data
output<-do.call(fr3PGDN,sitka)
results<-plotResults(output,ShortTS=F)
results[1]


## Calculate yield class from height
#output <- output%>%mutate(
#  yct = ((hdom/(1-exp(-0.033329*t.proj))^1.821054)-14.856317)/1.425397,
#  YC = ifelse(yct>24,24,ifelse(yct<4,4,yct))
#)

#get some previous run parameter estimates#
#codM<-getSample(out, start = 1000, coda = TRUE, thin = 1)[[1]]
##codM<-tail(as.data.frame(codM),5)
#codM<-data.table::transpose(data.frame(colMedians(codM)))
#names(codM)<-nm
#sitka[nm]<-codM
#
### Run the 3PGN model using the sitka parameters
##.GlobalEnv$interRad<-0
#output<-do.call(fr3PGDN,sitka)
output<-do.call(fr3PGDN,sitka)
plot(output$SWC_nr[c(1:2393)]~output$t[c(1:2393)],col="white")
lines(output$SWC_rz[c(1:2393)]~output$t[c(1:2393)],col="red")
lines(output$SWC_nr[c(1:2393)]~output$t[c(1:2393)],col="blue")

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
