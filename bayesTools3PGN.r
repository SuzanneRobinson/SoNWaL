## Load necessary packages
library(fr3PGDN,quietly=TRUE)
library(BayesianTools,quietly=TRUE)
library(tidyverse)

## Switch off annoying warnings
options(warn=-1)

## Custom function to sample the modelled data frame
sampleOutput<-function(df,sY,eY){


    m<-c(filter(df,Year>=sY&Year<=eY)$GPP,
         filter(df,Year>=sY&Year<=eY)$NPP,
         filter(df,Year>=sY&Year<=eY)$NEE,
         filter(df,Year>=sY&Year<=eY)$Reco,
         filter(df,Year>=sY&Year<=eY)$Rs,
         filter(df,Year>=sY&Year<=eY)$Etransp,
         filter(df,Year==2015&Month==8)$LAI,
         filter(df,Year==2018&Month==8)$LAI,
         filter(df,Year==2018&Month==8)$N,
         filter(df,Year==2018&Month==8)$dg,
         filter(df,Year==2015&Month==7)$Wr,
         filter(df,Year==2015&Month==7)$difRoots,
         filter(df,Year==2015&Month==7)$totC,
         filter(df,Year==2015&Month==7)$totN
         )
    
    m
    
    return(m)
}




## Custom function to sample the modelled data frame for shorter time-steps
sampleOutputTS<-function(df,sY,eY){
  
 df<- filter(df,Year>=sY&Year<=eY)
  
 df$week<- rep(1:52,4)
  
 aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[,3]
 
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Rs~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=mean)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       filter(df,Year==2015&Month==7)$Wr[1],
       filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1]
  )
  
  m
  
  return(m)
}


## Custom function to get the percentage accepted in a chain
acPerc<-function(chain){
    cL=nrow(chain)
    nAccepted = length(which(chain[,ncol(chain)]==1))##length(unique(pChain[,1]))
    acceptance = (paste(nAccepted, "out of ", cL, "candidates accepted ( = ",round(100*nAccepted/(cL-cL*0.1)), "%)"))
    print(acceptance)
}

## Years of data to use for calibration
startYear = 2015
endYear = 2017

library(dplyr)
library(lubridate)
years <- 2010:2012


##################needs putting into a neater function################
if(Sys.info()[1]=="Windows"){
  clm.df.full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\clm_df_full.csv")
  clm.df.daily<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\weather_day.csv")
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
  }else
{
  clm.df.full<-read.csv("/home/users/aaronm7/3pgData/clm_df_full.csv")
  clm.df.daily<-read.csv("/home/users/aaronm7/3pgData/weather_day.csv")
  data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
  
  
}


#Add date
clm.df.full$date<-as.Date(paste(clm.df.full$Year,"-",clm.df.full$Month,"-01",sep=""))
clm.df.full$week<-week(clm.df.full$date)
clm.df.daily$Date<-as.Date(clm.df.daily$DOY, origin = paste0(clm.df.daily$Year,"-01-01"))
clm.df.daily$week<-week(clm.df.daily$Date)

clm.df.dailyX<-clm.df.daily[!duplicated(clm.df.daily[c(1,12)]),]

weeklyRfall<-aggregate(clm.df.daily$Rain~clm.df.daily$week+clm.df.daily$Year,FUN=sum)
names(weeklyRfall)<-c("week","Year","Rain")

weeklyRfall$Month<-month(clm.df.dailyX$Date)
weeklyRfall<-cbind(weeklyRfall,clm.df.dailyX[c(3:5,7)])
weeklyRfall$FrostDays<-(aggregate(clm.df.daily$FrostHours~clm.df.daily$week+clm.df.daily$Year,FUN=sum)[,3])/24
weeklyRfall$MonthIrrig<-0

#Split into weekly data
clm.df.fullX<-NULL
clm_df<-clm.df.full
for(i in c(1:nrow(clm.df.full))){
  reps<-ifelse(clm.df.full$Month[i]!=12,clm.df.full[i+1,"week"]-clm.df.full[i,"week"],4)
  clm_df[i,]$Rain<-clm_df[i,]$Rain/reps
  clm_df[i,]$MonthIrrig<-clm_df[i,]$MonthIrrig/reps
  
  clm.df.fullX<-rbind(clm.df.fullX,do.call("rbind",(replicate(reps, clm_df[i,], simplify = FALSE))))
}
clm.df.fullX$week<-c(1:51)

weeklyRfall<-weeklyRfall[,c(names(clm.df.fullX[,c(1:9)]))]
################################################################################


###############################################################

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
            MaxASW = 300,
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
            SWconst0 = 0.55,
            SWpower0 = 6,
            m0 = 0,
            MaxAge = 400,
            nAge = 4,
            rAge = 0.95,
            gammaFx = 0.01888,
            gammaF0 = 0.001,
            tgammaF = 36,
            Rttover = 0.017, #same as gammaR?
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
            ##Water balance submodel
            waterBalanceSubMods =F, #Whether to run model using updated water balance submodels
            theta_wp = 0.1, #Wilting point in m^3/m^3? need to convert to mm per meter with rooting depth?
            theta_fc =0.29,#Field capacity
            K_s=1, #Soil conductivity
            V_nr=10, #Volume of non-rooting zone
            sigma_zR =0.2, #area/depth explored by 1kg of root biomass
            sigma_nr0=250, #SWC of non-rooting zone at time 0
            E_S1 =100, #Cumulitive evap threshold
            E_S2 =0.1, #how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
            MaxASW_state=50,
            timeStp = 12
            
)


nm<-c("theta_wp","theta_fc","K_s","V_nr","sigma_zR","E_S1","E_S2","sigma_nr0","pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR","mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond","Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")

## Load the parameter set from the fifth calibration run to use as starting values
#load("fifth_run_3pgn_all_data/par.Rdata")

#sitka[nm]<-ff[which(names(ff)!="fNn"&names(ff)!="fN0")]

f.decrease <- c(0.1,0.2,0.001,1,0.2,1,0.001,100,0.588503613257886,0.752929538228874,0.956131627577964,0.050456035523466,0.384021499609213,0.250229439327847,0.57408236899746,0.909666760291794,0.853276910139941,0.974961101217424,1,0.636422367959785,0.732916669791679,0.443930919848964,0.741758519667562,0.816463641720414,0.221779786451702,0.303779963365252,1,0.00141038795075,0.730688961031379,0.899808741360758,0.024817372196732,0.99632339563598,0.996373181003088,0.999649942946159,0.996388219783102,0.998203040988276,0.998245174258832,0.97983098579238,0.913069476259938,0.961283723717706,0.950056672692535,0.893875965852296,0.991080780202615,0.990457295759556)

f.increase <- c(0.16,0.4,10,5,0.9,25,0.2,1000,0.573973679288588,0.235352308855631,1.86098081013281,0.374136113325978,0.231957000781575,0.56202200140032,3.45793787115991,1.30349761255926,0.600615525746093,0.251944939128821,0.768680943537667,0.817888160201076,0.335416651041606,0.668207240453109,0.549448881994627,0.835363582795864,0.03762695139773,0.218385064110809,0.917925202458998,2.5949226033773,1.15448831174897,0.001912586392424,5.82627839462287,6.35320872803933,2.62681899691161,2.50057053840858,2.61178021689853,0.796959011723578,0.754825741168422,1.01690142076198,0.738610474801243,0.935813814114711,0.498299819223935,1.12248068295408,0.783843959477034,0.90854084808886)


pMaxima <- as.vector(unlist(sitka[nm])*(1+(f.increase)))
pMinima <- as.vector(unlist(sitka[nm])*(1-(f.decrease)))
pValues <- as.vector(unlist(sitka[nm]))

pMaxima[1:8] <- f.increase[1:8]
pMinima[1:8] <- f.decrease[1:8]

output<-do.call(fr3PGDN,sitka)

modelled <-sampleOutput(output,startYear,endYear)

observed <- c(data$gpp,                ## GPP
              data$npp,                ## NPP
              data$nee,                ## NEE
              data$reco,               ## Reco
              data$rs,                 ## Rs
              data$et,                 ## Etransp
              data$gs[2:nrow(data)],   ## CanCond
              5.7,5.56,                ## LAI
              1348,                    ## N
              24.1,                    ## dg
              4.88,                    ## Wr
              0.53,                    ## difRoots
              86.7,                    ## totC, see jarvis_total_soil.ods
              2.16                     ## totN, 40 C:N ratio
              )

dev <- c(rep(.001,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.001,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.01,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.01,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(0.01,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(0.01,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(0.5,(nrow(dplyr::filter(data,year>=startYear&year<=endYear))-1)),
         1.5,1.5,
         1,
         3,
         2,
         1,
         5,
         1
         )



## Likelihood function

NLL <- function(p){
    sitka[.GlobalEnv$nm]<-p
    output<-do.call(fr3PGDN,sitka)
    #use sampleOutputTS if using smaller time-steps
    modelled <-sampleOutput(output,.GlobalEnv$startYear,.GlobalEnv$endYear)

    NlogLik  <- sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=dev,log=T))
    return(NlogLik)
}


## pChain <- as.matrix(out$chain[,1:(ncol(as.matrix(out$chain))-3)])
## pMean <- apply(pChain,2,mean)
## pSD <- apply(pChain,2,sd)

## ~~~~~ ## 
## Prior ##
## ~~~~~ ##

## This is a uniform prior distribution
pMaxima[[27]]<-0.01
Uprior <- createUniformPrior(lower = pMinima, upper = pMaxima)

## This is a truncated normal distribution prior based on the chain produced
## from the fifth calibration run
## Nprior <- createTruncatedNormalPrior(mean = pMean, sd = pSD, lower = pMinima, upper = pMaxima)


## load("fifth_run_3pgn_all_data/out.Rdata")
## newPrior <- createPrior(out,lower=pMinima,upper=pMaxima)

## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Create the Bayesian setup ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
BS3PGDN <- createBayesianSetup(likelihood = NLL, prior = Uprior, names = nm, parallel = T, catchDuplicates = F )



## Use the previous calibration parameters as starting values
## mP<-matrix(rep(pValues,3),nrow=length(pValues),ncol=3)



## Choose the settings for the run

settings = list(
  iterations = 50000,
  ## Z = NULL,
  startValue = 30, #t(mP),#NULL, ## Use 5 chains instead of 3
  nrChains = 1,
  pSnooker = 0.5,
  burnin = 10,
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

##  Run the Monte Carlo Markov Chain

out <- runMCMC(bayesianSetup = BS3PGDN, sampler = "DEzs", settings = settings)

saveRDS(out,file="outx.RDS")


## rM <- function(p){
##     newSitka <- .GlobalEnv$sitka
##     newSitka[.GlobalEnv$nm] <- p
##     output <- do.call(fr3PGDN,newSitka)
##     modelled <-sampleOutput(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
##     return(modelled)
## }

## Get the median from the full chain
df <- as.data.frame(as.matrix(out$chain))
par <- as.vector(sapply(df,FUN=median))[1:length(nm)]

## Write the original parameter set and the newly calibrated parameter set
options(scipen=10)
parameters <- data.frame(names=nm,original=as.vector(unlist(sitka[nm])),calibrated=par)
save(parameters,file="par.Rdata")
