## Load necessary packages
library(fr3PGDN,quietly=TRUE)
library(BayesianTools,quietly=TRUE)
library(tidyverse)

#get climate data
clm_df_full<-getClimDat("monthly")
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))

#get parameter values 
sitka<-getParms()

##Set priors
priors<-createPriors_sitka(sitka=sitka)
pMaxima<-priors[[1]]
pMinima<-priors[[2]]
pMaxima[[30]]<-0.01
Uprior <- createUniformPrior(lower = pMinima, upper = pMaxima)

## Set observed calibration data and uncertainty
startYear = 2015
endYear = 2017
years <- 2010:2012
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

##from george's code
dev <- c(rep(3,nrow(filter(data,year>=startYear&year<=endYear))),
         rep(5,nrow(filter(data,year>=startYear&year<=endYear))),
         rep(5,nrow(filter(data,year>=startYear&year<=endYear))),
         rep(5,nrow(filter(data,year>=startYear&year<=endYear))),
         rep(5,nrow(filter(data,year>=startYear&year<=endYear))),
         rep(10,nrow(filter(data,year>=startYear&year<=endYear))),
         rep(0.05,(nrow(filter(data,year>=startYear&year<=endYear))-1)),
         1.5,1.5,
         100,
         3.0,
         2.0,
         1.0,
         30.0,
         1.0
)

#Initiate bayesian setup
BS3PGDN <- createBayesianSetup(likelihood = NLL, prior = Uprior, names = nm, parallel = T, catchDuplicates = F )
settings = list(
  iterations = 10000,
  ## Z = NULL,
  startValue = 8, #t(mP),#NULL, ## Use 5 chains instead of 3
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


#Run calibration
out <- runMCMC(bayesianSetup = BS3PGDN, sampler = "DEzs", settings = settings)

#Save output
#saveRDS(out,file="outx.RDS")


