## Load necessary packages
library(fr3PGDN)
library(BayesianTools)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)

##choose timestep size
timeStep<-"monthly"

#create filename based on time...maybe not necessary, just trying to avoid overwriting of outputs from diff sessions by JASMIN
fName=paste0("outx_",stringr::str_sub(Sys.time(), 0, -10),stringr::str_sub(Sys.time(), 15, -4),stringr::str_sub(Sys.time(), 18, -1),".RDS")
#get climate data
clm_df_full<-getClimDat(timeStep)
## Read Harwood data for Sitka spruce and mutate timestamp to POSIXct
if(Sys.info()[1]=="Windows"){
  data <- read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
}else
{
  data <- read.csv("/home/users/aaronm7/3pgData/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
}

#get parameter values 

sitka<-getParms(timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)



nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")



##Set priors
priors<-createPriors_sitka(sitka=sitka)

pMaxima<-priors[[1]]
pMinima<-priors[[2]]
pMaxima[[30]]<-0.1
Uprior <- createPrior(lower = pMinima, upper = pMaxima)
## Set observed calibration data and uncertainty
startYear = 2015
endYear = 2018


Uprior <- createTruncatedNormalPrior(mean = priors[[3]], sd=(pMaxima-pMinima)*0.1,
                     lower = pMinima, upper = pMaxima)

priorSDfunc<-function(Uprior){
  prVals<-Uprior$sampler(10000)
  prVals[,7]<-exp(prVals[,7])
  prVals[,8]<-exp(prVals[,8])
  prValsX<-as.data.frame(round(HPDinterval(as.mcmc(prVals)),3))
  prValsX$med<-colMedians(prVals)
  prValsX$names<-nm
}


observed <- c(data$gpp,                ## GPP
              data$npp,                ## NPP
              data$nee,                ## NEE
              data$reco,               ## Reco
              data$rs,                 ## Rs
              data$et,                 ## Etransp
              #  data$gs[2:nrow(data)],   ## CanCond
              5.7,5.56,                ## LAI
              1348,                    ## N - fairly well known
              24.1,                    ## dg
              #  4.88,                    ## Wr
              # 0.53,                    ## difRoots
              429.52,                    ## totC, see jarvis_total_soil.ods
              14.30,                     ## totN, 40 C:N ratio
              data$swc                ## SWC
              
              #measured values are from 2015 for totC and totN
              #possible starting value from other site using approx values
)


dev <- c(rep(.3,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.3,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.3,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.3,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(.1,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         rep(4,nrow(dplyr::filter(data,year>=startYear&year<=endYear))),
         # rep(0.5,(nrow(dplyr::filter(data,year>=startYear&year<=endYear))-1)),
         0.1,0.1,
         10,
         1,
         #  2,
         #  1,
         5,
         0.5,
         rep(.1,nrow(dplyr::filter(data,year>=startYear&year<=endYear)))
         
)


observed<-if(sitka$waterBalanceSubMods==FALSE)observed[-c(295:342)] else observed
dev<-if(sitka$waterBalanceSubMods==FALSE)dev[-c(295:342)] else dev

likelihoodFunc<-ifelse(sitka$waterBalanceSubMods==T,NLL,NLL_noHYD)

outX<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_pine_outx_2021-03-295236.RDS")

for (i in c(1:5)){
iters=100000
#Initiate bayesian setup
BS3PGDN <- createBayesianSetup(likelihood = likelihoodFunc, prior = Uprior, names = nm, parallel = 8, catchDuplicates = F )
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
out3 <- runMCMC(bayesianSetup = BS3PGDN, sampler = "DEzs", settings = settings)


#Save output
saveRDS(out2,file=paste0(timeStep,"_",fName))

}