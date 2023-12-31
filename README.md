Readme
================

# Installation

Open the file or load the package from within R studio "SoNWaL.proj"

run the following command to install the package:

``` r
devtools::install()
```

# Run the model and plot some quick outputs

## Setup

Firstly setup the data for running the model:

-   load in required packages
-   set start and end year
-   set required time-step
-   load in the climate and observed flux datasets

``` r
## Load necessary packages
library(SoNWaL)
library(tidyverse)
library(lubridate)
library(coda)
library(BayesianTools)
library(miscTools)
library(ggpubr)
library(matrixStats)
library(future)
library(furrr)
library(parallel)

# Years of data to use for calibration
startYear = 2015
endYear = 2018

# Time step to run SonWal with
timeStep<-"weekly"

# Directory where climate data is stored (default is "data" in SonWal folder)

climDir<-("data\\")

# read in and format climate data
clm_df_full<-data.frame(getClimDatX("weekly",climDir))%>%
  filter(Year<2019)
# Read Harwood data for Sitka spruce and mutate timestamp to POSIXct

  flxdata_daily <- read.csv("data\\harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
```

## Load in some example parameters

-   Initiating the parameter list with some default params
-   defining which parameters you wish to update
-   loading in a sample from an MCMC calibration ("exampParams.RDS") and in this instance taking the median values
-   update the parameter list with these values

``` r
# load in default parameters
sitka<-getParms(weather=clm_df_full,
                waterBalanceSubMods =T, #Whether to run model using updated water balance submodels
                timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365 #time step, 52 for weekly, 12 for monthly and 365 for daily
                )

# Names of fitted parameters  
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
        "SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")

# update with some example calibrated parameters (parameters are sample from full MCMC chain of Harwood fitting)
exampParams<-as.data.frame(readRDS("data//exampParams.RDS"))
exParms<-miscTools::colMedians(as.data.frame(exampParams))
names(exParms)<-nm
sitka[nm]<-exParms[nm]
```

# Run SoNWaL

To run SoNWaL you can then use the following command:

``` r
# run SoNWal
output<-do.call(SoNWaL,sitka)

# view some of the outputs, output contains time series for over 100 variables, only the first 10 shown here
head(output[, c(1:10)])
```

    ##       t.proj Year Month          t     hdom    N         G       dg        Vu
    ## 1 0.00000000   NA    NA 0.00000000 2.410997 2000 0.3068629 1.397694 0.3589744
    ## 2 0.01923077 1973     1 0.01923077 2.411089 2000 0.3068764 1.397725 0.3591332
    ## 3 0.03846154 1973     1 0.03846154 2.411069 2000 0.3068764 1.397725 0.3592696
    ## 4 0.05769231 1973     1 0.05769231 2.411048 2000 0.3068764 1.397725 0.3594059
    ## 5 0.07692308 1973     1 0.07692308 2.411028 2000 0.3068764 1.397725 0.3595420
    ## 6 0.09615385 1973     2 0.09615385 2.411007 2000 0.3068764 1.397725 0.3596779
    ##          LAI
    ## 1 0.01091896
    ## 2 0.01092526
    ## 3 0.01092217
    ## 4 0.01091892
    ## 5 0.01091549
    ## 6 0.01091191

## Quick plot

To get a quick plot, you can then use the quickPlot function with this output, however this does not include any credible intervals

``` r
# plotting QUICK PLOT - DOES NOT INCLUDE UNCERTAINTY! BUT QUICK :) -grouping aggregates the data, can be either "week" or "month"
quickPlot(flxdata_daily,output,grouping="month")
```

![](README_files/figure-markdown_github/quick_plot-1.png)

## Full plot

To include uncertainty, we need to sample from the posterior of the MCMC calibration, we can do this by sampling from the exampParams, which are in themselves a sample of the full MCMC chain.

``` r
# full plots - much slower but gives credible intervals and uncertainty of observed data - num samps is how many samples from posterior to use (>=500 ideal but 50-100 will give a pretty solid output for a quick checking)
results<-plotResultsNewMonthly(output,ShortTS=T,out=exampParams,numSamps = 25)
ggarrange(results[[1]],results[[2]],results[[8]],results[[3]],results[[5]],results[[4]],ncol=2,nrow=3)
ggarrange(results[[15]],results[[9]],results[[10]],results[[11]],results[[13]],results[[14]],ncol=2,nrow=3)
```

# Calibrations

To calibrate SoNWaL to new data

``` r
library(BayesianTools)
library(tidyverse)
library(dplyr)
library(coda)
library(miscTools)
library(SoNWaL)
```

If running on JASMIN or using a similar SLURM queue, the arguments can be stored in the batch file and passed to the script using the commandsArgs command if running locally these values can be passed manually

``` r
#read in arguments from batch file
args=(commandArgs(TRUE)) # if running from SLURM queue
args<-c("weekly","weekly_1_","1","T") # or if running locally
print(args) # print args so they are recorded in log file
timeStep=args[1] # time step to use e.g. weekly/monthly
chainID=args[2] # id number for the chain if doing multiple duplicate runs
chainNum=args[3] # numeric value for the chain number
```

an example of a SLURM batch file for running three repeats of the calibration script on JASMIN is as follows, with the last line being the arguments passed to the run script

    #!/bin/bash
    #SBATCH --partition=par-single
    #SBATCH --job-name=run1_cal3pgn
    #SBATCH -o w_%j.out
    #SBATCH -e w_%j.err
    #SBATCH --time=48:00:00
    #SBATCH -n 8
    #SBATCH --array=1-3
    module add jasr
    Rscript 3pgn_weekly_cal_SLURM_sitka.r weekly weekly_${SLURM_ARRAY_TASK_ID}_ $SLURM_ARRAY_TASK_ID T F

read in the climate data that will be used for calibration, climate data should always contain the Year, week, Month, Tmax, Tmin, Tmean, Rain, SolarRad, FrostDays, MonthIrrig (although often 0), VPD and RH, other data can also be stored here for passing to for example the likelihood function, here SWC is read in aswell.

**N.B. Climate data should be in the format of the time-step you are running SoNWaL in, e.g. weekly or monthly. SonWaL will run in whatever time-step the climate data is in and will adjust parameters used in the model from monthly rates on the fly**

``` r
#get climate data
#Directory where climate data is stored (default is "data" in SonWal folder)
climDir<-("data/")

## read in and format climate data
clm_df_full<-data.frame(getClimDatX(timeStep = "weekly",climDir))%>%
  filter(Year<2019)

head(clm_df_full)
```

    ##   Year week Month      Tmax       Tmin      Tmean Rain SolarRad FrostDays
    ## 1 1973    1     1 11.419694 -0.8494944  3.9919430    0 1.793574         2
    ## 2 1973    2     1 10.669222 -1.7075278  3.3847896    0 1.947364         4
    ## 3 1973    3     1  5.715594 -2.5673472  0.3225012    0 2.115103         6
    ## 4 1973    4     1  8.011911 -2.2194306  2.4256957    0 2.011973         4
    ## 5 1973    5     2  6.862892 -3.7515417 -0.2422386    0 2.221743         7
    ## 6 1973    6     2  6.880636 -1.9465639  2.1477353    0 6.183288         2
    ##   MonthIrrig        VPD       RH       SWC rainDays
    ## 1          0 0.05487941 93.02188 0.2798393        0
    ## 2          0 0.05629152 93.14685 0.2858093        0
    ## 3          0 0.05908644 91.70545        NA        0
    ## 4          0 0.02794263 96.16943        NA        0
    ## 5          0 0.07603817 86.90638        NA        0
    ## 6          0 0.05825735 91.90626 0.2801092        0

read in observed data, here we use flux data from Harwood and convert the data to single vectors of observed values and their standard deviations

``` r
flxdata <- read.csv("/data/harwood_data.csv")%>%mutate(timestamp=as.POSIXct(timestamp))
flxdata_daily <- read.csv("/data/harwood_daily.csv")%>%mutate(timestamp=as.POSIXct(timestamp))

observed<-observedVals(timeStep = "monthly",data=flxdata_daily,
                       sY=startYear,eY=endYear,swc=ifelse(args[4]=="T",TRUE,FALSE))[[1]]
dev<-observedVals(timeStep =  "monthly",data=flxdata_daily,
                  sY=startYear,eY=endYear, swc=ifelse(args[4]=="T",TRUE,FALSE))[[2]]
```

## likeihood function setup

-   A simple likelihood function:

-   here we firstly take the model output, and produce a vector of simulated data which exactly matches the vector of observed data extracted above in terms of order of outputs and length of time.

-   e.g. if your observed data is a vector of 3 years of GPP, NPP and NEE in weekly time-steps, this would be a vector 468 values in length (3 X 52 X 3), and your simulated data should also match this and be in the same order.

-   We then obtain a likelihood value for each simulated value against it's corresponding observed value, here we use a probability distribution proposed by Sivia, which reduces the weighting of outliers often found in flux data.

-   To run the likelihood function within bayesian tools, we have a single input *p*, this contains the list of proposed parameters which are updated by the MCMC algorithm at each iteration and for which we want to obtain a likelihood value for. **The output of SoNWaL is in tons per hectare per time step for flux data, so will need to be converted to whatever your observed data is, here the harwood data is in grams per m2 per day. in the sampleOutput function, tons per hectare per time step is converted using the modif coefficient, evapotranspiration is also converted to daily values. Depending on your observed data you will need to make sure the simulated data is converted to matching units!**

``` likelihood

## simulated data extraction function
#' @param sim output of SoNWaL model
#' @param sY start year of observed data
#' @param eY end of year observed data
#' @param swc whether to included soil water data in calibration
#' @return return vector of simulated data
sampleOutput<-function(sim,sY,eY,swc=T){
  #convert to average grams per m2 per day depending on timestep of model (using output length to quickly get time-step)
  modif<- if(nrow(sim)<1000) 1.6 else  7.142857
  nDays<- if(nrow(sim)<1000) 30 else  7
  
  sim<- filter(sim,Year>=sY&Year<=eY)

  sim<-sim%>%
    mutate(GPP=GPP*modif)%>%
    mutate(NPP=NPP*modif)%>%
    mutate(NEE=NEE*modif)%>%
    mutate(Rs=Rs*modif)%>%
    mutate(Reco=Reco*modif)
  
  
  m<-c(aggregate(sim$Rs~ sim$Month+sim$Year,FUN=mean)[,3],
       aggregate(sim$GPP~ sim$Month+sim$Year,FUN=mean)[,3],
       aggregate(sim$NEE~ sim$Month+sim$Year,FUN=mean)[,3],
       aggregate(sim$EvapTransp/nDays~ sim$Month+sim$Year,FUN=mean)[,3],
       filter(sim,Year==2015&Month==8)$LAI[1],
       filter(sim,Year==2018&Month==8)$LAI[1],
       filter(sim,Year==2018&Month==8)$N[1],
       filter(sim,Year==2018&Month==8)$dg[1],
       filter(sim,Year==2015&Month==7)$totC[1],
       filter(sim,Year==2015&Month==7)$totN[1],
      if(swc==T) {aggregate(sim$volSWC_rz~ sim$Month+sim$Year,FUN=mean)[,3]}
  )
  m
  return(m)
}


# sivia likelihood calculation
#' @param sims simulated data
#' @param data observed data values
#' @param data_s obsereved data standard deviation
flogL <- function(sims,data,data_s)
{ 
  Ri         <- (sims - data) / data_s
  i0         <- which( abs(Ri)<1.e-08 )
  
  logLi      <- log(1-exp(-0.5*Ri^2)) - log(Ri^2) - 0.5*log(2*pi) - log(data_s)
  logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])
  
  sum(logLi)
}



## Likelihood function
#' @param p proposed parameter values from MCMC algorithm
LL_sitka<- function(p){
  
  # get proposed parameters from MCMC algorthm and update parameter list for running SoNWaL
  p<-p*.GlobalEnv$param_scaler
  sitka[.GlobalEnv$nm]<-p
  
  #trycatch used incase model crashes or produces NA's return -Inf
   NlogLik <- tryCatch(
    {
      # run SoNWaL with proposed parameters
      output<-   do.call(SoNWaL,sitka)  
      #get vector of simualted values
      modelled <-sampleOutput(output,.GlobalEnv$startYear,.GlobalEnv$endYear,swc=sitka$waterBalanceSubMods)
      #obtain likihood of simualted values produced by proposed parameters
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,flogL(data=.GlobalEnv$observed,sims=modelled,data_s=.GlobalEnv$dev))
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}
```

## parameters

-   SoNWaL contains functions for creating parameter lists for sitka spruce and scots pine, for other species functions should be written which output matching lists, with corresponding parameter names and a list component which contains all the climate data for running the model

-   timeStp is an important parameter which tells SoNWaL the number of time-steps in a single year, e.g. for a weekly time-step it would be 52, this allows adjustments of monthly parameter rates etc. within the model

-   pseudo is for running the model using pseudo daily time-steps within the hydrological sub-model functions - see Almedia and Sands 2016 for details on the use of pseudo daily time-steps, in general i've found them not as good as using weekly time-steps but better than monthly - some strange behaviour has been noticed though so use with caution.

``` r
paramList<-getParms(
  waterBalanceSubMods=ifelse(args[4]=="T",TRUE,FALSE), timeStp = if (timeStep == "monthly") 12 else if (timeStep == "weekly") 52 else 365)

paramList$pseudo<-ifelse(args[5]=="F",FALSE,TRUE)
print(paste0("pseudo = ",args[5]))
```

## priors

-   Functions are available within SoNWaL which create priors for sitka spruce and scots pine based on previously published values, similar functions can be made for any new species

-   parameter scaler is a vector of values for scaling parameters to similar limits (this helps avoid problems potentially being caused by parameter values being on wildly different scales)

-   nm is a vector of the names of the parameters to calibrate

-   start and end year are the years the observed data goes to and from

``` r
startYear = 2015
endYear = 2018

##Set priors
priorVals<-createPriors_sitka(sitka=getParms(E_S1=1,E_S2 = 1))[[1]]
param_scaler<<-createPriors_sitka(sitka=getParms(E_S1=1,E_S2 = 1))[[2]]

nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
      "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
      "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
      "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k", "startN","startC","kF")

startYear = 2015
endYear = 2018
```

## run the calibration

-   As a calibration often takes around 48 hours or longer, I run the calibrations in a loop with 50000 or so MCMC iterations per loop iteration, and save as I go, this avoids problems of losing all progress if something goes wrong, or JASMIN times-out.

``` r
for (i in c(1:15)){
iters=50000
#Initiate bayesian setup
settings = list(
  iterations = iters,
  startValue = 7, # internal chain number (see diff evolution algorithm documentation)
  nrChains = 1, # Number of chains
  pSnooker = 0.5,
  burnin = round(iters/100*10), #10% burnin
  parallel = T,
  message = TRUE)

#check if files already exist and restart the chain or start from initial
if(file.exists(paste0("/home/users/aaronm7/",timeStep,"_",chainNum,"_",args[4],".RDS"))==TRUE){
print("previous file exists, restarting chain")
out<-readRDS(paste0("/home/users/aaronm7/",timeStep,"_",chainNum,"_",args[4],".RDS"))
}

#on JASMIN I found you need to create bayesian setup even if re-starting a chain, I think as this initiates the cluster needed to run in parallel
BS3PGDN <- createBayesianSetup(likelihood = LL_sitka, prior = priorVals, names = nm, parallel = 8, catchDuplicates = F )

#check if a run already exists and if so restart it
out<-suppressWarnings(if(file.exists(paste0("/home/users/aaronm7/",timeStep,"_",chainNum,"_",args[4],".RDS"))==TRUE) runMCMC(bayesianSetup =out, sampler = "DEzs", settings = settings) else runMCMC(bayesianSetup =BS3PGDN, sampler = "DEzs", settings = settings))
#Save output
saveRDS(out,file=paste0(timeStep,"_",chainNum,"_",args[4],".RDS"))

#stop cluster before restarting the loop otherwise JASMIN sometimes throws an error
stopParallel(BS3PGDN)
rm(BS3PGDN)


}
```

# Spatial runs on JASMIN
