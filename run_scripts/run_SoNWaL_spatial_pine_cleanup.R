library("lubridate")
library("RNetCDF")
library("zoo")
library("coda")
library("tidyr")
library("tidyverse")
library("purrr")
library("BayesianTools")
library("sensitivity")
library("dplyr")
library("future")
library("ncdf4" )
library("ggplot2")
library("httr")
library("furrr")
library("viridis")
library("SoNWaL")
library("tibble")
library("miscTools")
library("parallel")
library("sf")
library("rgdal")
library("raster")
library("ProfoundData")


args=(commandArgs(TRUE))
#args=c(1:1000)

print(args)
chunk=as.numeric(args[1])
print(chunk)

rcp=60


fileLocs<-paste0("/gws/nopw/j04/uknetzero/aarons_files/chess_manips/spatial_chunks_",rcp,"_01/")

fileName<-sub("\\/.*", "",list.files(path = fileLocs, pattern = "\\.RDS$", full.names = F,
                                     recursive = T))
fileName<-stringr::str_sub(fileName,14)
fileCheckLocs<-paste0("/gws/nopw/j04/uknetzero/suzies_files/chess_manips/spatOutput_pine_",rcp,"_01_unc")
fileCheck <- list.files(path = fileCheckLocs, pattern = "\\.RDS$", full.names = F,
                        recursive = T)
fileCheck<-stringr::str_sub(fileCheck,8)

missedFiles<- fileName[!fileName %in% fileCheck]


timeStep<-"weekly"
scape =T

# Read in and format climate data
climDir<<-("/home/users/suzanne/PRAFOR_3PG/data/")
clm_df_full<-getClimDatX("weekly",climDir)

#print("pine")
setwd("/gws/nopw/j04/uknetzero/suzies_files/spatialRunScriptsPine/")
setDB(db_name ="ProfoundData.sqlite")

# Load clm_df_pine
sites    <- c("peitz","soro" ,"hyytiala","solling_beech","solling_spruce")
##hyytiala
site = 3
clm_df_pine<-getClmPine(timeStep)
print(head(clm_df_pine))


# list of soil data locations, _fao files are the map files from https://esdac.jrc.ec.europa.eu/content/maps-indicators-soil-hydraulic-properties-europe
# soildataNZero2 is the data from Astley derived from BGS data
soil_dat<-c("/gws/nopw/j04/uknetzero/aarons_files/spatialRunScripts/Data/wp_fao.tif",
            "/gws/nopw/j04/uknetzero/aarons_files/spatialRunScripts/Data/fc_fao.tif",
            "/gws/nopw/j04/uknetzero/aarons_files/spatialRunScripts/Data/ks_fao_octop.tif",
            "/gws/nopw/j04/uknetzero/aarons_files/spatialRunScripts/Data/ths_fao_octop.tif",
            "/gws/nopw/j04/uknetzero/aarons_files/spatialRunScripts/Data/soildataNZero2.RDS")




#location of MCMC results file
paramsFile<-("/gws/nopw/j04/uknetzero/suzies_files/spatialRunScriptsPine/pineCals/output_26_08_21/weekly_pine_3_T.RDS")


simDat<-readRDS(paste0("/gws/nopw/j04/uknetzero/aarons_files/chess_manips/spatial_chunks_",rcp,"_01/spatialChunk_",missedFiles[chunk]))

# read in the mcmc results file to sample from the posterior
param_drawX<-readRDS(paramsFile)
param_drawX<-mergeChains(getSample(param_drawX, start = 100, coda = TRUE, thin = 1,numSamples = 500 ))#,colMedians(getSample(param_drawX, start = 1000, coda = TRUE, thin = 1 )[[2]])))

# Modified the conc_chains function script here to sample from the Finish Pine calibration
# This could be removed and the original function called under future spatial Pine runs using regional calibrations
conc_chains<-function(mcmcReg){

        # Concatonate site specific params from Finnish Pine MCMC chains
        #' @param mcmcReg: Finnish MCMC output
        #' @return concChains: concatonate chains

        # concatonate chains
        nm_all<-c("wiltPoint" , "fieldCap" ,    "satPoint"   ,  "K_s",          "V_nr",
        "sigma_zR",     "E_S1"        , "E_S2"  ,       "shared_area" , "maxRootDepth",
        "K_drain" ,     "pFS2"         ,"pFS20" ,       "aS"        ,   "nS",
        "pRx"      ,    "pRn"          ,"gammaFx",      "gammaF0"  ,    "tgammaF",
        "Rttover"   ,   "mF"           ,"mR"      ,     "mS"      ,     "SLA0",
        "SLA1"       ,  "tSLA"         ,"alpha"    ,    "Y"      ,      "m0",
        "MaxCond"     , "LAIgcx"       ,"CoeffCond" ,   "BLcond",       "Nf",
         "Navm"         ,"Navx"        , "klmax"      ,  "krmax",        "komax",
         "hc"          , "qir"        ,  "qil"         , "qh"  ,         "qbc",
        "el"          , "er"        ,   "SWconst0"     ,"SWpower0",     "Qa",
        "Qb"          , "MaxIntcptn")

        mChains<- if (is.null(nrow(mcmcReg))==T) as.data.frame((mcmcReg$chain[[1]])) else as.data.frame(mcmcReg)
        names(mChains)<-if (is.null(nrow(mcmcReg))==T) nm_all else nm_all[-c(200:202)]
        concChains<-mChains %>%
        pivot_longer(cols = starts_with(c("k","startN","startC")),
                 names_to = c(".value", "wpKey"), names_sep = "_Si") %>%
    dplyr::select(-wpKey)

        concChains<- if (is.null(nrow(mcmcReg))==F) distinct(concChains, pFS2, pFS20,gammaF0,tgammaF,Rttover, .keep_all=T ) else concChains

return(concChains)
}

# Return conc_chains for parameters
param_drawX<-conc_chains(param_drawX)

# Selecting the 1st 100/105 rows of params (so consistent with Spruce)
param_drawX<-head(param_drawX, 100)

# Load file with the merged MCMC Spruce parameters
# Assign missing columns startN, startC and K using Spruce parameters
# This is necessary to obtain soil data
paramsFileSpruce<-("/gws/nopw/j04/uknetzero/suzies_files/spatialRunScripts/spatial_spruce_params.RDS")
params_spruce<-readRDS(paramsFileSpruce)

#print(subset(params_spruce, select=c("startC", "startN")))


# Load file with the merged MCMC Spruce parameters
# Assign missing columns startN, startC and K using Spruce parameters and standard parameters
# This assignment is necessary to obtain soil data
paramsFileSpruce<-("/gws/nopw/j04/uknetzero/suzies_files/spatialRunScripts/spatial_spruce_params.RDS")
params_spruce<-readRDS(paramsFileSpruce)
#print(subset(params_spruce, select=c("startC", "startN")))

# Fill startN and startC with values from Spruce
param_drawX$startN <- params_spruce[c('startN')]
param_drawX$startC <- params_spruce[c('startC')]
# k: Extinction coefficient for absorption of PAR by canopy set from Xenakis et al. 2008
#param_drawX$k <- 0.52
# k: Extinction coefficient for absorption of PAR by canopy set by value in getParmsPine
param_drawX$k <- 0.6



param_draw<-as_tibble(1:nrow(param_drawX))
param_draw$pars<- tibble(list(split(param_drawX, 1:NROW(param_drawX))))%>%
  unnest_legacy()
names(param_draw)<-c("mcmc_id","pars")
param_draw<<-param_draw
param_draw$pars<-unname(param_draw$pars$`list(split(param_drawX, 1:NROW(param_drawX)))`)
# Line below is not uncommented in spatial run script
#simDat$site<-list(data.frame(from=(paste0( round(1971),"-01-01")),to=paste0( round(2016),"-30-12")))

# add soil data
simDat<-add_BGS_dat(simDat, soil_dat)

# predict some soil variables from Astleys data, landsberg and sands data and biosoil survey
simDat<-soil_regr(simDat)


simDat$site<-1

#add hazard years
hzYrsLoc<-paste0("/gws/nopw/j04/uknetzero/aarons_files/spatialRunScripts/Data/hazYrs_rcp",rcp,"_2050-2079.rds")
timePer<-as.numeric(sub('-.*', '', sub(paste0('.*hazYrs_rcp",rcp,"_'), '', hzYrsLoc)))-1
hzYrs<-do.call(rbind,readRDS(hzYrsLoc))%>%nest_by(x,y)
names(hzYrs)<-c("x","y","hzYrs")
simDat<-simDat %>% left_join(hzYrs, by=c("x" = "x", "y" = "y"))

# run spatial sonwal function - note fc, wp are either from astleys data "fc", "wp" or
# from european map data "fc_map", "wp_map", sp is always from european map data
outTemp<-mapply(SoNWaL_spat_run_pine, site = simDat$site, clm = simDat$clm,  hzYrs = simDat$hzYrs,
                wp=simDat$wp_map,fc=simDat$fc_map,sp=simDat$sp, plant_year = ifelse(scape==T,2020,1961),
                cond=simDat$soil_cond,carbon = simDat$SOC,soil_depth=simDat$soil_depth,N0=simDat$no,C0=simDat$co,
                grid_id=as.list(simDat$grid_id),MoreArgs = list(param_draw=param_draw, scape=scape),SIMPLIFY = F)

#bind into a single tibble
out<-as_tibble(data.table::rbindlist(outTemp,fill=T))

#add some metadata from original file
grF<-simDat[,c(1,3,4,7,8)]
out<-merge(out,grF,by.x="grid_id",by.y = "grid_id")


saveRDS(out,paste0("/gws/nopw/j04/uknetzero/suzies_files/chess_manips/spatOutput_pine_",rcp,"_01_unc/","SoNWal_",missedFiles[chunk]))
