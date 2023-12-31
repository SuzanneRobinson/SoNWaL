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
library(raster)


args=(commandArgs(TRUE))
print(args)
chunk=args[1]


# define location of data sets
clm_df_reg<-readRDS("/home/users/aaronm7/3pgData/regionalClmDat.RDS")
data_dir="/gws/nopw/j04/hydro_jules/data/uk/driving_data/chess/chess-met/daily"#Location of input data
output_dir="/home/users/aaronm7/"#output file for bricked rasters
save_file="/home/users/aaronm7/spatialChunks/"#save file for spatial chunks of data

# list of soil data locations, _fao files are the map files from https://esdac.jrc.ec.europa.eu/content/maps-indicators-soil-hydraulic-properties-europe
# soildataNZero2 is the data from Astley derived from BGS data
soil_dat<-c(
            "/home/users/aaronm7/3pgData/wp_fao.tif",
            "/home/users/aaronm7/3pgData/fc_fao.tif",
             "/home/users/aaronm7/3pgData/ks_fao_octop.tif",
            "/home/users/aaronm7/3pgData/ths_fao_octop.tif",
            "/home/users/aaronm7/3pgData/soildataNZero2.RDS")
fileLocs<-"/work/scratch-nopw/alm/chessSpatial/"

#location of MCMC results file
paramsFile<-("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\final_results\\sitka\\HD\\weekly\\mcmcOut\\weekly_3_final_T.RDS")

# list spatial climate files
files <- list.files(path = fileLocs, pattern = "\\.RDS$", full.names = TRUE, 
                    recursive = T)

fileName<-sub("\\/.*", "",list.files(path = fileLocs, pattern = "\\.RDS$", full.names = F, 
                                     recursive = T))[chunk]
fileName<-str_sub(fileName,14)

# read in spatial climate file corresponding to spatial chunk you are running (chunk comes from JASMIN array batch file)
simDat<-readRDS(files[chunk])


# read in the mcmc results file to sample from the posterior
param_drawX<-readRDS(paramsFile)
param_drawX<-mergeChains(getSample(param_drawX, start = 100, coda = TRUE, thin = 1,numSamples = 10 ))#,colMedians(getSample(param_drawX, start = 1000, coda = TRUE, thin = 1 )[[2]])))
param_drawX<-as.data.frame(param_drawX)
#param_drawX<-conc_chains(param_drawX,  clm_df_reg=clm_df_reg)
param_draw<-as_tibble(1:nrow(param_drawX))
param_draw$pars<- tibble(list(split(param_drawX, 1:NROW(param_drawX))))%>%
  unnest_legacy()
names(param_draw)<-c("mcmc_id","pars")
param_draw<<-param_draw
param_draw$pars<-unname(param_draw$pars$`list(split(param_drawX, 1:NROW(param_drawX)))`)
#simDat$site<-list(data.frame(from=(paste0( round(1971),"-01-01")),to=paste0( round(2016),"-30-12")))

# add soil data 
simDat<-add_BGS_dat(simDat, soil_dat)
mcmcReg<-readRDS(paramsFile)
# concatonate chains which were fit for individual sites
concChains<-conc_chains(mcmcReg,  clm_df_reg=clm_df_reg)

# predict some soil variables from Astleys data, landsberg and sands data and biosoil survey
simDat<-soil_regr(simDat,concChains)

simDat<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\simDatExamp.RDS")

#add hazard years
hzYrsLoc<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\hazardData\\hazYrs_V2_2051-2080.RDS"
hzYrs<-do.call(rbind,readRDS(hzYrsLoc))%>%nest_by(x,y)
names(hzYrs)<-c("x","y","hzYrs")
simDat<-simDat %>% left_join(hzYrs, by=c("x" = "x", "y" = "y"))



# run spatial sonwal function - note fc, wp are either from astleys data "fc", "wp" or
# from european map data "fc_map", "wp_map", sp is always from european map data
outTemp<-mapply(SoNWaL_spat_run, site = simDat$site, clm = simDat$clm, hzYrs = simDat$hzYrs,
                wp=simDat$wp_map,fc=simDat$fc_map,sp=simDat$sp, plant_year =1961,
                cond=simDat$soil_cond,carbon = simDat$SOC,soil_depth=simDat$soil_depth,N0=simDat$no,C0=simDat$co,
                grid_id=as.list(simDat$grid_id),MoreArgs = list(param_draw=param_draw),SIMPLIFY = F)

#bind into a single tibble
out<-as_tibble(data.table::rbindlist(outTemp,fill=T))

#add some metadata from original file
grF<-simDat[,c(1,3,4,7,8)]
out<-merge(out,grF,by.x="grid_id",by.y = "grid_id")

saveRDS(out,paste0("/work/scratch-nopw/alm/spatOutput/","SoNWal_",fileName))

hzYrs = simDat$hzYrs
site = simDat$site
clm = simDat$clm[[1]]
wp=simDat$wp_map
fc=simDat$fc_map
sp=simDat$sp
plant_year =2020
cond=simDat$soil_cond
carbon = simDat$SOC
soil_depth=simDat$soil_depth
N0=simDat$no
C0=simDat$co
