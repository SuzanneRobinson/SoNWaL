library("tidyr")
library("tidyverse")
library("purrr")
library("BayesianTools")
library("sensitivity")
library("dplyr")
library("future")
library("ncdf4)" )
library("raster)" )
library("ggplot2")
library("httr")
library("furrr")
library("viridis")
library("fr3PGDN")
library("tibble")
library("miscTools")
library("parallel")
library("sf")
library("rgdal")
library(raster)


args=(commandArgs(TRUE))
print(args)
chunk=args[1]


data_dir="/gws/nopw/j04/hydro_jules/data/uk/driving_data/chess/chess-met/daily"#Location of input data
output_dir="/home/users/aaronm7/"#output file for bricked rasters
save_file="/home/users/aaronm7/spatialChunks/"#save file for spatial chunks of data

output_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\"
data_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\models\\spatial_met_data\\CHESSscape\\daily"#Location of input data
save_file="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\spatialChunks\\"#save file for spatial chunks of data

soil_dat<-c("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\NZplus\\spatialData\\soildataNZero2.xlsx",
            #soilWP<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\wp_fao.tif")
            #soilFC<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\fc_fao.tif")
            "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\ks_fao_octop.tif",
            "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\ths_fao_octop.tif") 

fileLocs<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\spatialChunks\\"

paramsFile<-("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_24_T.RDS")

files <- list.files(path = fileLocs, pattern = "\\.RDS$", full.names = TRUE, 
                    recursive = T)

fileName<-sub("\\/.*", "",list.files(path = fileLocs, pattern = "\\.RDS$", full.names = F, 
                                     recursive = T))[chunk]
fileName<-str_sub(fileName,14)

simDat<-readRDS(files[chunk])
clm_df_full<<-getClimDat("weekly")

#Read in MCMC output to get parameter draws - this code needs tidying a bit
param_drawX<-readRDS(paramsFile)
param_drawX<-as_tibble(getSample(param_drawX, start = 100, coda = TRUE, thin = 1,numSamples = 25 )[[1]])#,colMedians(getSample(param_drawX, start = 1000, coda = TRUE, thin = 1 )[[2]])))
param_draw<-as_tibble(1:nrow(param_drawX))
param_draw$pars<- tibble(list(split(param_drawX, 1:NROW(param_drawX))))%>%
  unnest_legacy()
names(param_draw)<-c("mcmc_id","pars")
param_draw<<-param_draw
param_draw$pars<-unname(param_draw$pars$`list(split(param_drawX, 1:NROW(param_drawX)))`)
simDat$site<-list(data.frame(from=(paste0( round(1971),"-01-01")),to=paste0( round(2016),"-30-12")))
simDat<-add_BGS_dat(simDat, soil_dat)

#run grid squares in parallel - as running over grid squares should be very scaleable
#outTemp <-pmap(simDat$site, simDat$clm,as.list(simDat$grid_id), .f=~FR3PG_spat_run(site=..1, clm=..2,grid_id=..3,param_draw=param_draw),.progress = T)
simDat$wp<-simDat$wp/(simDat$soil_depth*10)#?

simDat$wp[is.na(simDat$wp)==T]<-0
simDat$fc[is.na(simDat$fc)==T]<-0
simDat$sp[is.na(simDat$sp)==T]<-0
simDat$soil_cond[is.na(simDat$soil_cond)==T]<-0
simDat$soil_depth[is.na(simDat$soil_depth)==T]<-0

outTemp<-mapply(FR3PG_spat_run, site = simDat$site, clm = simDat$clm,soil=simDat$soilTex,wp=simDat$wp,fc=simDat$fc,sp=simDat$sp,cond=simDat$soilCond,soil_depth=simDat$soil_depth,grid_id=as.list(simDat$grid_id),MoreArgs = list(param_draw=param_draw),SIMPLIFY = F)

#bind into a single tibble
out<-as_tibble(data.table::rbindlist(outTemp,fill=T))
#re-add grid_id values
#out$grid_id<-simDat$grid_id
grF<-simDat[,c(1,2,3,7,8,9)]
out<-merge(out,grF,by.x="grid_id",by.y = "grid_id")

saveRDS(out,paste0("/work/scratch-nopw/alm/spatOutput/","SoNWal_",fileName))


