library(tidyr)
library(tidyverse)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(dplyr)
library(future)
library(ncdf4) 
library(raster) 
library(ggplot2)
library(httr)
library(furrr)
library(viridis)
library(fr3PGDN)
library(tibble)
library(miscTools)
library(parallel)


dataDir="/gws/nopw/j04/hydro_jules/data/uk/driving_data/chess/chess-met/daily"#Location of input data
outputDir="/home/users/aaronm7/"#output file for bricked rasters
saveFile="/home/users/aaronm7/spatialChunks/"#save file for spatial chunks of data

outputDir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor"
dataDir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\models\\spatial_met_data\\CHESS\\daily"#Location of input data
saveFile="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\Prafor\\spatialChunks\\"#save file for spatial chunks of data




spatSplitX(dataDir=dataDir,outputDir=outputDir,saveFile=saveFile,dates=1971:2018,numChunks=91)
#create chunks of approx 10,000 grid cells in parallel, saves to file
#1:35 is approx scotland
plan(multisession, workers = coreNum - 1)

future_map(c(1:2), ~spatDatUKnc(chunk=.x,outputDir,saveFile),.progress = T)





args=(commandArgs(TRUE))
print(args)
chunk=args[1]

fileLocs<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\spatialChunks\\"
paramsFile<-("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\weekly_3_T.RDS")
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

#Create some random forests to match with the spatial climate data 

simDat$site<-list(data.frame(from=(paste0( round(1971),"-01-01")),to=paste0( round(2016),"-30-12")))


#run grid squares in parallel - as running over grid squares should be very scaleable
#outTemp <-pmap(simDat$site, simDat$clm,as.list(simDat$grid_id), .f=~FR3PG_spat_run(site=..1, clm=..2,grid_id=..3,param_draw=param_draw),.progress = T)

outTemp<-mapply(FR3PG_spat_run, site = simDat$site, clm = simDat$clm,grid_id=as.list(simDat$grid_id),MoreArgs = list(param_draw=param_draw),SIMPLIFY = F)
#outTemp<-do.call(rbind,outTemp)

#bind into a single tibble
out<-as_tibble(data.table::rbindlist(outTemp,fill=T))
#re-add grid_id values
#out$grid_id<-simDat$grid_id
grF<-simDat[,c(1,3,4)]
out<-merge(out,grF,by.x="grid_id",by.y = "grid_id")

saveRDS(out,paste0("/work/scratch-nopw/alm/spatOutput/","SoNWal_",fileName))










#function to find yield class - this has been changed - see updated function in shiny code etc!
YC.finder <- function(HT,AGE=59) 
{
  if(is.na(HT)==F){
  YC.site = ((HT/(1-exp(-0.033329*AGE))^1.821054)-14.856317)/1.425397
  if(YC.site>24) 24
  else if (YC.site<4) 4
  else unlist(sub("\\]","",unlist(strsplit(as.character(cut(YC.site,breaks=c(6,8,10,12,14,16,18,20,22,24),right=T)),split="\\,"))[2]))
} 
else return (NA)
  }

out$yc_value<-(out$yc_value-8)/(24-8)
out$yc_q95<-(out$yc_q95-8)/(24-8)
out$yc_q05<-(out$yc_q05-8)/(24-8)
#
#
out$yc_value<-ifelse(out$yc_value<0,0,out$yc_value)
out$yc_q95<-ifelse(out$yc_q95<0,0,out$yc_q95)
out$yc_q05<-ifelse(out$yc_q05<0,0,out$yc_q05)



#Plot results
g1<-out %>%
  mutate( range = yc_q95 - yc_q05) %>%
  dplyr::select( grid_id, mean = yc_value) %>%
  gather( variable, yc_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = yc_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Suitability")+
  coord_equal() + 
 scale_fill_viridis_c( "",option = "plasma", na.value="gray70",limits=c(0,1),
                       breaks=c(0,0.5,1),labels=c("Unsuitable","Suitable","V. Suitable"))


g2<-out %>%
  mutate( range = GPP_q95 - GPP_q05) %>%
  dplyr::select( grid_id, mean = GPP_value) %>%
  gather( variable, GPP_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = GPP_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("GPP")+
  coord_equal() + 
  scale_fill_viridis_c( expression(paste("GPP [tDM"," ",ha^-1,"]",sep="")),option = "plasma", na.value="gray70")

g3<-out %>%
  mutate( range = NPP_q95 - NPP_q05) %>%
  dplyr::select( grid_id, mean = NPP_value) %>%
  gather( variable, NPP_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = NPP_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("NPP")+
  coord_equal() + 
  scale_fill_viridis_c( expression(paste("NPP [tDM"," ",ha^-1,"]",sep="")),option = "plasma", na.value="gray70",
                        limits=c(1,4))

g4<-out %>%
  mutate( range = NEE_q95 - NEE_q05) %>%
  dplyr::select( grid_id, mean = NEE_value) %>%
  gather( variable, NEE_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = NEE_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("NEE")+
  coord_equal() + 
  scale_fill_viridis_c( expression(paste("NEE [tDM"," ",ha^-1,"]",sep="")),option = "plasma", na.value="gray70")
  
  
  
  
  ggarrange(g1,g2,g3,g4)
  
  
  
  
  
  
  