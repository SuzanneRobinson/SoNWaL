library(r3PG)
library(plyr)
library(tidyr)
library(tidyverse)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(dplyr)
library(future)
library(ncdf4) 
library(raster) 
library(rgdal) 
library(ggplot2)
library(httr)
library(furrr)
library(viridis)
library(fr3PGDN)
library(tibble)
library(miscTools)


#Get UK spatial data - needs updating to web scraping
simDat<-spatDatUK(dataDir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\monthly")
clm_df_full<<-getClimDat("monthly")

param.drawX<-readRDS('C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\monthly_outx_2021-02-265949.RDS')
param.drawX<-as_tibble(getSample(param.drawX, start = 250000, coda = TRUE, thin = 1,numSamples = 25 )[[1]])#,colMedians(getSample(param.drawX, start = 1000, coda = TRUE, thin = 1 )[[2]])))
param.draw<-as_tibble(1:nrow(param.drawX))
param.draw$pars<- tibble(list(split(param.drawX, 1:NROW(param.drawX))))%>%
  unnest_legacy()
names(param.draw)<-c("mcmc_id","pars")
param.draw<<-param.draw
param.draw$pars<-unname(param.draw$pars$`list(split(param.drawX, 1:NROW(param.drawX)))`)

#Create some random forests to match with the spatial climate data (will use actual forest data in future)
genSiteTb<-genRandSiteDat(from=1950,to=2009)

#take a sample of the spatial data (else it'll take a really long time to run not on a cluster)
spatSimDat <<- inner_join(genSiteTb, simDat, by = 'grid_id')#%>%
#sample_n(500) 
#Try and clump the data a bit?
#spatSimDat <- spatSimDat %>% group_by(ID) %>% sample_n(300)
#spatSimDat<-spatSimDat[spatSimDat$ID %in% round(rnorm(5,50,50)),]
#spatSimDat$dataAv<-sapply(spatSimDat$clm, function(x) sum(x))
#establish cluster
plan(multisession,workers = 7)
#run grid squares in parallel - as running over grid squares should be very scaleable
out2 <-future_map2(spatSimDat$site, spatSimDat$clm, ~FR3PG_spat_run(site=.x, clm=.y,param.draw=param.draw,clm_df_full=clm_df_full),.progress = T)
#bind into a single tibble
out<-as_tibble(data.table::rbindlist(out2,fill=T))
#re-add grid_id values
out$grid_id<-spatSimDat$grid_id
#Remove areas where there is no data (such as the ocean), for speed FR3PG_spat_run skips these and just adds NA's
#out<-drop_na(out)

#Write output to file
#saveRDS(out,"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\scotSpat_isle.rds")
#out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\scotSpat_clumped.rds")

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
##out<-na.omit(out)
#out$YC_value<-as.numeric(sapply(out$height_value,YC.finder))
#out$YC_q95<-as.numeric(sapply(out$height_q95,YC.finder))
#out$YC_q05<-as.numeric(sapply(out$height_q05,YC.finder))
#
out$yc_value<-(out$yc_value-8)/(24-8)
out$yc_q95<-(out$yc_q95-8)/(24-8)
out$yc_q05<-(out$yc_q05-8)/(24-8)
#
#
out$yc_value<-ifelse(out$yc_value<0,0,out$yc_value)
out$yc_q95<-ifelse(out$yc_q95<0,0,out$yc_q95)
out$yc_q05<-ifelse(out$yc_q05<0,0,out$yc_q05)

## out$an_dbh<-rep(aggregate(out$dg~out$Year,FUN=max)[,2],each=12)
#out$age<-rev(as.numeric(2018-out$Year))
##out$sVol<-2000*((out$dg/2)^2*pi*out$hdom)/100
#out$MAI<-(out$Vu)/out$age
#out$CAI<-c(rep(0,12),diff(out$Vu,lag=12))
#

#saveRDS(out2,file="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\map_NoDrought.RDS")

out2<-readRDS(out2,file="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\output\\map_Drought.RDS")

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
  
  
  
  
  
  
  