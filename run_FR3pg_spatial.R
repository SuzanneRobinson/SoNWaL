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

param.drawX<-read_rds('C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\outx.RDS')
param.drawX<-as_tibble(getSample(param.drawX, start = 100, coda = TRUE, thin = 1,numSamples = 1000 )[[1]])#,colMedians(getSample(param.drawX, start = 1000, coda = TRUE, thin = 1 )[[2]])))
param.draw<-as_tibble(1:nrow(param.drawX))
param.draw$pars<- tibble(list(split(param.drawX, 1:NROW(param.drawX))))%>%
  unnest()
names(param.draw)<-c("mcmc_id","pars")
param.draw<<-param.draw


#Create some random forests to match with the spatial climate data (will use actual forest data in future)
genSiteTb<-genRandSiteDat(from=1950,to=2009)

#take a sample of the spatial data (else it'll take a really long time to run not on a cluster)
spatSimDat <- inner_join(genSiteTb, simDat, by = 'grid_id')%>%
sample_n(25000) 
#Try and clump the data a bit?
#spatSimDat <- spatSimDat %>% group_by(ID) %>% sample_n(300)
#spatSimDat<-spatSimDat[spatSimDat$ID %in% round(rnorm(5,50,50)),]
#spatSimDat$dataAv<-sapply(spatSimDat$clm, function(x) sum(x))
#establish cluster
plan(multisession,workers = availableCores())
#run grid squares in parallel - as running over grid squares should be very scaleable
out2 <-future_map2(spatSimDat$site, spatSimDat$clm, ~FR3PG_spat_run(.x, .y,param.draw),.progress = T)
#bind into a single tibble
out<-as_tibble(data.table::rbindlist(out2))
#re-add grid_id values
out$grid_id<-spatSimDat$grid_id
#Remove areas where there is no data (such as the ocean), for speed FR3PG_spat_run skips these and just adds NA's
#out<-drop_na(out)

#Write output to file
#saveRDS(out,"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\scotSpat_isle.rds")
out<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\scotSpat_clumped.rds")


out[,c(3)]<-log(1+10*out[,c(3)])


#Plot results
out %>%
  mutate( range = q95 - q05) %>%
  dplyr::select( grid_id, mean = value, range) %>%
  gather( variable, GPP_value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = GPP_value )) + #*100 to get hectares from 1km grid squares
  geom_raster()+
  facet_wrap(~variable)+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle("Spatial analysis of FR3PGN model - Scotland")+
  coord_equal() + 
 scale_fill_viridis_c( expression(log~Wsbr~(tDM~ha^-1)),option = "plasma", na.value="gray70")


  
  
  
  
  
  
  
  
  
  
  
  