library(multidplyr)
library(r3PG)
library(tidyr)
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

#Get UK spatial data - needs updating to web scraping
simDat<-spatDatUK(dataDir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\monthly")

#Load in Swiss MCMC data to get a range of parameters - taken from Trotsiuk et al. (2020)
param.draw<-swissParams('C:\\R-Packages\\r3PG\\data\\solling.rda',numSamps=10)

#Create some random forests to match with the spatial climate data (will use actual forest data in future)
genSiteTb<-genRandSiteDat(from="1950-01",to="2009-01")

#take a sample of the spatial data (else it'll take a really long time to run not on a cluster)
spatSimDat <- inner_join(genSiteTb, simDat, by = 'grid_id')%>%
  sample_n(5000) 


#establish cluster
plan(multisession,workers = availableCores())
#run grid squares in parallel - as running over grid squares should be very scaleable
out2 <-future_map2(spatSimDat$site, spatSimDat$clm, ~FR3PG_spat_run(.x, .y),.progress = T)
#bind into a single tibble
out<-as_tibble(data.table::rbindlist(out2))
#re-add grid_id values
out$grid_id<-spatSimDat$grid_id
#Remove areas where there is no data (such as the ocean), for speed FR3PG_spat_run skips these and just adds NA's
out<-drop_na(out)

#Write output to file
saveRDS(out,"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\scotSpat.rds")

#Plot results
out %>%
  mutate( range = q95 - q05) %>%
  dplyr::select( grid_id, mean = value, range) %>%
  gather( variable, value, -grid_id) %>%
  inner_join( ., simDat, by = 'grid_id') %>%
  ggplot( aes(x, y, fill = value) ) +
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
  ggtitle("Spatial analysis of FR3PG model - Scotland")+
  coord_equal() +
  scale_fill_viridis_c( '', limits = c(0, 250),option = "B")


  
  
  
  
  
  
  
  
  
  
  
  
  