library(multidplyr)
library(r3PG)
library(tidyr)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(dplyr)
library(future)


##Using the spatial sweedish data - testing spatial analysis using the FR3pg model and calibrated parameters

#Load in sweedish data - taken from Trotsiuk et al. (2020)
load('C:\\R-Packages\\r3PG\\data\\grid_sim.rda')
load('C:\\R-Packages\\r3PG\\data\\grid_input.rda')
load('C:\\R-Packages\\r3PG\\data\\mcmc.rda')
load('C:\\R-Packages\\r3PG\\data\\solling.rda')


# default parameters for swedish data
par_def <- setNames(param_solling$default, param_solling$param_name)
err_def <- setNames(error_solling$default, error_solling$param_name)

# parameters for calibration and their ranges
param_morris.df <- bind_rows(param_solling, error_solling) %>% filter(!is.na(min))
par_cal_names <- param_morris.df$param_name
par_cal_min <- param_morris.df$min
par_cal_max <- param_morris.df$max
par_cal_best <- param_morris.df$default

par_def.df <- select(param_solling, parameter = param_name, piab = default)

#Take a sample of parameter sets from the MCMC chains to get 95% credible intervals - currently MCMC chains from sweden calibrations, will update for UK
param.draw <- getSample(mcmc_out, numSamples = 10, coda = F, whichParameters = 1:20) %>%
  as.data.frame() %>%
  mutate(mcmc_id = 1:n()) %>%
  nest_legacy(-mcmc_id, .key = 'pars')   %>%
  mutate( 
    pars = map(pars, unlist),
    pars = map(pars, ~tibble::enframe(.x, 'parameter', 'piab')),
    pars = map(pars, ~bind_rows(.x, filter(par_def.df, !parameter %in% par_cal_names))))


#Take a sample of 1x1km grid cells for the spatial runs, need to get data for UK
sim.grid <- inner_join(site.grid, climate.grid, by = 'grid_id')%>%
  sample_n(1000) 


#Spatial run function which takes a single 1x1km grid cell and associated climate variables, runs simulations for a range of MCMC posterior values
#Calculates the average and a 95% credible interval
#Currently there is some data manipulation going on so the sweedish data plays nicely with the FR3PG model
# - in particular the sweedish data just has as set of 12 average monthly climate variables, one for each month. 
FR3PG_spat_run <- function(site, forc){
  FR3pgRun<-function(params){
    
   
  #get default parameters
   sitka<- getParms()
    
    #Update parameters with proposals - In this data might be missing parameters which are named differently??
    sitka[names(sitka) %in% params$parameter]<-params$piab[params$parameter %in%  names(sitka)]
    #convert monthly averages to 45 years of data for sweden
    clmDat<-as.data.frame(forc[,1:6])
    names(clmDat)<-c("Tmin","Tmax","Tmean","Rain","SolarRad","FrostDays")
    clmDat$MonthIrrig<-0
    clmDat$Month<-c(1:12)
    sqLength<-year(as.Date(paste0(site$to,"-01")))-year(as.Date(paste0(site$from,"-01")))+1
    weatherDat<-as.data.frame(purrr::map_dfr(seq_len(sqLength), ~clmDat))
    weatherDat$Year<-rep(year(as.Date(paste0(site$from,"-01"))):year(as.Date(paste0(site$to,"-01"))),each=12) 
    weatherDat<-weatherDat[,c(names(clm.df.full[,c(1:9)]))]
    sitka$weather<-weatherDat
    #run model and output data
    out<-do.call(fr3PGDN,sitka)
    return(out)
  }
  
  #select stem and branch biomass or other outputs, take mean and 95% credible intervals
  site_out <- param.draw %>%
    mutate( sim = map(pars, FR3pgRun)) %>%
    select(mcmc_id, sim) %>%
    unnest_legacy() %>%
    summarise(
      q05 = quantile(Wsbr, 0.05, na.rm = T),
      q95 = quantile(Wsbr, 0.95, na.rm = T),
      value = quantile(Wsbr, 0.5, na.rm = T))
  
  return(site_out)
  
}


#establish cluster
  plan(multisession,workers = availableCores())
  #run grid squares in parallel - as running over grid squares should be very scaleable
  out2 <-future_map2(sim.grid$site, sim.grid$forc, ~FR3PG_spat_run(.x, .y),.progress = T)
  out<-as_tibble(data.table::rbindlist(out2))
  out$grid_id<-sim.grid$grid_id
  
  
  #Plot results
  out %>%
    mutate( range = q95 - q05) %>%
    select( grid_id, mean = value, range) %>%
    gather( variable, value, -grid_id) %>%
    inner_join( ., coord.grid, by = 'grid_id') %>%
    ggplot( aes(x, y, fill = value) ) +
    geom_raster()+
    facet_wrap(~variable)+
    theme_void()+
    coord_equal() +
    scale_fill_distiller( '', palette = 'Spectral', limits = c(0, 250))


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  