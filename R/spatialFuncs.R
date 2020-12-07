##UK spatial data extraction

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

##Need to fix web scraping of data!## Current issues with authentication cert
#filename<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\rainfall\\rainfall_hadukgrid_uk_1km_day_20111201-20111231.nc"
#fileUrl<-"http://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.2.1/1km/rainfall/day/v20200731/rainfall_hadukgrid_uk_1km_day_20070401-20070430.nc"
#certFile<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\creds.pem"
#
#options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))
#
#GET(fileUrl, config(cainfo=credtoken),
#    write_disk(filename, overwrite = T), timeout(60))
#
#x <- getURL(fileUrl, cainfo = system.file("CurlSSL", certFile, package = "RCurl"))
#
##make a list of dates to download files
#datseq <- function(t1, t2) { 
#  format(seq(as.Date(t1, "%Y%m%d"), 
#             as.Date(t2, "%Y%m%d"),by="month"), 
#         "%Y%m") 
#}
#
#fileDats<-datseq(as.Date(paste0("2015","01","01"), "%Y%m%d"),as.Date(paste0("2018","12","01"), "%Y%m%d"))
#urls<-paste0()
#
#Map(function(u, d) download.file(u, d, mode="wb"), urls, destinations)


##read in list of spatial climate file names##
#'@param dataDir directory which stores spatial HadUK climate data from ceda
#'@return tibble of site id key with associated dataframe of longitudinal climate data
spatDatUK<-function(dataDir){
  
files <- list.files(path = dataDir, pattern = "\\.nc$", full.names = TRUE, 
                    recursive = T)
fileNames <- sub("\\/.*", "",list.files(path = dataDir, pattern = "\\.nc$", full.names = F, 
                    recursive = T))


for(i in c(1:length(files))){
print(fileNames[i])
#read in files as rasters into list
hadUKRast <- lapply(files[i], function(x) { brick(x) })
#Crop each raster layer (currently approximately around scotland)
hadUKRast <- lapply(hadUKRast, function(x) {crop(x, extent(0.5e+05, 4.5e+05, 6e+05, 1000000)) })
#layer rasters into single raster
hadUKRast<-brick(hadUKRast)
#get climate values from raster 
rasValue=as.data.frame(raster::extract(hadUKRast,extent(0.5e+05, 4.5e+05, 6e+05, 1000000),cellnumbers=F))
#Transpose data before putting into table
rasValue<-rasValue %>% purrr::transpose()
#Convert transposed data for each cell into a dataframe
colNm<-fileNames[i]
rasValue2<-lapply(rasValue,function(x) setNames(data.frame(unlist(x)),fileNames[i]))
#add to tibble
if(i==1){
simDat <- tibble(id = c(1:nrow(coordinates(hadUKRast))),
             data = rasValue2)
} 
else {
simDat$data<-Map(cbind,simDat$data,rasValue2)
}

}
names(simDat)<-c("grid_id","clm")

#Get spatial coordinate data from rasters for plotting
simDat$x<-coordinates(hadUKRast)[,1]
simDat$y<-coordinates(hadUKRast)[,2]

hadUKRast<<-hadUKRast

return(simDat)
}


##create some random generic site data for each gridcell -obvs need to get some real data
#'@param from takes a random date between this value and +10 years as plantation establishment
#'@param to takes a random date between this value and +10 years for the end of plantation growth
genRandSiteDat<-function(from,to){
genSite<- list(data.frame(from=(paste0( round(runif(1, from, from+10)),"-01-01")),to=paste0( round(runif(1, to, to+10)),"-30-12")))
for(i in c(1:nrow(coordinates(hadUKRast)))) {
  genSite[[i]] <- data.frame(from=(paste0( round(runif(1, from, from+10)),"-01-01")),to=paste0( round(runif(1, to, to+10)),"-30-12"))
}

genSiteTb <- tibble(grid_id = c(1:nrow(coordinates(hadUKRast))),
                 site = genSite)
return(genSiteTb)
}



#Spatial run function which takes a single grid cell and associated climate variables, runs simulations for a range of MCMC posterior values
#Calculates the mean and a 95% credible interval
#Currently there is some data manipulation going on to spread a single years data over multiple years until all data is downloaded
#' @param site site data
#' @param clm climate data
FR3PG_spat_run <- function(site, clm,param.draw){
  
  if(is.na(clm$tasmax[1])==T) return (tibble::as_tibble(data.frame(q05=NA,q95=NA,value=NA)) ) else
    
    
    FR3pgRun<-function(params){
      
      #get default parameters
      sitka<- getParms()
      
      #Update parameters with proposals - In this data might be missing parameters which are named differently??
      sitka[names(sitka) %in% names(params)]<-params[names(params) %in%  names(sitka)]
      #convert monthly averages to 45 years of data for sweden
      clmDat<-as.data.frame(clm[,c(6,5,4,2,3,1)])
      names(clmDat)<-c("Tmin","Tmax","Tmean","Rain","SolarRad","FrostDays")
      clmDat$MonthIrrig<-0
      clmDat$Month<-c(1:12)
      
      ####################needs better solution##########################
      #####if using monthly data from CEDA, split into daily values######
      ######else with hydro submodels the land gets scorched :( #########
      clmDat$SolarRad<-clmDat$SolarRad/days_in_month(clmDat$Month)#######
      ###################################################################
      
      sqLength<-year(as.Date(site$to,"%Y-%d-%m"))-year(as.Date(site$from,"%Y-%d-%m"))+1
      weatherDat<-as.data.frame(purrr::map_dfr(seq_len(sqLength), ~clmDat))
      weatherDat$Year<-rep(year(as.Date(site$from,"%Y-%d-%m")):year(as.Date(site$to,"%Y-%d-%m")),each=12) 
      weatherDat<-weatherDat[,c(names(clm.df.full[,c(1:9)]))]
      sitka$weather<-weatherDat
      #run model and output data
      out<-do.call(fr3PGDN,sitka)
      return(out)
    }
  
  #select stem and branch biomass or other outputs, take mean and 95% credible intervals
  site_out <- param.draw %>%
    dplyr::mutate( sim = map(pars, FR3pgRun)) %>%
    dplyr::select(mcmc_id, sim) %>%
    tidyr::unnest_legacy() %>%
    dplyr::summarise(
      q05 = quantile(Wsbr, 0.05, na.rm = T),
      q95 = quantile(Wsbr, 0.95, na.rm = T),
      value = quantile(Wsbr, 0.5, na.rm = T))
  
  return(site_out)
  
}




###Load in swiss MCMC data to get some parameter values (will replace with UK once we can run the full calibrations)
swissParams<-function(MCMCfile,numSamps=10){
  load(MCMCfile)
  # default parameters for swedish data
  par_def <- setNames(param_solling$default, param_solling$param_name)
  err_def <- setNames(error_solling$default, error_solling$param_name)
  
  # parameters for calibration and their ranges
  param_morris.df <- bind_rows(param_solling, error_solling) %>% filter(!is.na(min))
  par_cal_names <- param_morris.df$param_name
  par_cal_min <- param_morris.df$min
  par_cal_max <- param_morris.df$max
  par_cal_best <- param_morris.df$default
  
  par_def.df <- dplyr::select(param_solling, parameter = param_name, piab = default)
  
  load('C:\\R-Packages\\r3PG\\data\\mcmc.rda')
  #Take a sample of parameter sets from the MCMC chains to get 95% credible intervals - currently MCMC chains from sweden calibrations, will update for UK
  param.draw <- getSample(mcmc_out, numSamples = numSamps, coda = F, whichParameters = 1:20) %>%
    as.data.frame() %>%
    dplyr::mutate(mcmc_id = 1:n()) %>%
    nest_legacy(-mcmc_id, .key = 'pars')   %>%
    mutate( 
      pars = map(pars, unlist),
      pars = map(pars, ~tibble::enframe(.x, 'parameter', 'piab')),
      pars = map(pars, ~bind_rows(.x, filter(par_def.df, !parameter %in% par_cal_names))))
  
  return(param.draw)
  
}



