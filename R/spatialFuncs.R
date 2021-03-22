##UK spatial data extraction

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
library(dplyr)
library(lubridate)

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


fx1<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\jansolrad\\w001001.adf")
fx2<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\febsolrad\\w001001.adf")
fx3<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\marsolrad\\w001001.adf")
fx4<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\aprsolrad\\w001001.adf")
fx5<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\maysolrad\\w001001.adf")
fx6<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\junsolrad\\w001001.adf")
fx7<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\julsolrad\\w001001.adf")
fx8<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\augsolrad\\w001001.adf")
fx9<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\sepsolrad\\w001001.adf")
fx10<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\octsolrad\\w001001.adf")
fx11<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\novsolrad\\w001001.adf")
fx12<-raster("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\prevDat\\scotland_spatial_data\\scotland_spatial_data\\Solar radiation\\decsolrad\\w001001.adf")

#hadUKRast <- crop(fj, extent(2.3e+05, 3.5e+5, 6.6e+05, 7.3e+5))
#plot(hadUKRast[[1]])


fj<-brick(fx1,fx2,fx3,fx4,fx5,fx6,fx7,fx8,fx9,fx10,fx11,fx12)
writeRaster(fj, "C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\monthly\\solRad\\pv_hadukgrid_uk_1km_mon_201601-201612.nc", format="CDF",overwrite=TRUE)

#left, ?,?,top

#plot(fj)

for(i in c(1:length(files))){
print(fileNames[i])
#read in files as rasters into list

hadUKRast <- lapply(files[i], function(x) { brick(x) })

#Crop each raster layer (currently approximately around scotland)

hadUKRast <- lapply(hadUKRast, function(x) {crop(x, fj) })
hadUKRast <- lapply(hadUKRast, function(x) {crop(x, extent(2.3e+05, 3.5e+5, 6.6e+05, 7.3e+5)) })

#crs_input<-crs("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +a=6377563.396 +b=6356256.909")
#crs_output<-crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")
#hadUKRast_decLatLong <- inborutils::reproject_coordinates(as.data.frame(coordinates(hadUKRast)),col_long = "x", col_lat = "y", crs_input=crs_input, crs_output=crs_output)

#layer rasters into single raster

hadUKRast<-raster::brick(hadUKRast)

#get climate values from raster 

rasValue=as.data.frame(raster::extract(hadUKRast,extent(2.3e+05, 3.5e+5, 6.6e+05, 7.3e+5),cellnumbers=F)) 

#Transpose data before putting into table
rasValue<-rasValue %>% purrr::transpose()
#Convert transposed data for each cell into a dataframe
colNm<-fileNames[i]
rasValue2<-lapply(rasValue,function(x) setNames(data.frame(unlist(x)),fileNames[i]))

#addSolar<-function(rasValue2,hadUKRast_decLatLong){
#  
#  hadUKRast <- lapply(rasValue2, function(x) {crop(x, extent(0.5e+05, 4.5e+05, 6e+05, 1000000)) })
#  
#  sirad::ap(as.Date("2016-01-16", "%Y-%m-%d"),  -8.043887,58.74152, extraT=NULL, A=NA, B=NA, 145.72)
#  
#  
#}


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
genSite<- list(data.frame(from=(paste0( round(from),"-01-01")),to=paste0( round(to),"-30-12")))
for(i in c(1:nrow(coordinates(hadUKRast)))) {
  genSite[[i]] <- data.frame(from=(paste0( round(from),"-01-01")),to=paste0( to,"-30-12"))
}

genSiteTb <- tibble(grid_id = c(1:nrow(coordinates(hadUKRast))),
                 site = genSite)
return(genSiteTb)
}

##calculate yield class
YCfunc<-function(out){
  MAI<-(aggregate(out$MAI~out$Year,FUN=mean))
  CAI<-(aggregate(out$CAI~out$Year,FUN=mean))
  names(MAI)<-c("x","y")
  names(CAI)<-c("x","y")
  MAIx<-MAI[-c(1:10),]
  CAIx<-CAI[-c(1:10),]
  
  MAI$x<-c(1:nrow(MAI))
  CAI$x<-c(1:nrow(CAI))
  
  CAI$y<-(predict(loess(CAI$y~CAI$x,span=1)))
  MAI[1,2]<-0
  MAI$y<-(predict(loess(MAI$y~MAI$x)))
  
  plot(CAI,col="white")
  lines(MAI,col="blue")
  lines(CAI,col="red")
  
  
  interSec<- tryCatch(
    {
      curve_intersect(MAI[-c(1:10),], CAI[-c(1:10),], empirical = TRUE, domain = NULL)[[2]]
      
    },
    error = function(cond){
      #CAIx[which(diff(sign(diff(CAIx)))==2)+1]
      CAIx$MAI<-MAIx$y
      CAIx[which(CAIx$y<CAIx$MAI),][1,3]
    })
  
  interSec<-if(length(interSec)>1) interSec[2] else interSec
  return(interSec)
  
}


#Spatial run function which takes a single grid cell and associated climate variables, runs simulations for a range of MCMC posterior values
#Calculates the mean and a 95% credible interval
#Currently there is some data manipulation going on to spread a single years data over multiple years until all data is downloaded
#' @param site site data
#' @param clm climate data
FR3PG_spat_run <- function(site, clm,param.draw,clm_df_full){
  library(lubridate)
  
    if(is.na(clm[1,3])==T) return (tibble::as_tibble(data.frame(
                                                                yc_q05=NA,yc_q95=NA,yc_value=NA,
                                                                Wsbr_q05=NA,Wsbr_q95=NA,Wsbr_value=NA,
                                                                height_q05=NA,height_q95=NA, height_value=NA,
                                                                MAI_q05=NA,MAI_q95=NA, MAI_value=NA,
                                                                CAI_q05=NA,CAI_q95=NA, CAI_value=NA,
                                                                GPP_q05=NA,GPP_q95=NA,GPP_value=NA,
                                                                NPP_q05=NA,NPP_q95=NA,NPP_value=NA,
                                                                NEE_q05=NA,NEE_q95=NA,NEE_value=NA,
                                                                Reco_q05=NA,Reco_q95=NA,Reco_value=NA)) ) else
    
    
    FR3pgRun<-function(params){
      #get default parameters
      sitka<- getParms(timeStp = 12)

      #Update parameters with proposals - In this data might be missing parameters which are named differently??
      nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
            "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
            "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
            "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er")
  
      sitka[nm]<-as.data.frame(params[nm])
      clmDat<-as.data.frame(clm[,c(6,5,4,2,3,1)])
      names(clmDat)<-c("Tmin","Tmax","Tmean","Rain","SolarRad","FrostDays")
      clmDat$MonthIrrig<-0
      clmDat$Month<-c(1:12)

      
      sqLength<-lubridate::year(as.Date(site$to,"%Y-%d-%m"))-lubridate::year(as.Date(site$from,"%Y-%d-%m"))+1
      weatherDat<-as.data.frame(purrr::map_dfr(seq_len(sqLength), ~clmDat))
      weatherDat$Year<-rep(lubridate::year(as.Date(site$from,"%Y-%d-%m")):lubridate::year(as.Date(site$to,"%Y-%d-%m")),each=12) 
      weatherDat<-weatherDat[,c(names(clm.df.full[,c(1:9)]))]
      weatherDat$date<-as.Date(paste0(weatherDat$Year,"-01-",weatherDat$Month),"%Y-%d-%m")
      weatherDat$FrostDays<-0
      weatherDat$week<-week(weatherDat$date)
      
      
      
      weatherDat$Rain<-weatherDat$Rain#*0.6
    #  weatherDat$Tmax<-weatherDat$Tmax*1.1
    #  weatherDat$Tmean<-weatherDat$Tmean*1.1
    #  weatherDat$Tmin<-weatherDat$Tmin*1.1
    #  weatherDat$SolarRad<-weatherDat$SolarRad*1.1
      
      sitka$weather<-weatherDat
      
    
      
      
      #run model and output data
      out<-do.call(fr3PGDN,sitka)
     # out<-out[-1,]
     # out$an_dbh<-rep(aggregate(out$dg~out$Year,FUN=max)[,2],each=12)
      # out$an_dbh<-rep(aggregate(out$dg~out$Year,FUN=max)[,2],each=12)
      out$age<-rev(as.numeric(2018-out$Year))
      #out$sVol<-2000*((out$dg/2)^2*pi*out$hdom)/100
      out$MAI<-(out$Vu)/out$age
      out$CAI<-c(rep(0,12),diff(out$Vu,lag=12))

      

      
      out$yc<-YCfunc(out)
      
      
      return(out)
    }
  
  #select stem and branch biomass or other outputs, take mean and 95% credible intervals
    #value to use? Wsbr or GPP etc?
 site_out <- tryCatch(
   {param.draw %>%
   
    dplyr::mutate( sim = map(pars, FR3pgRun)) %>%
    dplyr::select(mcmc_id, sim) %>%
    tidyr::unnest_legacy() %>%
    dplyr::summarise(
      yc_q05 = quantile(yc, 0.05, na.rm = T),
      yc_q95 = quantile(yc, 0.95, na.rm = T),
      yc_value = quantile(yc, 0.5, na.rm = T),
      
      Wsbr_q05 = quantile(dg, 0.05, na.rm = T),
      Wsbr_q95 = quantile(dg, 0.95, na.rm = T),
      Wsbr_value = quantile(dg, 0.5, na.rm = T),
      
      height_q05 = quantile(hdom, 0.05, na.rm = T),
      height_q95 = quantile(hdom, 0.95, na.rm = T),
      height_value = quantile(hdom, 0.5, na.rm = T),
      
      
      MAI_q05 = quantile(MAI, 0.05, na.rm = T),
      MAI_q95 = quantile(MAI, 0.95, na.rm = T),
      MAI_value = quantile(MAI, 0.5, na.rm = T),
      
      
      CAI_q05 = quantile(CAI, 0.05, na.rm = T),
      CAI_q95 = quantile(CAI, 0.95, na.rm = T),
      CAI_value = quantile(CAI, 0.5, na.rm = T),
      
      GPP_q05 = quantile(GPP, 0.05, na.rm = T),
      GPP_q95 = quantile(GPP, 0.95, na.rm = T),
      GPP_value = quantile(GPP, 0.5, na.rm = T),
      
      NPP_q05 = quantile(NPP, 0.05, na.rm = T),
      NPP_q95 = quantile(NPP, 0.95, na.rm = T),
      NPP_value = quantile(NPP, 0.5, na.rm = T),
      
      NEE_q05 = quantile(NEE, 0.05, na.rm = T),
      NEE_q95 = quantile(NEE, 0.95, na.rm = T),
      NEE_value = quantile(NEE, 0.5, na.rm = T),

      Reco_q05 = quantile(Reco, 0.05, na.rm = T),
      Reco_q95 = quantile(Reco, 0.95, na.rm = T),
      Reco_value = quantile(Reco, 0.5, na.rm = T)
      
    )
},
error = function(cond){
  site_out <- tibble::as_tibble(data.frame(Wsbr_q05=NA,Wsbr_q95=NA,Wsbr_value=NA,
                                           GPP_q05=NA,GPP_q95=NA,GPP_value=NA,
                                           NPP_q05=NA,NPP_q95=NA,NPP_value=NA,
                                           NEE_q05=NA,NEE_q95=NA,NEE_value=NA,
                                           Reco_q05=NA,Reco_q95=NA,Reco_value=NA))
})

    
    
  return(site_out)
  
}



years <- 2010:2012

if(Sys.info()[1]=="Windows"){
clm.df.full<-read.csv("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\PRAFOR_3PG\\data\\clm_df_full.csv")
}else
{
clm.df.full<-read.csv("/home/users/aaronm7/3pgData/clm_df_full.csv")
}

#Add date
clm.df.full$date<-as.Date(paste(clm.df.full$Year,"-",clm.df.full$Month,"-01",sep=""))
clm.df.full$week<-week(clm.df.full$date)