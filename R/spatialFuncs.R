##UK spatial data extraction
library(tidyr)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(future)
library(ncdf4) 
library(raster) 
library(ggplot2)
library(httr)
library(furrr)
library(dplyr)
library(lubridate)
library(parallel)
library(stringr)


##spatial splitting function - split spatial data into chunks
#'@param data_dir directory which stores the CHESS spatial data (under a filename called CHESS) - also should include a folder for the full time series to be written into called fullTS
#'@param save_file where to save the output files
#'@param startDate the year to start the processing from, CHESS goes back to 1961 but you may not need that far back (significantly increases time to run function the further back you go)
#'@param variable whether to process all CHESS climate variables or select individual ones to process (e.g. precip, temp etc. )
#'@param numChunks function creates chunks of 10,000 grid cells going down the UK, 1:39 covers scotland, higher values cover the rest of UK
#'@return climate data with coordinates
spat_split_X <-function(data_dir, 
                        save_file, 
                        output_dir, 
                        dates=c(1960:1990),
                        variable="all", 
                        num_Chunks=c(1:10)) {
  files <- list.files(path = data_dir, 
                      pattern = "\\.nc$", 
                      full.names = TRUE,
                      recursive = T)
  file_Names <- sub("\\/.*", "",list.files(path = data_dir, 
                                          pattern = "\\.nc$", 
                                          full.names = F, 
                                          recursive = T))
  
  #only run merging once as it's pretty slow, 
  #hopefully once done this wont need to be done again
  #Merge monthly files into longer time-series

    file_names <- file_names[which(str_sub(file_names, -11, -8) %in% dates) ]
    
    variable_names <-
      if (variable == "all")
        unique(str_match(file_names, "met_\\s*(.*?)\\s*_gb")[, 2]) else
      variable
    print(variable_names)
    print("merging layers to create single files with full climate time-series - may take a while :)")
    for (i in 1:unique(length(unique(variable_names)))) {
      print(unique(variable_names)[i])
      ifelse(!dir.exists(file.path(output_dir, "/fullTS")), 
             dir.create(file.path(output_dir, "fullTS")), 
             FALSE)
      filesTmp <-
        paste0(data_dir, "/", 
               file_names[grepl(variable_names[i], 
                                file_names) == T])
      splitter <- function(j) {
        topRow <- ifelse(j == 1, 1, (j - 1) * 6)
        rastLayer <- lapply(filesTmp, function(x) {
          brick(x)
        })
        
        rastLayer <-
          lapply(rastLayer, function(x)
            raster::crop(x, extent(x, topRow, j * 6, 1, 656))) #kershope extent 475, 476, 352, 353
        
        rastLayer <- raster::brick(rastLayer)
        rastLayerX <- list(getValues(rastLayer))
        rastLayerX$coords <- coordinates(rastLayer)
        
        saveRDS(rastLayerX,
                paste0(
                  output_dir,
                  "/fullTS/",
                  unique(variable_names)[i],
                  "_",
                  j,
                  ".RDS"
                ))
    }
      
      coreNum <- detectCores()
      
      if (coreNum > 1) {
        plan(multisession, workers = coreNum - 1)
        future_map(numChunks, ~ splitter(j = .x), 
                   .progress = T)
      } else {
        map(numChunks, ~ splitter(j = .x), 
            .progress = T)
    }
    
    
  }
  

}



##read in list of spatial climate data with grid coordinates from spatSplit function output and merge all climate vars into one table associated with each grid cell ##
#'@param chunk chunk number, 1:39 covers whole of scotland, can only go as high as files available from spatSplit
#'@param output_dir names of directory where files (outputs from spatSplit) to go merge are located
#'@param save_file location to save merged files
#'@return tibble of site id key with associated dataframe of longitudinal climate data
spat_dat_UKnc <- function(chunk = 1, output_dir,save_file) {
  library(stringr)
  library(dplyr)
  library(raster)
  
  files <-
    list.files(
      path = paste0(output_dir, "fullTS"),
      pattern = paste0("_", chunk, "\\.RDS$"),
      full.names = TRUE,
      recursive = T
    )

    #read in files as rasters into list
  mapFile <- lapply(files, function(x) {
    (readRDS(x))
  })
  
  #get file names
  file_names <- unique(str_match(files, "TS/\\s*(.*?)\\s*_")[, 2])
  #get names for ID vals
  gridIDNames <- unique(str_match(files, "TS/\\s*(.*?)\\s*.RDS")[, 2])
  
  
  #get climate values from rasters
  for (i in c(1:length(files))) {
    print(i)
    rasValue = as.data.frame(mapFile[[i]][[1]])
    #Transpose data before putting into table
    rasValue <- rasValue %>% purrr::transpose()
    #Convert transposed data for each cell into a dataframe
    colNm <- file_names[i]
    rasValue2 <-
      lapply(rasValue, function(x)
        setNames(data.frame(unlist(x)), unique(file_names)[i]))
    
    #add to tibble and give unique grid cell ID (grid ID poss redundant)
    if (i == 1) {
      simDat <- tibble(id = c(paste0(gridIDNames[i],"_",1:nrow(mapFile[[i]]$coords))),
                       data = rasValue2)
    }  else {
      simDat$data <- Map(cbind, simDat$data, rasValue2)
    }
    
  }
  names(simDat) <- c("grid_id", "clm")
  #Get spatial coordinate data from rasters
  simDat$x <- mapFile[[i]]$coords[, 1]
  simDat$y <- mapFile[[i]]$coords[, 2]
  
  #split into multiple files to avoid ultra large single files (could split by coordinates...probs still too large splitting by lat)
  simDatSp<-split(simDat, (seq(nrow(simDat))-1) %/% 45) 
  
  for(i in c(1:length(simDatSp))){
    saveRDS( simDatSp[[i]], paste0(save_file, "spatialChunk_",chunk,"_", i+1, ".RDS")) 
  }
  
}


##create some random generic site data for each gridcell -obvs need to get some "real" data
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



calcRH<-function(Tmean,Tref=273.16,p,q){
  TmeanK<-Tmean+273.15
  
  RH<- 0.263*p*q*(exp((17.67*(TmeanK-Tref))/(TmeanK-29.65)))^-1
  RH[RH>100]<-100
  RH[RH<0]<-0
  
  return(RH)
}


#Spatial run function which takes a single grid cell and associated climate variables, runs simulations for a range of MCMC posterior values
#Calculates the mean and a 95% credible interval
#Currently there is some data manipulation going on to spread a single years data over multiple years until all data is downloaded
#' @param site site data
#' @param clm climate data
#' @param param_draw parameter draws
FR3PG_spat_run <- function(site, clm,param_draw,grid_id,soil,soilDepth,wp,fc,sp,cond){
  library(lubridate)

  if(soilDepth==1){param_draw$pars <- lapply(param_draw$pars, function(x) {
    x$V_nr<- 0.25
    x$maxRootDepth<- 0.25
    return(x)})}
  if(soilDepth==2){param_draw$pars <- lapply(param_draw$pars, function(x) {
    x$V_nr<- 0.5
    x$maxRootDepth<- 0.5
    return(x)})}
  if(soilDepth==3){param_draw$pars <- lapply(param_draw$pars, function(x) {
    x$V_nr<- 0.8
    x$maxRootDepth<- 0.8
    return(x)})}
  if(soilDepth==4){param_draw$pars <- lapply(param_draw$pars, function(x) {
    x$V_nr<- 1
    x$maxRootDepth<- 1
    return(x)})}
  if(soilDepth==5){param_draw$pars <- lapply(param_draw$pars, function(x) {
    x$V_nr<- 1.5
    x$maxRootDepth<- 1.5
    return(x)})}
  
  if(soil!=0){
    param_draw$pars <- lapply(param_draw$pars, function(x) {
  x$wiltPoint<-wp
  x$fieldCap<-fc
  x$satPoint<-sp
  x$K_s<-cond
  
   return(x)})}
  
  if(soil==1) {param_draw$pars <- lapply(param_draw$pars, function(x) {
  x$E_S1<-0.05
  x$E_S2<-0.3
  x$K_drain<-0.66
  x$SWpower0<-8
  x$SWconst0<-0.65
  return(x)})}
  
  if(soil==0) {param_draw$pars <- lapply(param_draw$pars, function(x) {
  x$E_S1<-0.05
  x$E_S2<-0.3
  x$K_drain<-0.66
  x$SWpower0<-8
  x$SWconst0<-0.65
  return(x)})}
  
  if(is.na(soil)==T) {param_draw$pars <- lapply(param_draw$pars, function(x) {
  x$E_S1<-0.05
  x$E_S2<-0.3
  x$K_drain<-0.66
  x$SWpower0<-8
  x$SWconst0<-0.65
  return(x)})}
  
  if(soil==4) {param_draw$pars <- lapply(param_draw$pars, function(x) {
  x$E_S1<-0.1
  x$E_S2<-0.3
  x$K_drain<-0.5
  x$SWpower0<-5
  x$SWconst0<-0.5
  return(x)})}
  
  
  if(soil==9) {param_draw$pars <- lapply(param_draw$pars, function(x) {
  x$E_S1<-0.3
  x$E_S2<-0.6
  x$K_drain<-1
  x$SWpower0<-5
  x$SWconst0<-0.5
  return(x)})}
  
  
    if(is.na(clm[1,3])==T) return (tibble::as_tibble(data.frame(Year=NA,grid_id=grid_id,Wsbr_q05=NA,Wsbr_q95=NA,Wsbr_value=NA,
                                                                Rs_q05=NA,Rs_q95=NA,Rs_value=NA,
                                                                EvapTransp_q05=NA,EvapTransp_q95=NA,EvapTransp_value=NA,
                                                                volSWC_rz_q05=NA,volSWC_rz_q95=NA,volSWC_rz_value=NA,
                                                                yc_q05=NA,yc_q95=NA,yc_value=NA,
                                                                GPP_q05=NA,GPP_q95=NA,GPP_value=NA,
                                                                NPP_q05=NA,NPP_q95=NA,NPP_value=NA,
                                                                NEE_q05=NA,NEE_q95=NA,NEE_value=NA,
                                                                Reco_q05=NA,Reco_q95=NA,Reco_value=NA,
                                                                LAI_q05=NA,LAI_q95=NA,LAI_value=NA)) ) else
    
    
    FR3pgRun<-function(params){
      #get default parameters
      baseParms<- getParms(timeStp = 52, waterBalanceSubMods = T)

      #Update parameters with proposals 
      nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
            "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
            "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
            "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0","SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
      
      baseParms[nm]<-as.data.frame(params[nm])

      
      
      names(clm)<-c("dailyTmpRange","specHumid","Rain","surfPressure","SolarRadLW","SolarRad","wind","Tmean")
      clm$Tmean<-clm$Tmean-273.15 #convert to celsius from kelvin
      clm<-clm[1:(nrow(clm)-396),]
      #get min and max temp
      clm<-clm%>%
        mutate(Tmin=Tmean-(dailyTmpRange/2),Tmax=Tmean+(dailyTmpRange/2))
      #climate
      clm$RH<-calcRH(Tmean=clm$Tmean,Tref=273.16,p=clm$surfPressure,q=clm$specHumid)
      clm<-rownames_to_column(clm, var="Date")
      clm$Date<-str_sub(clm$Date, -10)
      clm$Date<-as.Date(clm$Date,"%Y.%m.%d")
      clm$Month<-month(clm$Date)
      clm$week<-week(clm$Date)
      clm$Year<-year(clm$Date)
      
      clm$VPD <-((( (0.61078 * exp(17.269 * (clm$Tmean)/(237.3 + clm$Tmean)))*(1-(clm$RH)/100))))
      

      clm<-clm%>%
        dplyr::group_by(Year,week)%>%
        dplyr::summarise(Month=median(Month),Tmax=max(Tmax),Tmin=min(Tmin),Tmean=mean(Tmean),Rain=sum(Rain*86400),SolarRad=mean((SolarRad*86400)/1e+6),FrostDays=0,MonthIrrig=0,VPD=mean(VPD))
      
      baseParms$weather<-as.data.frame(clm)
      
      #run model and output data
      out<-do.call(fr3PGDN,baseParms)
      out$age<-rev(as.numeric(max(out$Year,na.rm=T)-out$Year))
      out$MAI<-(out$Vu)/out$age
      out$CAI<-c(rep(0,52),diff(out$Vu,lag=52))

      out$yc<-YCfunc(out)
      

      
      return(out)
    }
  
  #select stem and branch biomass or other outputs, take mean and 95% credible intervals
    #value to use? Wsbr or GPP etc?
 site_out <- tryCatch(
   {param_draw %>%
    dplyr::mutate( sim = map(pars, FR3pgRun)) %>%
    dplyr::select(mcmc_id, sim) %>%
    tidyr::unnest_legacy() %>%
       group_by(Year,mcmc_id)%>%
       dplyr::summarise(
         grid_id=grid_id,
         dg = mean(dg,  na.rm = T),
         Rs = mean(Rs, na.rm = T),
         EvapTransp = mean(EvapTransp, na.rm = T),
         volSWC_rz = mean(volSWC_rz, na.rm = T),
         yc = mean(yc, na.rm = T),
         GPPsum = sum(GPP,  na.rm = T),
         NPPsum = sum(NPP, na.rm = T),
         NEEsum = sum(NEE, na.rm = T),    
         GPP = mean(GPP,  na.rm = T),
         NPP = mean(NPP, na.rm = T),
         NEE = mean(NEE, na.rm = T),
         Reco = mean(Reco, na.rm = T),
         LAI = mean (LAI, na.rm = T)
         
       )%>%    
       group_by(Year)%>%
         dplyr::summarise(
           grid_id=grid_id,
      Wsbr_q05 = quantile(dg, 0.05, na.rm = T),
      Wsbr_q95 = quantile(dg, 0.95, na.rm = T),
      Wsbr_value = quantile(dg, 0.5, na.rm = T),
      
      Rs_q05 = quantile(Rs, 0.05, na.rm = T),
      Rs_q95 = quantile(Rs, 0.95, na.rm = T),
      Rs_value = quantile(Rs, 0.5, na.rm = T),

      EvapTransp_q05 = quantile(EvapTransp, 0.05, na.rm = T),
      EvapTransp_q95 = quantile(EvapTransp, 0.95, na.rm = T),
      EvapTransp_value = quantile(EvapTransp, 0.5, na.rm = T),
      
      volSWC_rz_q05 = quantile(volSWC_rz, 0.05, na.rm = T),
      volSWC_rz_q95 = quantile(volSWC_rz, 0.95, na.rm = T),
      volSWC_rz_value = quantile(volSWC_rz, 0.5, na.rm = T),
      
      yc_q05 = quantile(yc, 0.05, na.rm = T),
      yc_q95 = quantile(yc, 0.95, na.rm = T),
      yc_value = quantile(yc, 0.5, na.rm = T),
      
      
      GPP_q05 = quantile(GPP, 0.05, na.rm = T),
      GPP_q95 = quantile(GPP, 0.95, na.rm = T),
      GPP_value = quantile(GPP, 0.5, na.rm = T),
      
      NPP_q05 = quantile(NPP, 0.05, na.rm = T),
      NPP_q95 = quantile(NPP, 0.95, na.rm = T),
      NPP_value = quantile(NPP, 0.5, na.rm = T),
      
      NEE_q05 = quantile(NEE, 0.05, na.rm = T),
      NEE_q95 = quantile(NEE, 0.95, na.rm = T),
      NEE_value = quantile(NEE, 0.5, na.rm = T),
      
      
      GPPsum_q05 = quantile(GPPsum, 0.05, na.rm = T),
      GPPsum_q95 = quantile(GPPsum, 0.95, na.rm = T),
      GPPsum_value = quantile(GPPsum, 0.5, na.rm = T),
      
      NPPsum_q05 = quantile(NPPsum, 0.05, na.rm = T),
      NPPsum_q95 = quantile(NPPsum, 0.95, na.rm = T),
      NPPsum_value = quantile(NPPsum, 0.5, na.rm = T),
      
      NEEsum_q05 = quantile(NEEsum, 0.05, na.rm = T),
      NEEsum_q95 = quantile(NEEsum, 0.95, na.rm = T),
      NEEsum_value = quantile(NEEsum, 0.5, na.rm = T),

      Reco_q05 = quantile(Reco, 0.05, na.rm = T),
      Reco_q95 = quantile(Reco, 0.95, na.rm = T),
      Reco_value = quantile(Reco, 0.5, na.rm = T),
      
      LAI_q05 = quantile(LAI, 0.05, na.rm = T),
      LAI_q95 = quantile(LAI, 0.95, na.rm = T),
      LAI_value = quantile(LAI, 0.5, na.rm = T)
      
    )
},
#add na values where there is no data, in the sea etc. 
error = function(cond){
  site_out <- tibble::as_tibble(data.frame(Year=NA,grid_id=grid_id,Wsbr_q05=NA,Wsbr_q95=NA,Wsbr_value=NA,
                                           Rs_q05=NA,Rs_q95=NA,Rs_value=NA,
                                           EvapTransp_q05=NA,EvapTransp_q95=NA,EvapTransp_value=NA,
                                           volSWC_rz_q05=NA,volSWC_rz_q95=NA,volSWC_rz_value=NA,
                                           yc_q05=NA,yc_q95=NA,yc_value=NA,
                                           GPP_q05=NA,GPP_q95=NA,GPP_value=NA,
                                           NPP_q05=NA,NPP_q95=NA,NPP_value=NA,
                                           NEE_q05=NA,NEE_q95=NA,NEE_value=NA,
                                           Reco_q05=NA,Reco_q95=NA,Reco_value=NA,
                                           LAI_q05=NA,LAI_q95=NA,LAI_value=NA))
})

    
    site_out<-site_out[!duplicated(site_out),]
  return(as.data.frame(site_out))
  
}





##kershop forest data extraction
#
#for (i in 2:unique(length(unique(variable_names)))) {
# print(unique(variable_names)[i])
# ifelse(!dir.exists(file.path(output_dir, "/fullTS2")), dir.create(file.path(output_dir, "fullTS2")), FALSE)
# filesTmp <-
#     paste0(data_dir, "/", file_names[grepl(variable_names[i], file_names) == T])
# splitter <- function(j) {
#     topRow <- ifelse(j == 1, 1, (j - 1) * 6)
#     rastLayer <- lapply(filesTmp, function(x) {
#         brick(x)
#       })
#     rastLayer <-
#         lapply(rastLayer, function(x)
#             raster::crop(x, extent(x, 475, 476, 352, 353))) #kershope extent 475, 476, 352, 353
#
#
#
#       rastLayer <- raster::brick(rastLayer)
#       rastLayerX <- list(getValues(rastLayer))
#       rastLayerX$coords <- coordinates(rastLayer)
#
#         saveRDS(rastLayerX,
#            paste0(
#                output_dir,
#                "/fullTS2/",
#                unique(variable_names)[i],
#                "_",
#                j,
#                ".RDS"
#              ))
#    }
#
#    coreNum <- detectCores()
#
#      coreNum <- detectCores()
#
#    if (coreNum > 1) {
#        plan(multisession, workers = coreNum - 1)
#        future_map(1, ~ splitter(j = .x), .progress = T)
#      } else {
#        map(1, ~ splitter(j = .x), .progress = T
#                     )
#
#           }
#
#
#  }
#
#
#
#
#
#
###read in list of spatial climate data with grid coordinates from spatSplit function output and merge all climate vars into one table associated with each grid cell ##
##'@param chunk chunk number, 1:39 covers whole of scotland, can only go as high as files available from spatSplit
##'@param output_dir names of directory where files (outputs from spatSplit) to go merge are located
##'@param save_file location to save merged files
##'@return tibble of site id key with associated dataframe of longitudinal climate data
spatDatUKnc2 <- function(chunk = 1, output_dir,save_file) {
  library(stringr)
  library(dplyr)
  library(raster)
  
  files <-
    list.files(
      path = paste0(output_dir, "/fullTS2"),
      pattern = paste0("_", chunk, "\\.RDS$"),
      full.names = TRUE,
      recursive = T
    )
  
  #read in files as rasters into list
  mapFile <- lapply(files, function(x) {
    (readRDS(x))
  })
  
  file_names <- unique(str_match(files, "TS2/\\s*(.*?)\\s*_")[, 2])
  
  #get climate values from raster
  for (i in c(1:length(files))) {
    rasValue = as.data.frame(mapFile[[i]][[1]])
    
    
    #Transpose data before putting into table
    rasValue <- rasValue %>% purrr::transpose()
    #Convert transposed data for each cell into a dataframe
    colNm <- file_names[i]
    rasValue2 <-
      lapply(rasValue, function(x)
        setNames(data.frame(unlist(x)), unique(file_names)[i]))
    
    #add to tibble
    if (i == 1) {
      simDat <- tibble(id = c(1:nrow(mapFile[[i]]$coords)),
                       data = rasValue2)
    }
    else {
      simDat$data <- Map(cbind, simDat$data, rasValue2)
    }
    
  }
  names(simDat) <- c("grid_id", "clm")
  #Get spatial coordinate data from rasters for plotting
  simDat$x <- mapFile[[i]]$coords[, 1]
  simDat$y <- mapFile[[i]]$coords[, 2]
  saveRDS(simDat, paste0(save_file, "spatialChunk_", chunk, ".RDS"))
  
}



data_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\chessReg\\"
files <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE, 
                    recursive = T)
siteList<-basename(files)
siteList<-unique(gsub("\\-.*","",siteList))


regioCom<-function(siteName){
  data_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\chessReg\\"
  
  file_names<-list.files(data_dir,pattern = siteName,full.names = T)
  reg<-read.csv(file_names[1])
  reg$siteName<-siteName
  for(i in c(2:length(file_names))){
  regX<-data.frame(read.csv(file_names[i]))
  regX2<-as.data.frame(regX[,3])
  names(regX2)<-names(regX)[3]
  reg<-cbind(reg,regX2)
  }
  return(reg)
  
}











##spatial splitting function - split spatial data into chunks
#'@param data_dir directory which stores the CHESS spatial data (under a filename called CHESS) - also should include a folder for the full time series to be written into called fullTS
#'@param save_file where to save the output files
#'@param startDate the year to start the processing from, CHESS goes back to 1961 but you may not need that far back (significantly increases time to run function the further back you go)
#'@param variable whether to process all CHESS climate variables or select individual ones to process (e.g. precip, temp etc. )
#'@param numChunks function creates chunks of 10,000 grid cells going down the UK, 1:39 covers scotland, higher values cover the rest of UK
#'@return climate data with coordinates
spatSplitXScape<-function(data_dir,
                          save_file,
                          output_dir,
                          dates=c(1960:1990),
                          variable="all",
                          numChunks=c(1:10)) {
  files <- list.files(path = data_dir,
                      pattern = "\\.nc$",
                      full.names = TRUE,
                      recursive = T)
  file_names <- sub("\\.*", "",list.files(path = data_dir,
                                         pattern = "\\.nc$",
                                         full.names = F,
                                         recursive = T))
  
  file_names <- sub(".*/", "",list.files(path = data_dir,
                                        pattern = "\\.nc$",
                                        full.names = F,
                                        recursive = T))
  
  #only run merging once as it's pretty slow, hopefully once done this wont need to be done again
  #Merge monthly files into longer time-series
  
  file_names <- file_names[which(str_sub(file_names, -11, -8)%in% dates) ]
  
  variable_names <-
    if (variable == "all")
      unique(str_match(file_names, "01_\\s*(.*?)\\s*_uk")[, 2]) else
        variable
  print(variable_names)
  print("merging layers to create single files with full climate time-series - may take a while :)")
  for (i in 1:unique(length(unique(variable_names)))) {
    print(unique(variable_names)[i])
    ifelse(!dir.exists(file.path(output_dir, "/fullTS")), dir.create(file.path(output_dir, "fullTS")), FALSE)
    filesTmp <-
      paste0(data_dir, "/",variable_names[i],"/", file_names[grepl(variable_names[i], file_names) == T])
    
    if(variable=="tas"){
    filesTmp<-grep("tasmax", filesTmp, invert=TRUE, value = TRUE)
    filesTmp<-grep("tasmin", filesTmp, invert=TRUE, value = TRUE)
    }
    
    splitter <- function(j) {
      topRow <- ifelse(j == 1, 1, (j - 1) * 6)
      rastLayer <- lapply(filesTmp, function(x) {
        brick(x)
      })
      rastLayer <-
        lapply(rastLayer, function(x)
          raster::crop(x, extent(x, topRow, j * 6, 1, 656))) #kershope extent 475, 476, 352, 353
      
      rastLayer <- raster::brick(rastLayer)
      rastLayerX <- list(getValues(rastLayer))
      rastLayerX$coords <- coordinates(rastLayer)
      
      saveRDS(rastLayerX,
              paste0(
                output_dir,
                "/fullTS/",
                unique(variable_names)[i],
                "_",
                j,
                ".RDS"
              ))
    }
    
    coreNum <- detectCores()
    
    if (coreNum > 1) {
      plan(multisession, workers = coreNum - 1)
      future_map(numChunks, ~ splitter(j = .x), .progress = T)
    } else {
      map(numChunks, ~ splitter(j = .x), .progress = T
      )
      
    }
    
    
  }
  
  
}
#










#library(spatialrisk)
#
#regClm<-do.call(rbind,lapply(siteList,regioCom))
#regClm$Tmin<-regClm$tas-regClm$dtr
#regClm$Tmax<-regClm$tas+regClm$dtr
#bMark<-read.csv("data/bMarkCalibrationData.csv")
#bMark<-unique(bMark[,c("SiteIdentification","lat","lon")])
#regClm<-merge(regClm,bMark,by.x="siteName",by.y="SiteIdentification",all=F)
#regClm$date<-(as.Date(regClm$doy, origin = paste0(regClm$year-1,'-12-31')))
#regClm$month<-month(regClm$date)
#
#regClmSoils<-regClm[!duplicated(regClm[,c('siteName')]),]
#soilDat<-readRDS("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\soilDataLocs.RDS")
#
##find closest soil values
#ff<-purrr::map2_dfr(regClmSoils$lat, regClmSoils$lon, 
#                ~spatialrisk::points_in_circle(soilDat, .y, .x, 
#                                               lon = lon, 
#                                               lat = lat, 
#                                               radius = 1e6)[1,])
#
#ff$site<-regClmSoils$siteName
#names(regClm)<-c("siteName","year","doy","tempRange_c","specHumid","precip_mm","psurf_pa","solarRad_MJ","surface_wind","Tmean","Tmin","Tmax","lat","lon","date","month")
#regClm<-merge(regClm,ff[,c(6:11,13)],by.x="siteName",by.y="site")
#
##regClm<-regClm%>%group_by(siteName,year,month)%>%
##  summarise(precip=sum(precip),Tmean=mean(tas),Tmin=min(Tmin),Tmax=max(Tmax),lat=median(lat))
#
#saveRDS(regClm,"regionalClmDat.RDS")
#
#
#library(SPEI)
#hargFunc<-function(siteName){
#  regClmSite<-regClm[regClm$siteName==siteName,]
#  regClmSite$harRes<- hargreaves(regClmSite$Tmin, regClmSite$Tmax,  lat = regClmSite$lat[1])
#  regClmSite$BAL <- regClmSite$precip-regClmSite$harRes
#  regClmSiteTS <- ts(regClmSite[,-c(1:3)], end=c(2017,12), frequency=12)
#  spei12 <- spei(regClmSiteTS[,'BAL'], scale=6)
#  return(spei12)
#}
#
#ff<-lapply(siteList,hargFunc)
#
#par(mfrow=c(3,1))
#for(i in c(1:length(siteList))){
#  plot(ff[[i]],main=paste0("Site: ",siteList[i]))
#}
#
#ggplot(data=regClm)+
#  geom_line(aes(x=doy,y=precip,col=siteName))
#