##UK spatial data extraction
library(tidyr)
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
library(dplyr)
library(lubridate)
library(parallel)
library(stringr)


##spatial splitting function - split spatial data into chunks
#'@param dataDir directory which stores the CHESS spatial data (under a filename called CHESS) - also should include a folder for the full time series to be written into called fullTS
#'@param saveFile where to save the output files
#'@param startDate the year to start the processing from, CHESS goes back to 1961 but you may not need that far back (significantly increases time to run function the further back you go)
#'@param variable whether to process all CHESS climate variables or select individual ones to process (e.g. precip, temp etc. )
#'@param numChunks function creates chunks of 10,000 grid cells going down the UK, 1:39 covers scotland, higher values cover the rest of UK
#'@return climate data with coordinates
spatSplitX<-function(dataDir,saveFile,outputDir,dates=c(1970:1980),variable="all",numChunks=c(1:10)){
  files <- list.files(path = dataDir, pattern = "\\.nc$", full.names = TRUE, 
                      recursive = T)
  fileNames <- sub("\\/.*", "",list.files(path = dataDir, pattern = "\\.nc$", full.names = F, 
                                          recursive = T))
  
  #only run merging once as it's pretty slow, hopefully once done this wont need to be done again
  #Merge monthly files into longer time-series

    fileNames <- fileNames[which(str_sub(fileNames, -11, -8)%in% dates) ]
    
    variableNames <-
      if (variable == "all")
        unique(str_match(fileNames, "met_\\s*(.*?)\\s*_gb")[, 2]) else
      variable
    print(variableNames)
    print("merging layers to create single files with full climate time-series - may take a while :)")
    for (i in 1:unique(length(unique(variableNames)))) {
      print(unique(variableNames)[i])
      ifelse(!dir.exists(file.path(outputDir, "/fullTS")), dir.create(file.path(outputDir, "fullTS")), FALSE)
      filesTmp <-
        paste0(dataDir, "/", fileNames[grepl(variableNames[i], fileNames) == T])
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
                  outputDir,
                  "/fullTS/",
                  unique(variableNames)[i],
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


##read in list of spatial climate data with grid coordinates from spatSplit function output and merge all climate vars into one table associated with each grid cell ##
#'@param chunk chunk number, 1:39 covers whole of scotland, can only go as high as files available from spatSplit
#'@param outputDir names of directory where files (outputs from spatSplit) to go merge are located
#'@param saveFile location to save merged files
#'@return tibble of site id key with associated dataframe of longitudinal climate data
spatDatUKnc <- function(chunk = 1, outputDir,saveFile) {
  library(stringr)
  library(dplyr)
  library(raster)
  
  files <-
    list.files(
      path = paste0(outputDir, "/fullTS"),
      pattern = paste0("_", chunk, "\\.RDS$"),
      full.names = TRUE,
      recursive = T
    )
  
  #read in files as rasters into list
  mapFile <- lapply(files, function(x) {
    (readRDS(x))
  })
  
  fileNames <- unique(str_match(files, "TS/\\s*(.*?)\\s*_")[, 2])
  
  #get climate values from raster
  for (i in c(1:length(files))) {
    print(i)
    rasValue = as.data.frame(mapFile[[i]][[1]])
    
    
    #Transpose data before putting into table
    rasValue <- rasValue %>% purrr::transpose()
    #Convert transposed data for each cell into a dataframe
    colNm <- fileNames[i]
    rasValue2 <-
      lapply(rasValue, function(x)
        setNames(data.frame(unlist(x)), unique(fileNames)[i]))
    
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
  
  simDatSp<-split(simDat, (seq(nrow(simDat))-1) %/% 45) 
  
  for(i in c(1:length(simDatSp))){
    saveRDS( simDatSp[[i]], paste0(saveFile, "spatialChunk_",chunk,"_", i+1, ".RDS")) 
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
#' @param param_draw parameter draws
FR3PG_spat_run <- function(site, clm,param_draw){
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

      #Update parameters with proposals 
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

      weatherDat$Rain<-weatherDat$Rain

      sitka$weather<-weatherDat
      
      #run model and output data
      out<-do.call(fr3PGDN,sitka)
      out$age<-rev(as.numeric(2018-out$Year))
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
#add na values where there is no data, in the sea etc. 
error = function(cond){
  site_out <- tibble::as_tibble(data.frame(Wsbr_q05=NA,Wsbr_q95=NA,Wsbr_value=NA,
                                           GPP_q05=NA,GPP_q95=NA,GPP_value=NA,
                                           NPP_q05=NA,NPP_q95=NA,NPP_value=NA,
                                           NEE_q05=NA,NEE_q95=NA,NEE_value=NA,
                                           Reco_q05=NA,Reco_q95=NA,Reco_value=NA))
})

    
    
  return(site_out)
  
}





##kershop forest data extraction
#
#for (i in 2:unique(length(unique(variableNames)))) {
# print(unique(variableNames)[i])
# ifelse(!dir.exists(file.path(outputDir, "/fullTS2")), dir.create(file.path(outputDir, "fullTS2")), FALSE)
# filesTmp <-
#     paste0(dataDir, "/", fileNames[grepl(variableNames[i], fileNames) == T])
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
#                outputDir,
#                "/fullTS2/",
#                unique(variableNames)[i],
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
##'@param outputDir names of directory where files (outputs from spatSplit) to go merge are located
##'@param saveFile location to save merged files
##'@return tibble of site id key with associated dataframe of longitudinal climate data
#spatDatUKnc2 <- function(chunk = 1, outputDir,saveFile) {
#  library(stringr)
#  library(dplyr)
#  library(raster)
#  
#  files <-
#    list.files(
#      path = paste0(outputDir, "/fullTS2"),
#      pattern = paste0("_", chunk, "\\.RDS$"),
#      full.names = TRUE,
#      recursive = T
#    )
#  
#  #read in files as rasters into list
#  mapFile <- lapply(files, function(x) {
#    (readRDS(x))
#  })
#  
#  fileNames <- unique(str_match(files, "TS2/\\s*(.*?)\\s*_")[, 2])
#  
#  #get climate values from raster
#  for (i in c(1:length(files))) {
#    rasValue = as.data.frame(mapFile[[i]][[1]])
#    
#    
#    #Transpose data before putting into table
#    rasValue <- rasValue %>% purrr::transpose()
#    #Convert transposed data for each cell into a dataframe
#    colNm <- fileNames[i]
#    rasValue2 <-
#      lapply(rasValue, function(x)
#        setNames(data.frame(unlist(x)), unique(fileNames)[i]))
#    
#    #add to tibble
#    if (i == 1) {
#      simDat <- tibble(id = c(1:nrow(mapFile[[i]]$coords)),
#                       data = rasValue2)
#    }
#    else {
#      simDat$data <- Map(cbind, simDat$data, rasValue2)
#    }
#    
#  }
#  names(simDat) <- c("grid_id", "clm")
#  #Get spatial coordinate data from rasters for plotting
#  simDat$x <- mapFile[[i]]$coords[, 1]
#  simDat$y <- mapFile[[i]]$coords[, 2]
#  saveRDS(simDat, paste0(saveFile, "spatialChunk_", chunk, ".RDS"))
#  
#}
#
#