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
#'@param dataDir directory which stores the CHESS spatial data (under a filename called CHESS) - also should include a folder for the full time series to be written into called fullTS
#'@param saveFile where to save the output files
#'@param startDate the year to start the processing from, CHESS goes back to 1961 but you may not need that far back (significantly increases time to run function the further back you go)
#'@param variable whether to process all CHESS climate variables or select individual ones to process (e.g. precip, temp etc. )
#'@param numChunks function creates chunks of 10,000 grid cells going down the UK, 1:39 covers scotland, higher values cover the rest of UK
#'@return climate data with coordinates
spatSplitX<-function(dataDir,saveFile,outputDir,dates=c(1960:1980),variable="all",numChunks=c(1:10)){
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
      path = paste0(outputDir, "fullTS"),
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
  
 # plot(CAI,col="white")
 # lines(MAI,col="blue")
  #lines(CAI,col="red")
  
  
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
FR3PG_spat_run <- function(site, clm,param_draw,grid_id,soil){
  library(lubridate)
  if(soil==2) {param_draw$pars <- lapply(param_draw$pars, function(x) {x$wiltPoint<-0.126
  x$fieldCap<-0.268
  x$satPoint<-0.461
 return(x)})}
  
  if(soil==1) {param_draw$pars <- lapply(param_draw$pars, function(x) {x$wiltPoint<-0.057
  x$fieldCap<-0.122
  x$satPoint<-0.46
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
      nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","E_S1","E_S2","shared_area","maxRootDepth","K_drain",
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
      out$age<-rev(as.numeric(2018-out$Year))
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