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


##spatial splitting function - split historical CHESS spatial data into chunks
#'@param data_dir directory which stores the CHESS spatial data (under a filename called CHESS) - also should include a folder for the full time series to be written into called fullTS
#'@param save_file where to save the output files
#'@param startDate the year to start the processing from, CHESS goes back to 1961 but you may not need that far back (significantly increases time to run function the further back you go)
#'@param variable whether to process all CHESS climate variables or select individual ones to process (e.g. precip, temp etc. )
#'@param numChunks function creates chunks of 10,000 grid cells going down the UK, 1:39 covers scotland, higher values cover the rest of UK
#'@return climate data with coordinates
#'@export
spat_split_hist <-function(data_dir, 
                        save_file, 
                        output_dir, 
                        dates=c(1960:1990),
                        variable="all", 
                        num_Chunks=c(1:10)) {
  files <- list.files(path = data_dir, 
                      pattern = "\\.nc$", 
                      full.names = TRUE,
                      recursive = T)
  file_names <- sub("\\/.*", "",list.files(path = data_dir, 
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
        future_map(num_Chunks, ~ splitter(j = .x), 
                   .progress = T)
      } else {
        map(num_Chunks, ~ splitter(j = .x), 
            .progress = T)
    }
    
    
  }
  

}


##read in list of historical CHESS spatial climate data with grid coordinates from spatSplit function
#output and merge all climate vars into one table associated with each grid cell ##
#'@param chunk chunk number, 1:39 covers whole of scotland, can only go as high as files available from spatSplit
#'@param output_dir names of directory where files (outputs from spatSplit) to go merge are located
#'@param save_file location to save merged files
#'@return tibble of site id key with associated dataframe of longitudinal climate data
#'@export
spat_dat_hist <- function(chunk = 1, output_dir,save_file) {
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



#library(spatialrisk)

#regioCom<-function(siteName){
#  data_dir="C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_met_data\\chessReg\\"
#  
#  file_names<-list.files(data_dir,pattern = siteName,full.names = T)
#  reg<-read.csv(file_names[1])
#  reg$siteName<-siteName
#  for(i in c(2:length(file_names))){
#    regX<-data.frame(read.csv(file_names[i]))
#    regX2<-as.data.frame(regX[,3])
#    names(regX2)<-names(regX)[3]
#    reg<-cbind(reg,regX2)
#  }
#  return(reg)
#  
#}

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

 data.frame(lon = c(-3.488368,-238439, -329229),
                        lat = c(53.0561,54.01876, 54.67536))

 site_locs <- data.frame(lon=
c(54.67536,
54.018758,
51.097886,
52.221398,
51.962955,
52.852497,
53.056101,
52.846386,
51.675991,
51.580907,
51.759179,
51.945555,
51.903267,
55.380781,
58.099751,
56.66524),
lat=
c(-3.2922868,
  -2.3843855,
  -3.112056,
  -3.7814199,
  -3.6890579,
  -3.6990238,
  -3.4883681,
  -3.8829362,
  -3.5916353,
  -3.7223068,
  -3.6936558,
  -4.2178236,
  -3.8602517,
  -5.6133599,
  -4.4234889,
  -3.4785394))
