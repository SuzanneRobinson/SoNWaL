
#' spat_dat_scape
#' @description read in list of CHESS-scape spatial
#' climate data with grid coordinates from spatSplit
#' function output and merge all climate vars into one
#' table associated with each grid cell
#'@param chunk chunk number, 1:39 covers whole of
#'scotland, can only go as high as files available from spatSplit
#'@param output_dir names of directory where files
#'(outputs from spatSplit) to go merge are located
#'@param save_file location to save merged files
#'@return tibble of site id key with associated dataframe
#'of longitudinal climate data
#'@export
spat_dat_scape <- function(chunk = 1, output_dir, save_file) {
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

  # read in files as rasters into list
  map_file <- lapply(files, function(x) {
    (readRDS(x))
  })

  # get file names
  file_names <- unique(str_match(files, "TS/\\s*(.*?)\\s*_")[, 2])
  # get names for ID vals
  grid_ID_names <- unique(str_match(files, "TS/\\s*(.*?)\\s*.RDS")[, 2])


  # get climate values from rasters
  for (i in c(seq_along(files))) {
    print(i)
    ras_value <- as.data.frame(map_file[[i]][[1]])
    # Transpose data before putting into table
    ras_value <- ras_value %>% purrr::transpose()
    # Convert transposed data for each cell into a dataframe
    ras_value_tmp <-
      lapply(ras_value, function(x)
        setNames(data.frame(unlist(x)), unique(file_names)[i]))
    
    # add to tibble and give unique grid cell ID (grid ID poss redundant)
    if (i == 1) {
      sim_dat <- tibble(id = c(paste0(grid_ID_names[i],
                                      "_", seq_along(map_file[[i]]$coords[, 1]))),
                        data = ras_value_tmp)
    }  else {
      sim_dat$data <- Map(cbind, sim_dat$data, ras_value_tmp)
    }

}
  names(sim_dat) <- c("grid_id", "clm")
  # Get spatial coordinate data from rasters
sim_dat$x <- map_file[[i]]$coords[, 1]
sim_dat$y <- map_file[[i]]$coords[, 2]

# split into multiple files to avoid ultra large single files
# (could split by coordinates...probs still too large splitting by lat)
  sim_dat_sp <- split(sim_dat,
                      (seq(nrow(sim_dat)) - 1) %/% 45)

  for (i in c(seq_along(sim_dat_sp))) {
    saveRDS(sim_dat_sp[[i]],
            paste0(save_file,
                   "spatialChunk_", chunk, "_", i + 1, ".RDS"))
  }
return(print("manip complete"))
}



#'spat_split_scape 
#'@description spatial splitting function - split spatial
#' CHESS-scape data into chunks
#'@param data_dir directory which stores the CHESS spatial
#'data (under a filename called CHESS) - also should include
#'a folder for the full time series to be written into called fullTS
#'@param save_file where to save the output files
#'@param startDate the year to start the processing from, CHESS goes
#'back to 1961 but you may not need that far back
#'(significantly increases time to run function the further back you go)
#'@param variable whether to process all CHESS climate variables or select
#'individual ones to process (e.g. precip, temp etc. )
#'@param num_chunks function creates chunks of 10,000
#'grid cells going down the UK, 1:39 covers scotland,
#'higher values cover the rest of UK
#'@return climate data with coordinates
#'@export
spat_split_scape <- function(data_dir,
                           save_file,
                           output_dir,
                           dates=c(1960:1990),
                           variable="all",
                           num_chunks=c(1:10)) {
  files <- list.files(path = data_dir,
                      pattern = "\\.nc$",
                      full.names = TRUE,
                      recursive = T)
  file_names <- sub("\\.*", "", list.files(path = data_dir,
                                          pattern = "\\.nc$",
                                          full.names = F,
                                          recursive = T))

  file_names <- sub(".*/", "", list.files(path = data_dir,
                                         pattern = "\\.nc$",
                                         full.names = F,
                                         recursive = T))

  # only run merging once as it's pretty slow,
  # hopefully once done this wont need to be done again
  # Merge monthly files into longer time-series

  file_names <- file_names[which(str_sub(file_names, -11, -8) %in% dates)]
  variable_names <-
    if (variable == "all")
      unique(str_match(file_names, "01_\\s*(.*?)\\s*_uk")[, 2]) else
        variable
  print(variable_names)
  print("merging layers to create single files with full
        climate time-series - may take a while :)")
  for (i in 1:unique(length(unique(variable_names)))) {
    print(unique(variable_names)[i])
    ifelse(!dir.exists(file.path(output_dir, "/fullTS")),
           dir.create(file.path(output_dir, "fullTS")), FALSE)
    files_tmp <-
      paste0(data_dir, "/", variable_names[i], "/",
             file_names[grepl(variable_names[i], file_names) == T])

    if (variable == "tas") {
      files_tmp <- grep("tasmax", files_tmp, invert = TRUE, value = TRUE)
      files_tmp <- grep("tasmin", files_tmp, invert = TRUE, value = TRUE)
    }

    splitter <- function(j) {
      top_row <- ifelse(j == 1, 1, (j - 1) * 6)
      rast_layer <- lapply(files_tmp, function(x) {
        brick(x)
      })

      rast_layer <-
        lapply(rast_layer, function(x)
          raster::crop(x, extent(x, top_row,
                                 j * 6, 1, 656))) #kershope extent 475, 476, 352, 353

rast_layer <- raster::brick(rast_layer)
rast_layer_tmp <- list(getValues(rast_layer))
rast_layer_tmp$coords <- coordinates(rast_layer)

      saveRDS(rast_layer_tmp,
              paste0(
                output_dir,
                "/fullTS/",
                unique(variable_names)[i],
                "_",
                j,
                ".RDS"
              ))
}

core_num <- parallel::detectCores()

    if (core_num > 1) {
      plan(multisession, workers = core_num - 1)
      furrr::future_map(num_chunks, ~ splitter(j = .x), .progress = T)
    } else {
      map(num_chunks, ~ splitter(j = .x), .progress = T
)

    }
  }
}

#########################################################add soil data##################################################################


#' chess_xy_latlon 
#' @description convert xy vals from CHESS to lat lon
#' @param x CHESS x coord
#' @param y CHESS y coord
#' @return dataframe of lat lon
#' @export
chess_xy_latlon<-function(x,y){
  proj4string <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"
  xy<-data.frame(x=x,y=y)
  pj <- proj4::project(xy, proj4string, inverse=TRUE)
  return(data.frame(lat=pj$y, lon=pj$x))
}


#' soil_regr
#' @description use regression models to predict some unkown soil texture parameters from known parameters
#' @param simDat current known parameters/variables
#' @param mcmcReg output from MCMC for fitting regional models
#' @return simDat but with predicted soil texture params, es1,es2,swconst,swpower
#' @export
soil_regr<-function(simDat,concChains){
# data from landsberg & sands book
  texT<-data.frame(sand=c(30 ,33,42, 82, 92, 52, 60, 65, 7, 7, 10, 20),
                 clay =c(50, 34 ,18, 6, 5, 42, 28 ,10, 6, 47, 34, 20),
                 silt=c(20, 33 ,40, 12, 3, 6, 12, 25, 87, 46, 56, 60),
                 no=c(3,5,6,8,9,4,5.5,7,7.5,3.5,4.5,6.5),
                 co=c(0.4,0.5,0.55,0.65,0.7,0.45,0.525,0.6,0.625,0.425,0.475,0.575)
)


# models for predicting approximate drainage values 
 mod_K_drain<-lm(K_drain~wiltPoint+fieldCap+V_nr,data=concChains)
 mod_K_s<-lm(K_s~satPoint+wiltPoint+fieldCap+V_nr,data=concChains)

# models for predicting constant and power fertility modifiers from soil texture (driven by clay content)
modCo<-lm(co~clay,data=texT)
modNo<-lm(no~clay,data=texT)

# predict drainage values from mapped soil data
 #simDat$K_drain<-as.vector(predict(mod_K_drain,data.frame(wiltPoint=simDat$wp/(simDat$depth*10),fieldCap=simDat$fc/(simDat$depth*10),V_nr=simDat$depth/100)))
 #simDat$K_s<-as.vector(predict(mod_K_s,data.frame(satPoint=simDat$sp,wiltPoint=simDat$wp/(simDat$soil_depth*10),fieldCap=simDat$fc/(simDat$soil_depth*10),V_nr=simDat$soil_depth/100)))

# predict fertility modifiers from mapped soil data
simDat$co <- as.vector(predict(modCo, data.frame(sand = simDat$Tsand, clay=simDat$Tclay, silt=simDat$Tsilt)))
simDat$no <- as.vector(predict(modNo, data.frame(sand = simDat$Tsand, clay=simDat$Tclay, silt=simDat$Tsilt)))

# derive wp fc etc from values in mm and soil depth in cm
simDat$wp<-simDat$wp/(simDat$soil_depth*10)
simDat$fc<-simDat$fc/(simDat$soil_depth*10)
simDat$wp[simDat$wp>=1]<-1
simDat$fc[simDat$fc>=1]<-1

simDat$wp[is.na(simDat$wp)==T]<-0
simDat$fc[is.na(simDat$fc)==T]<-0
simDat$sp[is.na(simDat$sp)==T]<-0
simDat$soil_cond[is.na(simDat$soil_cond)==T]<-0
simDat$soil_depth[is.na(simDat$soil_depth)==T]<-0
return(simDat)
}



#' add_BGS_dat 
#' @description add soil data from Astley (BGS data)
#' @param spat_chunk spatial chunk to add soil data to
#' @param soil_dat location of BGS soil data files
#' @return updated spat_chunk dataframe with additional soil data
#' @export
add_BGS_dat<-function(spat_chunk, soil_dat){
  
  #identify if climate data is  present for each grid cell
  spat_chunk$clm_pres<-!is.na(lapply(spat_chunk$clm,sum))
  
  spat_chunk<-spat_chunk%>%
    tibble::add_column(chess_xy_latlon(spat_chunk$x,spat_chunk$y))
  
  #' ext_eu_soil
  #' @description function to extract soil data from EU maps and match with closest grid ref of climate
  #' @param map_location file location for eu map (wp,sp,cond etc.)
  #' @param spat_chunk spatial chunk file containing grid square locations, climate and metadata
  #' @param val_name name to use for extracted data
  #' @param max_distance if closest grid square in soil data to climate data is empty (e.g. water, no data) 
  #' how far to look for next grid square before defaulting to NA
  #' @param rast is the input file one of the eu raster files or astleys bgs data
  ext_eu_soil<-function(data_location, spat_chunk, val_name, max_distance, rast =T){
    
    if(rast ==T){
    eu_soil_dat<-raster(data_location)
    eu_soil_dat<-raster::crop(eu_soil_dat, extent(eu_soil_dat, 1000,
                                                  2500, 3200, 3900))
    sr<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    eu_soil_dat<-projectRaster(eu_soil_dat, crs = sr)
    eu_soil_dat<-data.frame(eu_soil_dat=getValues(eu_soil_dat),coordinates(eu_soil_dat))
     
    } else 
    {
      eu_soil_dat<- readRDS(data_location)
    }
    
    eu_soil_dat<-eu_soil_dat%>%
      rename(any_of(c("lon" = "x", "lat" = "y", "lon" = "long", "fc" = "FC", "wp" = "WP", "soil_depth" = "depth")))
    
    closest_vals<-purrr::map2_dfr(spat_chunk$lat, spat_chunk$lon, 
                                  ~na.omit(spatialrisk::points_in_circle(eu_soil_dat, .y, .x, 
                                                                         lon = lon, 
                                                                         lat = lat, 
                                                                         radius = 1e5))[1,])

    closest_vals[is.na(closest_vals$distance_m) ==T,]<-0
    if(rast ==T) closest_vals[,1][closest_vals$distance_m>max_distance]<-NA
    if(rast ==T) names(closest_vals)[1]<-val_name
    if(rast ==F) closest_vals[closest_vals$distance_m > max_distance, !names(closest_vals) %in% c("lon","lat")]<-NA
    closest_vals[closest_vals$lat==0,]<-NA
    
    return(closest_vals)
  }
  
  # extract eu map soil data
  closest_vals_wp<-ext_eu_soil(soil_dat[1],
                                 spat_chunk,"wp_map", 1000)
  
  closest_vals_fc<-ext_eu_soil(soil_dat[2],
                            spat_chunk,"fc_map", 1000)
  print("soilDat_read in")
  closest_vals_cond<-ext_eu_soil(soil_dat[3],
                                 spat_chunk,"soil_cond", 1000) %>%
    mutate(soil_cond=10^soil_cond/10)#convert from log10 cm/day
  print("sat point read in")
  
  closest_vals_sat<-ext_eu_soil(soil_dat[4],
                                 spat_chunk,"sp", 1000)
  
  closest_vals<-ext_eu_soil(soil_dat[5],
                            spat_chunk,"NA", 1000, F)
  print("cond read in")
  
# use closest value known if center of grid square is in the sea or data is missing, over 1km however may be too far?
  
  spat_chunk<-spat_chunk%>%
    tibble::add_column(dplyr::select(closest_vals,-lat,-lon))%>%
  tibble::add_column(dplyr::select(closest_vals_cond,-lat,-lon,-distance_m))%>%
  tibble::add_column(dplyr::select(closest_vals_sat,-lat,-lon,-distance_m)) %>%
    tibble::add_column(dplyr::select(closest_vals_wp,-lat,-lon,-distance_m)) %>%
  tibble::add_column(dplyr::select(closest_vals_fc,-lat,-lon,-distance_m))
  
  

  
  return(spat_chunk)
}













