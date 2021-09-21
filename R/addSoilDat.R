###UK spatial data extraction
#library(tidyr)
#library(purrr)
#library(BayesianTools)
#library(sensitivity)
#library(future)
#library(ncdf4) 
#library(raster) 
#library(ggplot2)
#library(httr)
#library(furrr)
#library(dplyr)
#library(lubridate)
#library(parallel)
#library(stringr)
#library(proj4)
#
#
#clm<-nc_open("C:/Users//aaron.morris//OneDrive - Forest Research//Documents//Projects//PRAFOR//models//spatial_met_data//CHESS//daily//chess-met_dtr_gb_1km_daily_19610101-19610131.nc")
##soilTex<-nc_open("C:/Users//aaron.morris//OneDrive - Forest Research//Documents//Projects//PRAFOR//models//spatial_soil_data//soilzLatLon2.nc")
#lon <- ncvar_get(clm, "lon")
#lat <- ncvar_get(clm, "lat", verbose = F)
#t <- ncvar_get(clm, "time")
#ndvi.array <- ncvar_get(clm, "dtr") # store the data in a 3-dimensional array
#dim(ndvi.array) 
#ndvi.slice <- ndvi.array[, , 1] 
#
#r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
#r <- flip(r, direction='y')
#clmRast <- as.data.frame(getValues(r))
#clmRast$y <- coordinates(r)[,2]
#clmRast$x <- coordinates(r)[,1]
#names(clmRast)<-c("vals","lat","lon")
#clmRast$x<-ncvar_get(clm, "x")
#clmRast$y<-ncvar_get(clm, "y")
#clmRast<-clmRast[,-1]
#
#
###Read in texture data
#soilDatExt<-function(fileLoc,varName="surf_text_dom"){
#soilTex<- brick(fileLoc,
#                crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
#latLon<-coordinates(soilTex)
#proj4string <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
#+ellps=GRS80 +units=m +no_defs"
## Source data
#xy <- data.frame(x=latLon[,1], y=latLon[,2])
## Transformed data
#pj <- proj4::project(xy, proj4string, inverse=TRUE)
#soilDB <- data.frame(lat=pj$y, lon=pj$x)
#soilVals<-getValues(soilTex)[,1]
#soilDB$vals<-soilVals
#names(soilDB)<-c("lat","lon",varName)
#return(soilDB)
#}
#
#soilTexDom<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\txsrfse_directory_dom_stu\\txsrfse\\w001001.adf"
#soilDB<-soilDatExt(soilTexDom,varName="soilTexDom")
#
#
#
#clmRast$lat<-round(clmRast$lat,2)
#clmRast$lon<-round(clmRast$lon,2)
#soilDB$lat<-round(soilDB$lat,2)
#soilDB$lon<-round(soilDB$lon,2)
#
#soilDB2<-left_join(clmRast, soilDB, by = c("lat" = "lat", "lon" = "lon"))
#
#
#
#
#soilTexDom<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\txsrfse_directory_dom_stu\\txsrfse\\w001001.adf"
#soilDB<-soilDatExt(soilTexDom,varName="soilTexDom")
#
#soilTexDom<-"C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\txsrfse_directory_dom_stu\\txsrfse\\w001001.adf"
#soilDB<-soilDatExt(soilTexDom,varName="soilTexDom")

