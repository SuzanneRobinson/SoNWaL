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
#clm<-nc_open("C:/Users//aaron.morris//OneDrive - Forest Research//Documents//Projects//PRAFOR//models//spatial_met_data//CHESS//daily//chess-met_huss_gb_1km_daily_19610101-19610131.nc")
##soilTex<-nc_open("C:/Users//aaron.morris//OneDrive - Forest Research//Documents//Projects//PRAFOR//models//spatial_soil_data//soilzLatLon2.nc")
#lon <- ncvar_get(clm, "lon")
#lat <- ncvar_get(clm, "lat", verbose = F)
#t <- ncvar_get(clm, "time")
#ndvi.array <- ncvar_get(clm, "huss") # store the data in a 3-dimensional array
#dim(ndvi.array) 
#ndvi.slice <- ndvi.array[, , 1] 
#
#r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
#r <- flip(r, direction='y')
#plot(r)
#clmRast <- as.data.frame(getValues(r))
#clmRast$y <- coordinates(r)[,2]
#clmRast$x <- coordinates(r)[,1]
#names(clmRast)<-c("vals","lat","lon")
#clm<-brick("C:/Users//aaron.morris//OneDrive - Forest Research//Documents//Projects//PRAFOR//models//spatial_met_data//CHESS//daily//chess-met_huss_gb_1km_daily_19610101-19610131.nc")
#
#clmRast$x<-coordinates(clm)[,1]
#clmRast$y<-coordinates(clm)[,2]
##clmRast<-clmRast[,-1]
#
#
#soilTex<- brick("C:\\Users\\aaron.morris\\OneDrive - Forest Research\\Documents\\Projects\\PRAFOR\\models\\spatial_soil_data\\txsrfdo_directory_dom_stu\\txsrfdo\\w001001.adf",
#                crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
#
#soilTex <- projectRaster(soilTex, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0",method = "ngb") 
#sg<-crop(soilTex,r)
#sg2 <- projectRaster(sg, r,method = "ngb")
#mergedSoil<-stack(r,sg2)
#clmRast$soil<-getValues(mergedSoil$txsrfdo)
#
#
#ggplot() +
#  geom_raster(data = clmRast , aes(x = x, y = y, fill = as.factor(soil)))+
#  theme_bw()+
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank())+
#  ggtitle(i)+
#  coord_equal() + 
#  scale_fill_viridis_d("Chunk",option = "D")
