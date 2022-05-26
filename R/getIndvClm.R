#' @description function for getting climate data for individual sites extracted for NZ+ project
#' @export
getIndvClm<-function(scape=F, siteName = "Harwood", soil_dat, simDatLoc, baseParms){
  
  simDat<-readRDS(simDatLoc)%>%
    filter(site==siteName)
  
  
  clm<-simDat[,-1]
  
  if(scape==F){
    clm$pyear<-1961
    
    
    clm$Tmax<-clm$tas+(clm$dtr/2)
    clm$Tmin<-clm$tas-(clm$dtr/2)
    clm$RH<-relative_humidity_calc(Tmean=clm$tas-273.15,pp=clm$psurf,spec_hum=clm$huss)
    
    clm<-clm%>%dplyr::select(date, RH,psurf,precip,rlds,rsds,sfcWind,tas,Tmax,Tmin,pyear)
    
    names(clm) <-
      c("date",
        "RH",
        "specHumid",
        "Rain",
        "SolarRadLW",
        "SolarRad",
        "wind",
        "Tmean",
        "Tmax",
        "Tmin",
        "pyear"
      )
  }
  
  if(scape==T){
    clm$pyear<-2018
    
    clm<-clm%>%
      dplyr::select("date", "hurs",    "huss",    "pr",      "rlds",    "rsds",    "sfcWind", "tas",
                    "tasmax",  "tasmin", "pyear")
    
    names(clm) <-
      c("date",
        "RH",
        "specHumid",
        "Rain",
        "SolarRadLW",
        "SolarRad",
        "wind",
        "Tmean",
        "Tmax",
        "Tmin",
        "pyear"
      )
  }
  #convert to celsius from kelvin add date, weeks and months
  clm<-clm%>%
    mutate(Tmean = Tmean-273.15, Tmax = Tmax - 273.15, Tmin = Tmin - 273.15) %>% 
    mutate(Date = as.Date(date,"%Y-%m-%d"))%>%
    mutate(Month = month(Date), week = week(Date), Year = year(Date))
  
  
  #calc vpd
  clm$VPD <-
    ((((0.61078 * exp(17.269 * (clm$Tmean) / 
                        (237.3 + clm$Tmean))) * (1 - (clm$RH) / 100))))
  
  
  #get climate data into correct units
  clm <- clm %>%
    dplyr::group_by(Year, week) %>%
    dplyr::summarise(
      Month = median(Month),
      Tmax = max(Tmax),
      Tmin = min(Tmin),
      Tmean = mean(Tmean),
      Rain = sum(Rain * 86400),
      SolarRad = mean((SolarRad * 86400) / 1e+6),
      FrostDays = 0,
      MonthIrrig = 0,
      VPD = mean(VPD)
    )
  
  
  

  lat<-ifelse(siteName=="AH",51.1536,55.2159)
  lon<-ifelse(siteName=="AH",-0.8582,-2.0235)
  
  spatChunk<-data.frame(lat=lat, lon=lon)
  spatChunk$wp<-ext_eu_soil(soil_dat[1], spatChunk, val_name = "wp_map", 1000)[[1]]
  spatChunk$fc<-ext_eu_soil(soil_dat[2], spatChunk, val_name = "fc_map", 1000)[[1]]
  spatChunk$soil_cond<-10^(ext_eu_soil(soil_dat[3], spatChunk, val_name = "cond", 1000)[[1]])/10
  spatChunk$sp<-ext_eu_soil(soil_dat[4], spatChunk, val_name = "sp", 1000)[[1]]
  spatChunk<-cbind(spatChunk,ext_eu_soil(soil_dat[5], spatChunk, val_name = "sat", 1000, rast = F)[,-c(1:4)])
  simDat<-cbind(spatChunk,soil_regr(simDat= spatChunk,concChains=paramsFile)[,c(30,31)])
  
  
  
  
  baseParms$wiltPoint=simDat$wp
  baseParms$fieldCap=simDat$fc
  baseParms$satPoint=simDat$sp
  baseParms$K_s=simDat$soil_cond
  #baseParms$startC = simDat$SOC
  baseParms$V_nr=simDat$soil_depth/100
  baseParms$maxRootDepth=simDat$soil_depth/100
  baseParms$SWpower0=simDat$no
  baseParms$SWconst0=simDat$co
  #baseParms$startN<-baseParms$startC/10
  
  baseParms$weather<-as.data.frame(clm)
  
  return(baseParms)
}
