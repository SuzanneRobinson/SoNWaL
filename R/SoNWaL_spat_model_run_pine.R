#' SoNWaL_spat_model_run_pine execute run of SoNWaL for spatial grid cell with single parameter set
#' @param params vector of parameter set to run model with
#' @param clm climate data to run model with
#' @return dataframe with model outputs, columns are output variables, rows are time-steps
#' @export
SoNWaL_spat_model_run_pine <- function(params,clm, scape=T) {
  #get default parameters
  baseParms <-getParmsPine(timeStp = 52, waterBalanceSubMods = T)

  #parameter names to update
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
        "SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")

  #Update parameters with proposals
  params<-as.data.frame(params)
  baseParms[nm] <- as.data.frame(params[nm])

  if(scape==F){
    clm$pyear<-1961

    clm$Tmax<-clm$tas+(clm$dtr/2)
    clm$Tmin<-clm$tas-(clm$dtr/2)
    clm$RH<-relative_humidity_calc(Tmean=clm$tas-273.15,pp=clm$psurf,spec_hum=clm$huss)

    clm<-clm%>%dplyr::select(RH,huss,precip,rlds,rsds,sfcWind,tas,Tmax,Tmin,pyear)

    names(clm) <-
      c(
        "RH",
        "specHumid",
        "Rain",
        "SolarRadLW",
        "SolarRad",
        "wind",
        "Tmean",
        "Tmax",        "Tmin",
        "pyear"
      )
  }

  if(scape==T){
    clm<-clm%>%
      dplyr::select("hurs",    "huss",    "pr",      "rlds",    "rsds",    "sfcWind", "tas",
                    "tasmax",  "tasmin", "pyear")

    names(clm) <-
      c(
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
    tibble::rownames_to_column(var = "Date") %>%
    mutate(Date = if(scape==F) as.Date(str_sub(Date,start = 2),"%Y.%m.%d") else as.Date(str_sub(Date,end = -10, start = 2),"%Y.%m.%d")) %>%
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

  # add climate data
  baseParms$weather <- as.data.frame(clm)

  #run model and calculate MAI and CAI for yield class
  out<-  do.call(SoNWaL, baseParms) %>%
    mutate(age = rev(as.numeric(max(Year, na.rm = T) - Year))) %>%
    mutate(MAI = (Vu / age), CAI = c(rep(0, 52), diff(Vu, lag = 52))) %>%
    mutate(yc=yield_class_calc(.))

  return(out)
}
