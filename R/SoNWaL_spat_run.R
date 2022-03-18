

# Spatial run function which takes a single grid cell and associated climate variables, 
# runs simulations for a range of MCMC posterior values
# Calculates the mean and a 95% credible interval
# Currently there is some data manipulation going on to spread a single years data
# over multiple years until all data is downloaded
#' @param site site data
#' @param clm climate data
#' @param param_draw parameter draws
SoNWaL_spat_run <-
  function(site,
           clm,
           param_draw,
           grid_id,
           soil_depth,
           wp,
           fc,
           sp,
           carbon,
           cond,
           N0,
           C0,
           plant_year) {
    
   clmY<- clm%>%
      tibble::rownames_to_column()
    
    clmY<-clmY[,1]
    clm$pyear<-stringr::str_sub(clmY,start=2, end =5)
    clm<-filter(clm,pyear>=plant_year)
    

    library(lubridate)

    # update parameters with fixed values from soil data etc.
      param_draw$pars <- lapply(param_draw$pars, function(x) {
        x$V_nr <- soil_depth/100
        x$maxRootDepth <- soil_depth/100
        x$startC <- carbon
        x$startN <- carbon/10
        x$wiltPoint <- wp
        x$fieldCap <- fc
        x$satPoint <- sp
        # currently K_s and drainage use same value from maps, could be improved?
        x$K_s <- cond
        x$K_drain <- cond
        x$K_drain_nrz <- cond
        #######################
        x$SWpower0<-N0
        x$SWconst0<-C0
        return(x)
      })
      
      
    # return empty cell output if there is no climate data (no land mass)
    if (is.na(clm[1, 2]) == T) {
      print("no climate data")
      return (tibble::as_tibble(
        data.frame(
          Year = NA,
          grid_id = grid_id,
          Wsbr_q05 = NA,
          Wsbr_q95 = NA,
          Wsbr_value = NA,
          Rs_q05 = NA,
          Rs_q95 = NA,
          Rs_value = NA,
          EvapTransp_q05 =
            NA,
          EvapTransp_q95 = NA,
          EvapTransp_value = NA,
          volSWC_rz_q05 =
            NA,
          volSWC_rz_q95 = NA,
          volSWC_rz_value = NA,
          yc_q05 = NA,
          yc_q95 = NA,
          yc_value = NA,
          GPP_q05 =
            NA,
          GPP_q95 = NA,
          GPP_value = NA,
          NPP_q05 =
            NA,
          NPP_q95 = NA,
          NPP_value = NA,
          NEE_q05 =
            NA,
          NEE_q95 = NA,
          NEE_value = NA,
          Reco_q05 =
            NA,
          Reco_q95 = NA,
          Reco_value = NA,
          LAI_q05 =
            NA,
          LAI_q95 = NA,
          LAI_value = NA
        )
      )) 
      } else {


      site_out <-
        tryCatch({

          param_draw %>%
            dplyr::mutate(sim = mapply(SoNWaL_spat_model_run, pars, MoreArgs = list(clm),SIMPLIFY = F)) %>%
            dplyr::select(mcmc_id, sim) %>%
            tidyr::unnest_legacy() %>%
            group_by(Year, mcmc_id) %>%
            dplyr::summarise(
              grid_id = grid_id,
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
              
            ) %>%
            group_by(Year) %>%
            dplyr::summarise(
              grid_id = grid_id,
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
        error = function(cond) {
          message(cond)
          
          site_out <-
            tibble::as_tibble(
              data.frame(
                Year = NA,
                grid_id = grid_id,
                Wsbr_q05 = NA,
                Wsbr_q95 = NA,
                Wsbr_value = NA,
                Rs_q05 =
                  NA,
                Rs_q95 = NA,
                Rs_value = NA,
                EvapTransp_q05 =
                  NA,
                EvapTransp_q95 = NA,
                EvapTransp_value = NA,
                volSWC_rz_q05 =
                  NA,
                volSWC_rz_q95 = NA,
                volSWC_rz_value = NA,
                yc_q05 =
                  NA,
                yc_q95 = NA,
                yc_value = NA,
                GPP_q05 =
                  NA,
                GPP_q95 = NA,
                GPP_value = NA,
                NPP_q05 =
                  NA,
                NPP_q95 = NA,
                NPP_value = NA,
                NEE_q05 =
                  NA,
                NEE_q95 = NA,
                NEE_value = NA,
                Reco_q05 =
                  NA,
                Reco_q95 = NA,
                Reco_value = NA,
                LAI_q05 =
                  NA,
                LAI_q95 = NA,
                LAI_value = NA
              )
            )
        })
    
    
    site_out <-
      site_out[!duplicated(site_out), ]
    return(as.data.frame(site_out))
    
  }

}




# execute run of SoNWaL for spatial grid cell with single parameter set
#' @param params vector of parameter set to run model with
#' @param clm climate data to run model with
#' @return dataframe with model outputs, columns are output variables, rows are time-steps
SoNWaL_spat_model_run <- function(params,clm) {
  #get default parameters
  baseParms <- getParms(timeStp = 52, waterBalanceSubMods = T)
  
  #parameter names to update
  nm<-c("wiltPoint","fieldCap","satPoint","K_s","V_nr","sigma_zR","shared_area","maxRootDepth","K_drain",
        "pFS2","pFS20","aS","nS","pRx","pRn","gammaFx","gammaF0","tgammaF","Rttover","mF","mR",
        "mS","SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
        "Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
        "SWpower0","Qa","Qb","MaxIntcptn","k","startN","startC")
  
  #Update parameters with proposals
  params<-as.data.frame(params)
  baseParms[nm] <- as.data.frame(params[nm])
  
  
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
  
  #convert to celsius from kelvin add date, weeks and months
  clm%>%
    mutate(Tmean = Tmean-273.15, Tmax = Tmax - 273.15, Tmin = Tmin - 273.15) %>% 
    tibble::rownames_to_column(var = "Date") %>%
    mutate(Date = as.Date(str_sub(Date,end = -10, start = 2),"%Y.%m.%d")) %>%
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
 out<-  do.call(fr3PGDN, baseParms) %>% 
    mutate(age = rev(as.numeric(max(Year, na.rm = T) - Year))) %>%
    mutate(MAI = (Vu / age), CAI = c(rep(0, 52), diff(Vu, lag = 52))) %>%
    mutate(yc=yield_class_calc(.))
  
  return(out)
}
