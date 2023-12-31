#     ____     _ __  _      __     __            _____          __           __   #
#    / __/__  (_) / | | /| / /__ _/ /____ ____  / ___/__  ___  / /____ ___  / /_  #
#   _\ \/ _ \/ / /  | |/ |/ / _ `/ __/ -_) __/ / /__/ _ \/ _ \/ __/ -_) _ \/ __/  #
#  /___/\___/_/_/   |__/|__/\_,_/\__/\__/_/    \___/\___/_//_/\__/\__/_//_/\__/   #

#'soilWC soil water content calculator
#'@description Calculates soil water content for rooting zone
#'@param parms global 3PG parameter values
#'@param weather global weather values
#'@param state state values
#'@return SWC_rz SWC of rooting zone
#'@export

soilWC<-function(parms,weather,state,K_s,SWC,soilVol){
  
  #size of time-step. 
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
  
  #Volume of initial non-rooting zone
  V_nr = parms[["V_nr"]]
  
  #soil water in rooting zone at t0 (start of time-step)
  SWC_rz0 = state[["SWC_rz"]]

    #rooting depth / volume - Almedia describes this as depth, assumed to be proportional in paper to biomass
  z_r = max(min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["maxRootDepth"]]),0.05) # can't go deeper than non-rooting zone/max root depth

  V_rz = z_r #Almedia and Sands paper suggests volume of root zone is equivalent to z_r
  
  volSWC_rz<-min((SWC_rz0 /(z_r*1000)),parms[["satPoint"]])
  
  #Shared area, area is in m^2, so area around the tree?
  A = parms[["shared_area"]]
  
  #Non-root zone decreases as root zone increases, V_nr is max non-root zone volume
  V_nrx<-max(V_nr-V_rz,0)

  #Time constant - based on rate of movement soil water between zones
  t_s0 = (V_rz * V_nrx) / (K_s * A * (V_rz + V_nrx))
  
  #Current state of soil water content in non-rooting zone (at time-0 for month)
  SWC_nr0 = state[["SWC_nr"]]
  
  #State of soil water content in rooting zone at the end of the time step
  SWC_rz = (((SWC_rz0 * V_nrx - SWC_nr0 * V_rz) / (V_rz + V_nrx)) * exp(-t /t_s0)) +

    (V_rz / (V_rz + V_nrx) * (SWC_rz0 + SWC_nr0)) 
  # print(paste0("SWC_rz= ", SWC_rz ))
  return(SWC_rz)
}



#'drainageFunc calculate volume of water draining out of soil
#'@description Calculates drainage out of a soil zone
#'@param parms global 3PG parameter values
#'@param weather global weather values
#'@param SWC soil water content of zone to drain
#'@return drainage
#'@export
drainageFunc<-function(parms,weather,SWC,soilVol,K_drain){
  
  #size of time-step. This is for monthly time steps
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
 
  #Volumetric field capacity
  volSWC_fc= parms[["fieldCap"]]

  ##calculate drainage, convert soilVol to mm
  Qd<-(SWC-(volSWC_fc*soilVol*1000))*(1-exp(-K_drain*t))
  return(Qd)
}


#      ____     _ __  ____                             __  _          # 
#     / __/__  (_) / / __/  _____ ____  ___  _______ _/ /_(_)__  ___  # 
#    _\ \/ _ \/ / / / _/| |/ / _ `/ _ \/ _ \/ __/ _ `/ __/ / _ \/ _ \ # 
#   /___/\___/_/_/ /___/|___/\_,_/ .__/\___/_/  \_,_/\__/_/\___/_//_/ # 
#                               /_/                                   # 

#'
#'@description Calculates soil evaporation
#'@param parms global parameter values
#'@param weather global weather values
#'@param state state values
#'@param t_S1 duration of phase 1 evaporation
#'@param e_S1 potential rate of evaporation, equal to wet surface evaporation - could be driven by temp?
#'@param pseudo whether to use pseudo time-steps 
#'@return soil evaportation data
#'@export
soilEvap <-
  function(parms,
           weather,
           state,
           interRad,
           h,
           throughFall,
           pseudo = F) {
    e20 <- parms[["e20"]]
    rhoAir <- parms[["rhoAir"]]
    lambda <-
      parms[["lambda"]]#Volumetric latent heat of vaporization. Energy required per water volume vaporized (J/kg-1) (see penman monteith)
    #VPD significantly reduced at the soil level due to increase in below canopy humidity
    VPD <- weather[["VPD"]] * exp(state[["LAI"]] * (-log(4)) / 5)
    E_S1 = (parms[["E_S1"]])
    E_S2 = (parms[["E_S2"]])
    
    #size of time-step. This is for monthly time steps
    if (parms[["timeStp"]] == 12)
      t =   days_in_month(weather[["Month"]])
    if (parms[["timeStp"]] == 52)
      t =   7
    if (parms[["timeStp"]] == 365)
      t =  1
    if (is.numeric(t) == F)
      print("unsupported time step used")
    
    soilBoundaryCond <- 0.01
    soilCond <- 1e+10
    interRad <- max(interRad, 1e-5)
    Delta = 145 #Rate of change of saturation specific humidity with air temperature. (Pa/K−1)
    Cp = 1004 #specific heat cap of air J/kg-1
    Pa = 1.204 #Dry air density kg/m-3
    gamma = 66.1 # phsychrometric constant Pa/K-1
    rainDays <- max(weather[["rainDays"]], 1)
    rain <- throughFall / rainDays
    
    
    #whether to use pseudo time-steps for soil evaporation (more useful if using longer overall time-steps such as monthly)
    if (pseudo == T) {
      daysApart <- (365 / parms[["timeStp"]]) / rainDays
      eX = 0
      if (weather[["rainDays"]] > 0) {
        for (i in c(1:rainDays)) {
          E_S0 = max(E_S0 - rain, 0)
          t_S1 = E_S1 / e0
          
          if (E_S0 <= E_S1) {t0 = E_S0 / e0} else
            {
            t0 = as.numeric((((-2 * E_S0 * E_S1) + (E_S0 ^ 2) + (2 * E_S0) +
                                (E_S1 ^ 2) - (2 * E_S1) + 1 + (2 * e0 * E_S2 * t_S1) - (E_S2 ^ 2)) / (2 * e0 * E_S2)))}
          E_Sum = if ((daysApart + t0) <= t_S1) e0 * (daysApart + t0) - E_S0 else (E_S1 + E_S2 * (sqrt(1 + 2 * (e0 / E_S2) * ((daysApart + t0) - t_S1) - 1)))- E_S0
          
          E_S0 = E_S0 + E_Sum
          eX = eX + E_Sum}
      }
      return(list(E_S0, eX))
    }
    if (pseudo == F){
      e0 <-max(h * (soilCond * (Delta * interRad + soilBoundaryCond * Pa * Cp * ((VPD * 1000))) /(lambda * ((gamma + Delta) * soilCond + gamma * soilBoundaryCond))), 0)
      e0 <- e0 * 365 / parms[["timeStp"]]
      
      E_S0 = state[["E_S"]]
      
      E_S <- E_S0 + ifelse(E_S0 <= E_S1, e0, e0 / (1 + (E_S0 - E_S1) / E_S2))
      
      return(list(E_S, e0))
    }
  }



#' SWGmod soil water modifier (see landsberg and sands)
#' @param SWconst
#' @param SWpower
#' @param MoistRatio
#' @return soil water modifier
#' @export
SWGmod<-function(SWconst,SWpower,MoistRatio){
  f_theta<-(1-(1-MoistRatio)^SWpower)/(1+((1-MoistRatio)/SWconst)^SWpower)
#  f_theta<-(1-(1-MoistRatio)^SWpower)/(1+(1-2*SWconst^SWpower)*((1-MoistRatio)/SWconst)^SWpower)
  
  
  return(f_theta)
}

