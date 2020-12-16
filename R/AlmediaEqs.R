########################################################################################         
#     _    _                    _ _              _           _     _____      _        #
#    / \  | |_ __ ___   ___  __| (_) __ _    ___| |_    __ _| |   | ____|__ _( )___    #
#   / _ \ | | '_ ` _ \ / _ \/ _` | |/ _` |  / _ \ __|  / _` | |   |  _| / _` |// __|   #
#  / ___ \| | | | | | |  __/ (_| | | (_| | |  __/ |_  | (_| | |_  | |__| (_| | \__ \_  #
# /_/   \_\_|_| |_| |_|\___|\__,_|_|\__,_|  \___|\__|  \__,_|_(_) |_____\__, | |___(_) #
#                                                                          |_|         #
########################################################################################          



#     ____     _ __  _      __     __            _____          __           __   #
#    / __/__  (_) / | | /| / /__ _/ /____ ____  / ___/__  ___  / /____ ___  / /_  #
#   _\ \/ _ \/ / /  | |/ |/ / _ `/ __/ -_) __/ / /__/ _ \/ _ \/ __/ -_) _ \/ __/  #
#  /___/\___/_/_/   |__/|__/\_,_/\__/\__/_/    \___/\___/_//_/\__/\__/_//_/\__/   #

#'@description Calculates soil water content for rooting zone
#'@param parms global 3PG parameter values
#'@param weather global weather values
#'@param state state values

#Internal function variable descriptions
#'@param t length of time-step
#'@param V_nr volume of non-rooting zone
#'@param V_rz volume of rooting zone
#'@param A shared area
#'@param z_r depth (or volume?) of rooting zone
#'@param sigma_zR area/depth explored by 1kg of root biomass
#'@param sigma_nr0 SWC of non-rooting zone at time 0
#'@param sigma_rz0 SWC of rooting zone at time 0

#'@return sigma_rz SWC of rooting zone

soilWC<-function(parms,weather,state){
  
  #size of time-step. 
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
  
  #Volume of initial non-rooting zone
  V_nr = parms[["V_nr"]]
  
  #soil water in rooting zone at t0 (start of time-step)
  sigma_rz0 = state[["ASW"]]
  
  #Soil conductivity - see Landsberg book for more details on this 
  K_s = parms[["K_s"]]
  
  #rooting depth / volume - Almedia describes this as depth, but sigma_zR could also be used to derive volume from root biomass?
  z_r = min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),V_nr) # can't go deeper than non-rooting zone
  
  V_rz = z_r # not sure if this is equivalent, or whether there needs to be a conversion from depth to volume?
  #Shared area, area is in m^2, so area around the tree?
  A = 5
  
  #Non-root zone decreases as root zone increases, V_nr is max non-root zone volume
  V_nrx<-max(V_nr-V_rz,0)
  #Time constant - based on rate of movement soil water between zones
  t_s0 = (V_rz * V_nrx) / (K_s * A * (V_rz + V_nrx))
  
  #Current state of soil water content in non-rooting zone (at time-0 for month)
  sigma_nr0 = state[["sigma_nr0"]]# should be dynamic? but not sure how to implement just yet...
  
  
  #State of soil water content in rooting zone at the end of the time step
  sigma_rz = (((sigma_rz0 * V_nrx - sigma_nr0 * V_rz) / (V_rz + V_nrx)) * exp(-t /t_s0)) +
    (V_rz / (V_rz + V_nrx) * (sigma_rz0 + sigma_nr0)) 
 
  return(sigma_rz)
}



#'@description Calculates soil water content for non-rooting zones
#'@param parms global 3PG parameter values
#'@param weather global weather values
#'@param state state values

#Internal function variable descriptions
#'@param t length of time-step
#'@param V_nr volume of non-rooting zone
#'@param V_rz volume of rooting zone
#'@param A shared area
#'@param z_r depth (or volume?) of rooting zone
#'@param sigma_zR area/depth explored by 1kg of root biomass
#'@param sigma_nr0 SWC of non-rooting zone at time 0
#'@param sigma_rz0 SWC of rooting zone at time 0

#'@return soilWC_NRZ SWC of non-rooting zone

soilWC_NRZ<-function(parms,weather,state){
  #size of time-step. This is for monthly time steps
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
  
  #Volume of initial non-rooting zone
  V_nr = parms[["V_nr"]]
  #soil water in rooting zone at t0, equiv ASW? depends if rain is assumed to occur at begining or end of month,
  #need to be careful so as not to mess up biomass allocation function which comes in after in the same time step (see RunModel.r)
  sigma_rz0 = state[["ASW"]]
  #Soil conductivity - see Landsberg book for more details on this
  K_s = parms[["K_s"]]
  #rooting depth / volume - Almedia describes this as depth, but sigma_zR could also be used to derive volume from root biomass?
  z_r = min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),V_nr) # can't go deeper than non-rooting zone
  V_rz = z_r # not sure if this is equivalent, or whether there needs to be a conversion from depth to volume?
  #Shared area
  A = 5
  #Non-root zone decreases as root zone increases, V_nr is max non-root zone volume
  V_nrx<-max(V_nr-V_rz,0)
  
  ##calculate max SWC of non-rooting zone so SWC can not go over this
  theta_fc_nr= parms[["theta_fc"]]*V_nrx*1000

  
  # Re-arrange equation A.14 to get value for non-rooting zone
  t_s0 = (V_rz * V_nrx) / (K_s * A * (V_rz + V_nrx))
  sigma_nr0 =max(((state[["ASW"]]*V_rz)+(state[["ASW"]]*V_nrx)-(exp(-t/t_s0)*sigma_rz0*V_nrx)-(sigma_rz0*V_rz))/((1-exp(-t/t_s0))*V_rz),0)

  #SWC of non-root zone can't be above field capacity of non-root zone
  sigma_nr0=min(sigma_nr0,theta_fc_nr)
  
  
  return(sigma_nr0)
}

#      ____     _ __  ____                             __  _          # 
#     / __/__  (_) / / __/  _____ ____  ___  _______ _/ /_(_)__  ___  # 
#    _\ \/ _ \/ / / / _/| |/ / _ `/ _ \/ _ \/ __/ _ `/ __/ / _ \/ _ \ # 
#   /___/\___/_/_/ /___/|___/\_,_/ .__/\___/_/  \_,_/\__/_/\___/_//_/ # 
#                               /_/                                   # 

#'@description Calculates soil evaporation
#'@param parms global parameter values
#'@param weather global weather values
#'@param state state values
#'@param t_S1 duration of phase 1 evaporation
#'@param e_S1 potential rate of evaporation, equal to wet surface evaporation - could be driven by temp?
#- calc from penman montieth - also see landsdown book and discussion on methods by Choudhury, Zhang etc.
#'@param E_S1 threshold - based on soil structure
#'@param E_S2 how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
#'@param E_S0 initial condition based on E_S
#'@return E_S Accumulated soil evaporation during the time step 
#Calc using penman monteith eq. g_c is infinite - could be better modified by soil dryness etc.


soilEvap<-function(parms,weather,state,interRad){
  
  e20 <- parms[["e20"]]
  rhoAir <- parms[["rhoAir"]]
  lambda <- parms[["lambda"]]
  VPDconv <- parms[["VPDconv"]]
  soilCond <- 0.01 #lower than canopy conductivity
  VPD <- weather[["VPD"]]
  E_S1 = parms[["E_S1"]]
  E_S2 = parms[["E_S2"]]



  
  #size of time-step. This is for monthly time steps
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
  
  

  
  #Potential wet surface evaporation rate - using Penman Monteith eq. mm per day loss
  e0 <- max(((e20 * interRad + rhoAir * lambda * VPDconv * VPD *
                  soilCond) / (1 + e20 + soilCond / Inf)), 1e-6) # 1e-6 to avoid division by zero...not sure it'd ever go below this though

  #E_S0 is E_S at the start of the time-step
  E_S0 = max(state[["E_S"]] , 0)
  

  #Duration of phase 1 evaporation
  t_S1 = E_S1 / e0
  
  #Solved for t equation A.10 in Almedia, to get equivalent t for E_S0 value
  t0 = as.numeric(round(((-2 * E_S0 * E_S1) + (E_S0 ^ 2) + (2 * E_S0) +
                           (E_S1 ^ 2) - (2 * E_S1) + 1 + (2 * e0 * E_S2 * t_S1) - (E_S2 ^ 2)
  ) / (2 * e0 * E_S2)))
  
  
  #Integrate equation A.9 to get value at time t (assuming t is number of days in month)
  #and Calc E_S using t+t0 to get amount of evaporation between t0 and t

  E_S = if (t <= t_S1)
    e0 * t
  else
    (E_S1 + E_S2 * (sqrt(1 + 2 * (e0 / E_S2) * ((t + t0) - t_S1) - 1)))-E_S0
  

  return(E_S)
  
}



##Soil water growth modifier
SWGmod<-function(SWconst,SWpower,MoistRatio){
  f_theta<-(1-(1-MoistRatio)^SWpower)/(1+(1-2*SWconst^SWpower)*(1-MoistRatio/SWconst)^SWpower)
  return(f_theta)
}