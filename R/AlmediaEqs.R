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
#'@param SWC_nr0 SWC of non-rooting zone at time 0
#'@param SWC_rz0 SWC of rooting zone at time 0

#'@return SWC_rz SWC of rooting zone

soilWC<-function(parms,weather,state){
  
  #size of time-step. 
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
  
  #Volume of initial non-rooting zone
  V_nr = parms[["V_nr"]]
  
  #soil water in rooting zone at t0 (start of time-step)
  SWC_rz0 = state[["SWC_rz"]]
  
  #Soil conductivity - see Landsberg book for more details on this 
  K_s = parms[["K_s"]]
  
  #rooting depth / volume - Almedia describes this as depth, assumed to be proportional in paper to biomass
  z_r = min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["maxRootDepth"]]) # can't go deeper than non-rooting zone/max root depth
  z_r=min(z_r,V_nr)
  
  V_rz = z_r #Almedia and Sands paper suggests volume of root zone is equivalent to z_r
  
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
 
  return(SWC_rz)
}


#'@description Calculates drainage out of a soil zone
#'@param parms global 3PG parameter values
#'@param weather global weather values
#'@param SWC soil water content of zone to drain
#Internal function variable descriptions
#'@param t length of time-step
#'@return amount of drainage

drainageFunc<-function(parms,weather,SWC,soilVol){
  
  #size of time-step. This is for monthly time steps
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
 
  #Volumetric field capacity
  volSWC_fc= parms[["fieldCap"]]
  
  #Drainage parameter based on soil texture
  K_drain<-parms[["K_drain"]]
  
  ##calculate drainage, convert soilVol to mm
  Qd<-(SWC-(volSWC_fc*soilVol*1000))*(1-exp(-K_drain*t))
  return(Qd)
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


soilEvap<-function(parms,weather,state,interRad,h){
  
  e20 <- parms[["e20"]]
  rhoAir <- parms[["rhoAir"]]
  lambda <- parms[["lambda"]]#Volumetric latent heat of vaporization. Energy required per water volume vaporized (see penman monteith)
  VPDconv <- parms[["VPDconv"]]
  VPD <- weather[["VPD"]]
  E_S1 = (parms[["E_S1"]])
  E_S2 = (parms[["E_S2"]])

  
  #size of time-step. This is for monthly time steps
  if (parms[["timeStp"]] ==12) t =   days_in_month(weather[["Month"]]) 
  if (parms[["timeStp"]] ==52) t =   7
  if (parms[["timeStp"]] ==365) t =  1
  if (is.numeric(t)==F) print ("unsupported time step used")
  #parameter values from landsberg and sands book (7.2.1)
 # s<-145
 # gb_s<-0.01#soil boundary layer conductance
 # g_C<-Inf
 # P_a = 1.204
 # C_pa = 1004
 # gamma= 66.1
 # s=145
 # D=VPD*100
 # lambda=2.45
  
  soilBoundaryCond<-0.01
  soilCond<-Inf
  #convert joules m^2 per hour to watts m^2 (1 Wm^2 = 1 J m^2 per second)

  #get potential evaporation rate using penman-monteith with soil specific params - following eq in landsberg and sands (7.2.1), adjusted to output evap volume rate
  #e0<-h*(((s*interRad+gb_s*P_a*C_pa*D)/(s+gamma*(1+gb_s/g_C)))/lambda)
 # print(paste0("eox= ",e0x))
  e0<- (e20 * interRad + rhoAir * lambda  * VPD * 
         soilBoundaryCond)/(1 + e20 + soilBoundaryCond/soilCond)
  e0<-max(e0/lambda * h,0)
 # 
  #print(paste0("eo= ",e0))
  
  #e0<-max(h*(soilCond*(145*interRad+soilBoundaryCond*1.204*1004*((VPD*100)*VPDconv))/(2.45*((66.1+145)*soilCond+66.1*soilBoundaryCond))),1)
  #E_S0 is E_S at the start of the time-step
  E_S0 = state[["E_S"]]

  #Duration of phase 1 evaporation
 t_S1 = E_S1 / e0
 #Solved for t equation A.10 in Almedia, to get equivalent t for E_S0 value
 t0 = as.numeric(round(((-2 * E_S0 * E_S1) + (E_S0 ^ 2) + (2 * E_S0) +
                          (E_S1 ^ 2) - (2 * E_S1) + 1 + (2 * e0 * E_S2 * t_S1) - (E_S2 ^ 2)
 ) / (2 * e0 * E_S2)))
 #Integrate equation A.9 to get value at time t (assuming t is number of days in month)
 #and Calc E_S using t+t0 to get amount of evaporation between t0 and t
 
 
 
 E_S = if ((t + t0) <= t_S1)
   e0 * (t + t0)
 else
   (E_S1 + E_S2 * (sqrt(1 + 2 * (e0 / E_S2) * ((t + t0) - t_S1) - 1)))-E_S0
  
 # E_S<-E_S0+ifelse(E_S0<=E_S1,e0,e0/(1+(E_S0-E_S1)/E_S2))

  return(list(E_S,e0))
  
}



##Soil water growth modifier
SWGmod<-function(SWconst,SWpower,MoistRatio){
  f_theta<-(1-(1-MoistRatio)^SWpower)/(1+(1-2*SWconst^SWpower)*(1-MoistRatio/SWconst)^SWpower)
  return(f_theta)
}