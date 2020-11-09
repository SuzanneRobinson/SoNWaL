UpdateASW <-
function (state, weather, site, parms, general.info) #requires leaffall and leafgrow
{
    Rain <- weather[["Rain"]]
    MonthIrrig <- weather[["MonthIrrig"]]
    ASW <- state[["ASW"]]
    LAI <- state[["LAI"]] #calculate LAI with is.dormant - maybe not, because LAI was already calculated with 'predictVariablesInterest'
    ASWrain <- ASW + Rain + MonthIrrig
    LAImaxIntcptn <- parms[["LAImaxIntcptn"]]
    MaxIntcptn <- parms[["MaxIntcptn"]]
    RainIntcptn <- MaxIntcptn * ifelse(LAImaxIntcptn <= 0, 1, 
        min(1, LAI/LAImaxIntcptn)) * Rain
    Qa <- parms[["Qa"]]
    Qb <- parms[["Qb"]]
    month <- weather[["Month"]]
    latitude <- site[["latitude"]]
    RAD.day <- state[["RAD.day"]]
    h <- GetDayLength(Lat = latitude, month = month)
    netRad <- Qa + Qb * (RAD.day * 1e+06/h)
    MinCond <- parms[["MinCond"]]
    MaxCond <- parms[["MaxCond"]]
    LAIgcx <- parms[["LAIgcx"]]
    fCg0 <- parms[["fCg0"]]
    CO2 <- site[["CO2"]]
    PhysMod <- state[["PhysMod"]]
    #LAI <- state[["LAI"]] #calculate LAI with is.dormant - maybe not, because LAI was already calculated with 'predictVariablesInterest'
    fCg <- fCg0/(1 + (fCg0 - 1) * CO2/350)
    CanCond <- (MinCond + (MaxCond - MinCond) * (min(1, LAI/LAIgcx))) * 
        PhysMod * fCg
    CanCond <- ifelse(CanCond == 0, 1e-04, CanCond)
    e20 <- parms[["e20"]]
    rhoAir <- parms[["rhoAir"]]
    lambda <- parms[["lambda"]]
    VPDconv <- parms[["VPDconv"]]
    BLcond <- parms[["BLcond"]]
    VPD <- weather[["VPD"]]
    Etransp <- (e20 * netRad + rhoAir * lambda * VPDconv * VPD * 
        BLcond)/(1 + e20 + BLcond/CanCond)
    CanTransp <- Etransp/lambda * h
    Transp <- general.info$daysinmonth[month] * CanTransp
    EvapTransp <- min(Transp + RainIntcptn, ASWrain)
    MaxASW<-site[["MaxASW"]]

    
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
  
    #'@param t length of time-step
    #'@param V_nr volume of non-rooting zone
    #'@param V_rz volume of rooting zone
    #'@param A shared area
    #'@param z_r depth (or volume?) of rooting zone
    #'@param sigma_zR area/depth explored by 1kg of root biomass
    #'@param O_nr0 SWC of non-rooting zone at time 0
    #'@param O_rz0 SWC of rooting zone at time 0
    #'@param O_rz SWC of rooting zone
    
    #size of time-step. This is for monthly time steps
    t =  as.numeric(days_in_month(month)) 
    #Volume of initial non-rooting zone
    V_nr = parms[["V_nr"]]
    
    #soil water in rooting zone at t0, equiv ASW? depends if rain is assumed to occur at begining or end of month,
    #need to be careful so as not to mess up biomass allocation function which comes in after in the same time step (see RunModel.r)
    O_rz0 = state[["ASW"]]
    
    #Soil conductivity - see Landsberg book for more details on this
    K_s = parms[["K_s"]]
    
    #rooting depth / volume - Almedia describes this as depth, but sigma_zR could also be used to derive volume from root biomass?
    sigma_zR = parms[["sigma_zR"]]
    z_r = (0.1 * sigma_zR * state[["Wr"]])
    
    V_rz = z_r # not sure if this is equivalent, or whether there needs to be a conversion from depth to volume?
    
    #Shared area, equal to both volumes combined? Need to check this
    A = V_nr + V_rz
    
    #Time constant
    t_s0 = (V_rz * V_nr) / (K_s * A * (V_rz + V_nr))
    
    #Current state of soil water content in non-rooting zone (at time-0 for month)
    O_nr0 = state[["O_nr0"]]# should be dynamic? but not sure how to implement just yet...
    
    #State of soil water content in rooting zone at the end of the time step
    O_rz = (((O_rz0 * V_nr - O_nr0 * V_rz) / (V_rz + V_nr)) * exp(-t /t_s0)) +
      (V_rz / (V_rz + V_nr) * (O_rz0 + O_nr0)) #re-arrange equation to get O_nz?
    

#      ____     _ __  ____                             __  _          # 
#     / __/__  (_) / / __/  _____ ____  ___  _______ _/ /_(_)__  ___  # 
#    _\ \/ _ \/ / / / _/| |/ / _ `/ _ \/ _ \/ __/ _ `/ __/ / _ \/ _ \ # 
#   /___/\___/_/_/ /___/|___/\_,_/ .__/\___/_/  \_,_/\__/_/\___/_//_/ # 
#                               /_/                                   # 
       
       
    #'@param t_S1 duration of phase 1 evaporation
    #'@param e_S1 potential rate of evaporation, equal to wet surface evaporation - could be driven by temp?
    #- calc from penman montieth - also see landsdown book and discussion on methods by Choudhury, Zhang etc.
    #'@param E_S1 threshold - based on soil structure
    #'@param E_S2 how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
    #'@param E_S0 initial condition based on E_S
    #Calc using penman monteith eq. g_c is infinite - could be better modified by soil dryness etc.
    E_S1 = parms[["E_S1"]]
    E_S2 = parms[["E_S2"]]
    
    
    #Wet surface evaporation - using Penman Monteith eq.
    e_S1 <- max(((e20 * netRad + rhoAir * lambda * VPDconv * VPD *
                    BLcond) / (1 + e20 + BLcond / Inf)
    ), 0) #not sure about canopy transpiration?
    
    
    #E_S0 is E_S minus the amount of rainfall not intercepted and irigation
    #(assumes all rainfall occurs at the start of the month - need to implement finer scale timesteps!)
    E_S0 = max(site[["E_S"]] , 0)
    
    #Duration of phase 1 evaporation
    t_S1 = E_S0 / e_S1
    
    #Solved for t equation A.10 in Almedia, to get equivalent t for E_S0 value
    t0 = as.numeric(round(((-2 * E_S0 * E_S1) + (E_S0 ^ 2) + (2 * E_S0) +
                             (E_S1 ^ 2) - (2 * E_S1) + 1 + (2 * e_S1 * E_S2 * t_S1) - (E_S2 ^ 2)
    ) / (2 * e_S1 * E_S2)))
    
    
    #Integrate equation A.9 to get value at time t (assuming t is number of days in month)
    #and Calc E_S using t+t0 to get amount of evaporation between t0 and t
    E_S = if (t <= t_S1)
      e_S1 * t
    else
      (E_S1 + E_S2 * (sqrt(1 + 2 * (e_S1 / E_S2) * ((t + t0) - t_S1) - 1)))
    -
      E_S0
    
    Evaporation = ifelse(E_S <= 0, 0, E_S)
    
    #Minus monthly rainfall and irrigation from cumulative E_S value
    E_S <-
      E_S + RainIntcptn - Rain - MonthIrrig # rainfall needs ot be added after as biomass allocation function requires months rainfall
    
    site[["E_S"]] = ifelse(E_S <= 0, 0, E_S)
    
    #Calculate any excess soil water
    excessSW <- max(ASWrain - EvapTransp - MaxASW - Evaporation, 0)
    
    #Update available soil water
    state[["ASW"]] <-  max(O_rz + ASWrain - EvapTransp - excessSW - Evaporation,0)
    scaleSW <- EvapTransp/(Transp + RainIntcptn)
    GPP <- state[["GPP"]]
    NPP <- state[["NPP"]]
    state[["GPP"]] <- GPP * scaleSW
    state[["NPP"]] <- NPP * scaleSW
    state[c("RainIntcptn", "netRad", "CanCond", "Etransp", "CanTransp", 
        "Transp", "EvapTransp", "excessSW", "scaleSW")] <- c(RainIntcptn, 
        netRad, CanCond, Etransp, CanTransp, Transp, EvapTransp, 
        excessSW, scaleSW)
    return(state)
}
