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
    netRad <- Qa + Qb * (RAD.day * 1e+06/h)*(1-exp(-parms[["k"]]*LAI)) # 1e+06 converts mega-joules to joules
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
    #less accurate but easier to have flexible time-steps, could include set useable time-steps and use if statements?
    Transp <- CanTransp*365/parms[["timeStp"]] #general.info$daysinmonth[month] * CanTransp
    EvapTransp <- min(Transp + RainIntcptn, ASWrain)
    
   


    #non-Intercepted radiation
    interRad<-(RAD.day * 1e+06/h)-max(netRad,0)
#    .GlobalEnv$interRad <- rbind(.GlobalEnv$interRad,interRad)
#    .GlobalEnv$Rad <- rbind(.GlobalEnv$Rad,(RAD.day * 1e+06/h))
    
    ####Run using water balance sub-models from Almedia et al.####
    if (parms[["waterBalanceSubMods"]] == T) {
      
      #derive theta values, rooting zone can't be larger than non-rooting zone
      theta_wp= parms[["theta_wp"]]*min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["V_nr"]])*1000
      theta_fc= parms[["theta_fc"]]*min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["V_nr"]])*1000
      MaxASW <- theta_fc-theta_wp

      
      #Run soil water content function to get root zone SWC
      sigma_rz <- soilWC(parms, weather, state)
      
      #Calculate accumulated soil evaporation for the month using soilEvap function
      E_S <- soilEvap(parms, weather, state, interRad)
      
      #Minus monthly rainfall and irrigation from cumulative E_S value
      E_S <- E_S + RainIntcptn - Rain - MonthIrrig # rainfall needs ot be added after as biomass allocation function requires months rainfall
     
      Evaporation = ifelse(E_S <= 0, 0, E_S)
      state[["E_S"]] = ifelse(E_S <= 0, 0,  E_S)
      
      #Calculate any excess soil water use max ASW of field cap?
      excessSW <- max(sigma_rz + Rain + MonthIrrig - EvapTransp - theta_fc - Evaporation, 0)
      
      #Update here, to have soil water go from rooting zone to non-rooting zone, use K_s and rooting zone depth?
      #state[["sigma_nr0"]]<- max(state[["sigma_nr0"]] + Rain + MonthIrrig - excessSW - Evaporation, 0)
      state[["MaxASW_state"]] <-MaxASW

      #Update value for non-rooting zone SWC - assuming rainfall, evaporation and evapotranspiration occurs in the rooting zone
      state[["sigma_nr0"]] <-soilWC_NRZ(parms, weather, state)
      
      print("-")
      print(paste0("rz = ", state[["ASW"]]))
      print(paste0("nrz = ", state[["sigma_nr0"]]))
      
      #Update available soil water -CHECK IF MIN OF 0 NEEDED
      state[["ASW"]] <- max(sigma_rz + Rain + MonthIrrig - EvapTransp - excessSW - Evaporation, 0)
      print(paste0("rz2 = ", state[["ASW"]]))
      print("-")
      
    } else
    {
      MaxASW<- site[["MaxASW"]]
      excessSW <- max(ASWrain - EvapTransp - MaxASW, 0)
      state[["ASW"]] <-  max(ASWrain - EvapTransp - excessSW, 0)
    }
    
    
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
