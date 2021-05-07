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
    netRad <- Qa + Qb * (RAD.day * 1e+06/h)*(1-exp(-parms[["k"]]*LAI)) # 1e+06 converts mega-joules to joules netrad is in W m^-2

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
    Etransp <- ((e20 * netRad + rhoAir * lambda * VPDconv * VPD * 
        BLcond)/(1 + e20 + BLcond/CanCond))*12/parms[["timeStp"]]
    CanTransp <- Etransp/lambda * h
    #less accurate but easier to have flexible time-steps?
    Transp <- CanTransp*365/parms[["timeStp"]] #general.info$daysinmonth[month] * CanTransp
    EvapTransp <- min(Transp + RainIntcptn, state[["SWC_rz"]]+Rain-RainIntcptn) # NOT SURE IF THIS NEEDS A MINIUM VALUE???? if so, ASW+rain or SWC+rain? remember ASW is now volumetric
    #evap from soil and drainage occur after transpiration so max is just plus rain...?
    
    #non-Intercepted radiation
    interRad<-(RAD.day * 1e+06/h)-max((RAD.day * 1e+06/h)*(1-exp(-parms[["k"]]*LAI)),0)
    state[["soilRad"]]<-interRad
    state[["totalRad"]]<-(RAD.day * 1e+06/h)
    
        ####Run using water balance sub-models from Almedia et al.####
    if (parms[["waterBalanceSubMods"]] == T) {
      
      #root depth assumed to be proportional to root biomass (see almedia and sands)
      z_r = min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["maxRootDepth"]])
      
      #derive the VOUMETRIC SWC in mm (theta_x vals in Almedia and Sands) for the rooting-zone soil profile at wp, fc and fsat in mm
      volSWC_wp= (parms[["wiltPoint"]])
      volSWC_fc= (parms[["fieldCap"]])
      volSWC_sat= (parms[["satPoint"]])
      
      #Run soil water content function to get root zone SWC
      SWC_rz <- soilWC(parms, weather, state)
      
      #Calculate accumulated soil evaporation for the month using soilEvap function
      evapRes <- soilEvap(parms, weather, state, interRad,h)
      E_S<-evapRes[[1]]
      state[["potentialEvap"]]<-evapRes[[2]]
      
      #evaporation Minus monthly rainfall and irrigation from cumulative E_S value
      E_S <- max(E_S + RainIntcptn - Rain - MonthIrrig,0) # rainfall needs to be added after as biomass allocation function requires months rainfall
      state[["E_S"]] = E_S
     
      #Calculate drainage from root zone to non-root zone and out of non-root zone, where SWC of root zone exceeds field capacity of the soil zone (eq A.11)
      rz_nrz_drain <- max(drainageFunc(parms, weather, SWC=state[["SWC_rz"]],soilVol=z_r),0)
      nrz_out_drain <-max(drainageFunc(parms, weather, SWC=state[["SWC_nr"]],soilVol=parms[["V_nr"]]-z_r),0)
      
      #volume of water moving from non-root zone to root zone is diff between current state of root zone SWC and updated root zone SWC
      rz_nrz_recharge<-state[["SWC_rz"]]-SWC_rz
      
      #update SWC by adding drainage from root zone and removing drainage out from non-root zone
      state[["SWC_nr"]] <- max(state[["SWC_nr"]]+rz_nrz_drain-nrz_out_drain+rz_nrz_recharge,0)
      #Update root zone SWC with the addition of rainfall, irrigation, minus evap and drainage into non-root zone
      state[["SWC_rz"]] <- max(SWC_rz + (Rain + MonthIrrig - EvapTransp - rz_nrz_drain - E_S),0)
      

      #Volumetric SWC of rooting zone (z_r converted to mm as SWC in mm), equiv to SWC fraction in observed data - water to soil ratio
      volSWC_rz<-(SWC_rz/(z_r*1000))
      
      #ASW calculated as SWC in the root zone divided by depth (mm) minus volumetric SWC of soil profile at wilting point
      # See Almedia and Sands eq A.2 and landsdown and sands 7.1.1
      state[["ASW"]] <-max(volSWC_rz-volSWC_wp,0)
      state[["volSWC_rz"]]<-volSWC_rz
      
    } else
    {
      MaxASW<- site[["MaxASW"]]
      excessSW <- max(ASWrain - EvapTransp - MaxASW, 0)
      rz_nrz_drain<-excessSW
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
        rz_nrz_drain, scaleSW)
    return(state)
}
