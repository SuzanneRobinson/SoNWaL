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
    #Penman-monteith 
    Etransp <- ((e20 * netRad + rhoAir * lambda  * VPD * 
       BLcond)/(1 + e20 + BLcond/CanCond))#*12/parms[["timeStp"]]
    #CanTransp <- h*(CanCond*(145*netRad+BLcond*1.204*1004*(VPD*100))/(2.45*((66.1+145)*CanCond+66.1*BLcond))) #Etransp/lambda * h
    CanTransp <- Etransp/lambda * h
    
    #CanTransp<-h*(CanCond*(145*netRad+BLcond*1.204*1004*(VPD*1000))/(2.45*((66.1+145)*CanCond+66.1*BLcond)))
    #less accurate but easier to have flexible time-steps?
    Transp <- CanTransp*365/parms[["timeStp"]] #general.info$daysinmonth[month] * CanTransp

    #non-Intercepted radiation
    interRad<-(RAD.day * 1e+06/h)-max((RAD.day * 1e+06/h)*(1-exp(-parms[["k"]]*LAI)),0)

    state[["soilRad"]]<-interRad
    state[["totalRad"]]<-(RAD.day * 1e+06/h)
    
        ####Run using water balance sub-models from Almedia et al.####
    if (parms[["waterBalanceSubMods"]] == T) {
      
      #root depth assumed to be proportional to root biomass (see almedia and sands)
      z_r = min((0.1 * parms[["sigma_zR"]] * state[["Wr"]])/0.2,parms[["maxRootDepth"]])
      
      #derive the VOUMETRIC SWC in mm (theta_x vals in Almedia and Sands) for the rooting-zone soil profile at wp, fc and fsat in mm
      volSWC_wp= (parms[["wiltPoint"]])
      volSWC_fc= (parms[["fieldCap"]])
      volSWC_sat= (parms[["satPoint"]])
      
      #Run soil water content function to get root zone SWC
      SWC_rz <- soilWC(parms, weather, state,K_s=0.5,SWC=state[["SWC_rz"]],soilVol=z_r)
      
      #Calculate accumulated soil evaporation for the month using soilEvap function
      evapRes <- soilEvap(parms, weather, state, interRad,h)
      state[["potentialEvap"]]<-evapRes[[2]]
      
      #soil evaporation Minus monthly rainfall and irrigation gives cumulative E_S value - does not include drainage as that is loss from bottom of profile, this is topy layer evaporation until drying at surface
      E_S <- max(evapRes[[1]] + RainIntcptn - Rain - MonthIrrig,0) 
      state[["E_S"]] <- E_S
      
      #to get total evaptranspiration use evapRes not E_S as this is modified by rainfall and so is amount soil has lost but not amount "evaporating" from soil
      EvapTransp <- min(Transp +  evapRes[[1]] , state[["SWC_rz"]]) 
      EvapTransp<-ifelse(EvapTransp<0,0,EvapTransp)
      #fsmod for sonwal app (not used in fitting or anything)
      #EvapTransp <- if(state[["t"]]>46) EvapTransp*parms[["fsMod"]] else EvapTransp
      
      #Calculate drainage from root zone to non-root zone and out of non-root zone, where SWC of root zone exceeds field capacity of the soil zone (eq A.11)
      
  
      #Drainage parameter based on soil texture
      K_drain_rz=max((18.6*((state[["SWC_rz"]]/(z_r*1000))/parms[["satPoint"]])^13.2)*24,0.00001)
      K_drain_nrz=max((114*((state[["SWC_nr"]]/(z_r*1000))/parms[["satPoint"]])^15)*24,0.00001)
      
      rz_nrz_drain <- max(drainageFunc(parms, weather, SWC=state[["SWC_rz"]],soilVol=z_r,K_drain=K_drain_rz),0)
      nrz_out_drain <-max(drainageFunc(parms, weather, SWC=state[["SWC_nr"]],soilVol=parms[["V_nr"]]-z_r,K_drain=K_drain_nrz),0)
      
      #volume of water moving from non-root zone to root zone is diff between current state of root zone SWC and updated root zone SWC
      rz_nrz_recharge<-state[["SWC_rz"]]-SWC_rz
      
      #update SWC by adding drainage from root zone and removing drainage out from non-root zone
      state[["SWC_nr"]] <- max(state[["SWC_nr"]]+rz_nrz_drain-nrz_out_drain+rz_nrz_recharge,0)
      #Update root zone SWC with the addition of rainfall, irrigation, minus evap and drainage into non-root zone
      state[["SWC_rz"]] <- max(SWC_rz + Rain + MonthIrrig - EvapTransp - rz_nrz_drain-RainIntcptn,0)
      

      #Volumetric SWC of rooting zone (z_r converted to mm as SWC in mm), equiv to SWC fraction in observed data - water to soil ratio
      volSWC_rz<-min(( state[["SWC_rz"]] /(z_r*1000)),volSWC_fc)
      # 
     # print(paste0("volSWC_rz= ",volSWC_rz))
      #ASW calculated as SWC in the root zone divided by depth (mm) minus volumetric SWC of soil profile at wilting point
      # See Almedia and Sands eq A.2 and landsdown and sands 7.1.1
      state[["ASW"]] <-max(volSWC_rz-volSWC_wp,0)
      state[["volSWC_rz"]]<-volSWC_rz
      scaleSW <-1
      
    } 
    if (parms[["waterBalanceSubMods"]] == F) {
      EvapTransp <- min(Transp + RainIntcptn, ASWrain)
      
      MaxASW<- site[["MaxASW"]]
      excessSW <- max(ASWrain - EvapTransp - MaxASW, 0)
      rz_nrz_drain<-excessSW
      state[["ASW"]] <-  max(ASWrain - EvapTransp - excessSW, 0)
      scaleSW <-EvapTransp/(Transp + RainIntcptn)
      
      }
    
    
    #print(scaleSW)
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
