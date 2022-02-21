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

    RainIntcptn <- MaxIntcptn * ifelse(LAImaxIntcptn <= 0, 1, min(1, LAI/LAImaxIntcptn)) * Rain

    Qa <- parms[["Qa"]]
    Qb <- parms[["Qb"]]
    month <- weather[["Month"]]
    latitude <- site[["latitude"]]
    RAD.day <- state[["RAD.day"]]
    h <- GetDayLength(Lat = latitude, month = month)
    netRad <- max(Qa + Qb * (RAD.day * 1e+06 / h) * (1 - exp(-parms[["k"]] * LAI)), 0) # 1e+06 converts mega-joules to joules netrad is in W m^-2
    
    MinCond <- parms[["MinCond"]]
    MaxCond <- parms[["MaxCond"]]
    LAIgcx <- parms[["LAIgcx"]]
    fCg0 <- parms[["fCg0"]]
    CO2 <- site[["CO2"]]
    PhysMod <- state[["PhysMod"]]
    fCg <- fCg0/(1 + (fCg0 - 1) * CO2/350)
    CanCond <- (MinCond + (MaxCond - MinCond) * (min(1, LAI/LAIgcx))) * PhysMod * fCg
    CanCond <- ifelse(CanCond == 0, 1e-04, CanCond)
    #evaporation off leaves not affected by physical mods and stomatal response etc. 
    CanEvap <- (MinCond + (MaxCond - MinCond) * (min(1, LAI/LAIgcx)))
    CanEvap <- ifelse(CanEvap == 0, 1e-04, CanEvap)
    BLcond <- parms[["BLcond"]]
    VPD<-weather[["VPD"]] * exp(LAI * (-log(2)) / 5)

    Delta = 145 #Rate of change of saturation specific humidity with air temperature. (Pa/K−1)
    Cp = 1004 #specific heat cap of air J/kg-1
    Pa = 1.204 #Dry air density kg/m-3
    gamma = 66.1 # phsychrometric constant Pa/K-1
    lambda <- parms[["lambda"]]#energy for vapourisation of water J/kg-1

    CanTransp <- h*((CanCond*(Delta*netRad+BLcond*Pa*Cp*(VPD*1000)))/(lambda*((gamma+Delta)*CanCond+gamma*BLcond))) #Etransp/lambda * h
    
    Transp <- min(CanTransp*365/parms[["timeStp"]],state[["SWC_rz"]]) #general.info$daysinmonth[month] * CanTransp
    
    #intercepted rain evaporates based on solar radiation hitting leaves
    CanEvapSol <- h*((CanEvap*(Delta*netRad+BLcond*Pa*Cp*(VPD*1000)))/(lambda*((gamma+Delta)*CanEvap+gamma*BLcond))) #Etransp/lambda * h
    rIntercptEvap <- min((CanEvapSol*365/parms[["timeStp"]]),RainIntcptn) #general.info$daysinmonth[month] * CanTransp
    
    #non-Intercepted radiation
    interRad<-max((RAD.day * 1e+06/h),0)-max((RAD.day * 1e+06/h)*(1-exp(-parms[["k"]]*LAI)),0)
    
    state[["soilRad"]]<-interRad
    state[["totalRad"]]<-(RAD.day * 1e+06/h)
    
    ####Run using water balance sub-models from Almedia et al.####
    if (parms[["waterBalanceSubMods"]] == T) {

      
      #root depth assumed to be proportional to root biomass (see almedia and sands) .2 is probability of stones
      z_r = max(min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["maxRootDepth"]]),0.05)
      
      
      #calc volume of non-root zone
      vnr<-max(parms[["V_nr"]]-z_r,0)
      
      #derive the VOUMETRIC SWC in mm (theta_x vals in Almedia and Sands) for the rooting-zone soil profile at wp, fc and fsat in mm
      volSWC_wp= (parms[["wiltPoint"]])
      volSWC_sat= (parms[["satPoint"]])
      pseudo=parms[["pseudo"]]

      #Calculate accumulated soil evaporation for the month using soilEvap function
      throughFall<-Rain -rIntercptEvap
      evapRes <- soilEvap(parms, weather, state, interRad,h, throughFall,pseudo=pseudo)
      state[["potentialEvap"]]<-evapRes[[1]]
      
if(pseudo==T){
      evapSoil<-min(evapRes[[2]], (state[["SWC_rz"]]))
      #soil evaporation Minus monthly rainfall and irrigation gives cumulative E_S value - does not include drainage as that is loss from bottom of profile, this is top layer evaporation until drying at surface
      E_S <- max(evapRes[[1]]+evapRes[[2]] -throughFall - MonthIrrig,0) 
      state[["E_S"]] <- E_S}
if (pseudo==F){
     evapSoil<-min(evapRes[[1]], (state[["SWC_rz"]]))
     #soil evaporation Minus monthly rainfall and irrigation gives cumulative E_S value - does not include drainage as that is loss from bottom of profile, this is top layer evaporation until drying at surface
     E_S <- max(evapSoil -throughFall - MonthIrrig,0) 
     state[["E_S"]] <- E_S
}
      
      #Calculate drainage from root zone to non-root zone and out of non-root zone, where SWC of root zone exceeds field capacity of the soil zone (eq A.11)
      #Drainage parameter based on soil texture *24 to get from hourly to daily rates - use daily rates in these functions as they calc drainage of t days. 
      K_drain_rz=parms[["K_drain"]]#max((4.56*((state[["SWC_rz"]]/(z_r*1000))/parms[["satPoint"]])^18.4)*24,0.00001)
      K_drain_nrz=0.8#max((114*((state[["SWC_nr"]]/(z_r*1000))/parms[["satPoint"]])^15)*24,0.00001)

 
      rz_nrz_drain <- max(drainageFunc(parms, weather, 
                                       SWC=state[["SWC_rz"]],soilVol=z_r,
                                       K_drain=K_drain_rz), 0)
      nrz_out_drain <-max(drainageFunc(parms, weather,
                                       SWC=state[["SWC_nr"]],
                                       soilVol=vnr, K_drain=K_drain_nrz), 0)
      
      #Run soil water content function to get root zone SWC 
      SWC_rz <- soilWC(parms, weather, state,
                       K_s=parms[["K_s"]],
                       SWC=state[["SWC_rz"]],
                       soilVol=z_r)
      
      #volume of water moving from root zone to non-root zone is diff between current state of root zone SWC and updated root zone SWC
      rz_nrz_recharge<-state[["SWC_rz"]]-SWC_rz

      
      #update SWC by adding drainage from root zone and removing drainage out from non-root zone
      state[["SWC_nr"]] <- min(max(state[["SWC_nr"]] + rz_nrz_drain - nrz_out_drain + rz_nrz_recharge, 0),
                               volSWC_sat*vnr*1000)
      
      #Update root zone SWC with the addition of rainfall, irrigation, minus evap and drainage into non-root zone
      state[["SWC_rz"]] <- min(max(SWC_rz + Rain + MonthIrrig - rIntercptEvap - evapSoil-Transp - rz_nrz_drain,0),
                               volSWC_sat*z_r*1000)

      #to get total evaptranspiration use evapRes not E_S as this is modified by rainfall and so is amount soil has lost but not amount "evaporating" from soil
      #rainfall actually evaporated is not all that is intercepted during this time-step - some drips down to the soil over time
      EvapTransp <- min(Transp  + evapSoil,SWC_rz) +rIntercptEvap
      EvapTransp<-ifelse(EvapTransp < 0, 0, EvapTransp)

      volSWC_rz<-( state[["SWC_rz"]]/(z_r * 1000))
      
      #ASW calculated as SWC in the root zone divided by depth (mm) minus volumetric SWC of soil profile at wilting point
      # See Almedia and Sands eq A.2 and landsdown and sands 7.1.1
      state[["ASW"]] <-max(volSWC_rz - volSWC_wp, 0)
      state[["volSWC_rz"]]<-volSWC_rz
      
      
      scaleSW <-1#ifelse( (Transp+evapSoil)  + RainIntcptn== 0 ,1,EvapTransp/((Transp+evapSoil) + RainIntcptn))

      Etransp <- z_r
    } 
    if (parms[["waterBalanceSubMods"]] == F) {
      EvapTransp <- min(Transp + RainIntcptn, ASWrain)
      
      MaxASW<- site[["MaxASW"]]
      excessSW <- max(ASWrain - EvapTransp - MaxASW, 0)
      rz_nrz_drain<-excessSW
      state[["ASW"]] <-  max(ASWrain - EvapTransp - excessSW, 0)
      scaleSW <-EvapTransp/(Transp + RainIntcptn)
      Etransp <- 1
      nrz_out_drain<-1
      rz_nrz_drain<-1
    }
    
    GPP <- state[["GPP"]]
    NPP <- state[["NPP"]]
    state[["GPP"]] <- GPP * scaleSW
    state[["NPP"]] <- NPP * scaleSW
    state[c("RainIntcptn", "netRad", "CanCond", "Etransp", "CanTransp", 
            "Transp", "EvapTransp", "excessSW", "scaleSW","rz_nrz_recharge","excessSW_nr")] <- c(rIntercptEvap, 
                                                                 netRad, CanCond, Etransp, CanTransp, Transp, EvapTransp, 
                                                                 rz_nrz_drain, scaleSW,rz_nrz_recharge,nrz_out_drain)
    return(state)
  }
