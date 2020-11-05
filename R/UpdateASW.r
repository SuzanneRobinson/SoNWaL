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
   #site[["MaxASW"]] <- 400-(parms[["sigma_2R"]]*state[["Wr"]]) 
     MaxASW <- site[["MaxASW"]] # calc from root biomass etc. thus MaxASW = 0.1*rbm*state[["Wr"]]
     
     #Almedia equations:
     sigma_zR=0.5
     V_nr=10
     A=20
     O_nr0= 500
     O_rz0=state[["ASW"]]
     K_s =0.1
     t = 31 #size of time-step?
     #rooting depth 
     z_r = (0.1*sigma_zR*state[["Wr"]])
     
     V_rz = z_r # not sure if this is equivalent
     
     #V_nr = A-z_r # shared area minus root zone?
     
     t_s0= (V_rz*V_nr)/(K_s*A*(V_rz+V_nr))
     
     O_rz= (((O_rz0*V_nr - O_nr0*V_rz)/(V_rz+V_nr))*exp(-t/t_s0))+
       ( V_rz/(V_rz+V_nr)*(O_rz0+O_nr0)) #re-arrange equation to get O_nz?
     
    excessSW <- max(ASWrain - EvapTransp - MaxASW, 0)
    state[["ASW"]] <- O_rz + ASWrain - EvapTransp - excessSW
    print(state[["ASW"]])
    #ASW is rainfall minus evap and excess SW
    #ASWrain is current ASW + rainfall + irrigation
    #moisture ratio is ASW/MaxASW
    #Lower maxASW means moisture ratio is higher (more water available) but max water moisture is lower?
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
