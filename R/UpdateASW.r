UpdateASW <-
function (state, weather, site, parms, general.info) #requires leaffall and leafgrow
{
    Rain <- weather[["Rain"]]
    MonthIrrig <- weather[["MonthIrrig"]]
    ASW <- state[["ASW"]]
    LAI <- state[["LAI"]] #calculate LAI with is.dormant - maybe not, because LAI was already calculated with 'predictVariablesInterest'
    ASWrain <- ASW + Rain + MonthIrrig
    NRrain <- state[["O_nr0"]] + Rain + MonthIrrig
    
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
    print(Etransp)
    CanTransp <- Etransp/lambda * h
    Transp <- general.info$daysinmonth[month] * CanTransp
    EvapTransp <- min(Transp + RainIntcptn, ASWrain)
    EvapTranspNR <- min(Transp + RainIntcptn, NRrain)
    MaxASW<-site[["MaxASW"]]

     #Almedia equations:
     t =  as.numeric(days_in_month(month)) #size of time-step. This is for monthly time steps
     V_nr=parms[["V_nr"]]
     
     #soil water in rooting zone equiv ASW?
     O_rz0=state[["ASW"]]
     K_s =parms[["K_s"]]
     
     
     #rooting depth 
     sigma_zR=parms[["sigma_zR"]]
     z_r = (0.1*sigma_zR*state[["Wr"]])
     
     V_rz = z_r # not sure if this is equivalent, or whether there needs to be a conversion from depth to volume?
     
     #Shared area equal to both volumes combined? 
     A=V_nr+V_rz
     
     #Time constant
     t_s0= (V_rz*V_nr)/(K_s*A*(V_rz+V_nr))
     
     #Current state of soil water content in non-rooting zone (at time-0 for month)
     O_nr0 = state[["O_nr0"]]
     
     #State of soil water content in rooting zone at the end of the time step
     O_rz= (((O_rz0*V_nr - O_nr0*V_rz)/(V_rz+V_nr))*exp(-t/t_s0))+
       ( V_rz/(V_rz+V_nr)*(O_rz0+O_nr0)) #re-arrange equation to get O_nz?
     
    
     
     #soil water content at the end of month for non-rooting zone, to be used as starting point for next time step
   # state[["O_nr0"]] =((O_rz*V_rz)+(O_rz*V_nr)-(exp(-t/t_s0)*O_rz0*V_nr)-(O_rz0*V_rz))/((1-exp(-t/t_s0))*V_rz) #Not too sure about recharge of non-rooting area, re-arranged the equation for O_rz
     

   #  state[["O_nr0"]]=# min(O_nr0+(excessSW-O_rz),150)#ASWrain - EvapTransp - excessSW+(((O_nr0*V_rz-O_rz0*V_nr )/(V_rz+V_nr))*exp(-t/t_s0))+
      #( V_nr/(V_nr+V_rz)*(O_rz0+O_nr0))
     
#Evaporation
     #t_S1 = duration of phase 1 evaporation
     #e_S1 = potential rate of evaporation, equal to wet surface evaporation - could be driven by temp?
     #- calc from penman montieth - also see landsdown book and discussion on methods by Choudhury, Zhang etc.
     #E_S1 = threshold - based on soil structure
     #E_S2 = how quickly evaporation rate declines with accumulated phase 2 evaporation - based on soil structure
    #E_S0 = initial condition based on E_S
     #Calc using penman monteith eq. g_c is infinite - could be better modified by soil dryness etc.
     E_S1 =parms[["E_S1"]]
     E_S2 =parms[["E_S2"]]

     
     #Wet surface evaporation - using Penman Monteith eq. 
    e_S1 <- max(((e20 * netRad + rhoAir * lambda * VPDconv * VPD * 
              BLcond)/(1 + e20 + BLcond/Inf)),0) #not sure about canopy transpiration?


    #E_S0 is E_S minus the amount of rainfall not intercepted and irigation 
    E_S0 = max(site[["E_S"]]+RainIntcptn -Rain - MonthIrrig ,0)
    
    #Duration of phase 1 evaporation
    t_S1 = E_S0/e_S1
    
    #inverted equation A.10 in Almedia, solved for t, to get equivalent t from E_S0 value
    t0=as.numeric(round(((-2*E_S0*E_S1)+(E_S0^2)+(2*E_S0)+(E_S1^2)-(2*E_S1)+1+(2*e_S1*E_S2*t_S1)-(E_S2^2))/(2*e_S1*E_S2)))
    

    #Integrate equation A.9 to get value at time t and Calc E_S using t+t0 to get amount of evaporation between t0 and t
   E_S = if(t<=t_S1) e_S1*t else (E_S1 +E_S2*(sqrt(1+2*(e_S1/E_S2)*((t+t0)-t_S1)-1)))-E_S0
    Evaporation = ifelse(E_S<=0,0,E_S)
    site[["E_S"]]=ifelse(E_S<=0,0,E_S)
   #state[["O_nr0"]] =O_nrz CHECK TIME FOR WETTING EVENTS???
    
    excessSW <- max(ASWrain - EvapTransp - MaxASW-Evaporation, 0)
    excessSWNR <- max(NRrain - EvapTranspNR - MaxASW-Evaporation, 0)
     
    state[["ASW"]] <-  max(O_rz+ ASWrain - EvapTransp - excessSW - Evaporation,0)
    print( paste0("O_nr0 = ",site[["E_S"]], " ASW = ",state[["ASW"]]))
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
