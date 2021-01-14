CalculateModifiers <-
function (state, weather, site, parms, general.info) 
{
    Tmin <- parms[["Tmin"]]
    Tmax <- parms[["Tmax"]]
    Topt <- parms[["Topt"]]
    Tav <- weather[["Tmean"]]
    if (Tav < Tmin | Tav > Tmax) {
        fT <- 0
    }
    else {
        fT <- ((Tav - Tmin)/(Topt - Tmin)) * ((Tmax - Tav)/(Tmax - 
            Topt))^((Tmax - Topt)/(Topt - Tmin))
    }
    kF <- parms[["kF"]]
    FrostDays <- weather[["FrostDays"]]
    fF <- 1 - kF * (FrostDays/parms[["timeStp"]])
    CO2 <- site[["CO2"]]
    fCalphax <- parms[["fCalphax"]]
    fCalpha <- fCalphax * CO2/(350 * (fCalphax - 1) + CO2)
    ## fN0 <- parms[["fN0"]]
    ## FR <- state[["FR"]]
    ## fN <- FR #fN0 + (1 - fN0) * FR
    t <- state[["t"]]
    RelAge <- t/parms[["MaxAge"]]
    rAge <- parms[["rAge"]]
    nAge <- parms[["nAge"]]
    fAge <- 1/(1 + (RelAge/rAge)^nAge)
    CoeffCond <- parms[["CoeffCond"]]
    VPD <- weather[["VPD"]]
    fVPD <- exp(-CoeffCond * VPD)
    

    
    parms.soil <- general.info$parms.soil
    parms.sw.site <- parms.soil[which(parms.soil$soilclass == 
        site[["soilclass"]]), ]
    SWconst <- parms.sw.site[["SWconst"]]
    SWpower <- parms.sw.site[["SWpower"]]
    ASW <- state[["ASW"]]

    ##change moistratio and soil water growth mod if using updated sub-models
    if(parms[["waterBalanceSubMods"]]==T){

     # z_r = min((0.1 * parms[["sigma_zR"]] * state[["Wr"]]),parms[["maxRootDepth"]])
      #calc rooting zone soil profile VOLUMETRIC SWC at wp and fc (mm)
      volSWC_wp= parms[["wiltPoint"]]
      volSWC_fc= parms[["fieldCap"]]
      MaxASW <- (volSWC_fc-volSWC_wp)
      
      #calc moist ratio
      MoistRatio<- ASW/MaxASW
      #modify MoistRatio if numerators are above or below certain values (see Landsberg and waring)
      MoistRatio<-ifelse(ASW>=0,MoistRatio,0)
      MoistRatio<-ifelse(ASW>MaxASW,1,MoistRatio)
      fSW<-SWGmod(SWconst,SWpower,MoistRatio)
    
    }else
    {
      MaxASW <- site[["MaxASW"]]
      MoistRatio<-ASW/MaxASW
      fSW <- 1/(1 + ((1 - MoistRatio)/SWconst)^SWpower)
    }
    
    
    PhysMod <- fAge * min(fVPD, fSW)
    ## state[c("fT", "fF", "fCalpha", "fN", "fAge", "fVPD", "fSW", 
    ##     "PhysMod")] <- c(fT, fF, fCalpha, fN, fAge, fVPD, fSW, 
    ##                      PhysMod)
    state[c("fT", "fF", "fCalpha", "fAge", "fVPD", "fSW", 
            "PhysMod")] <- c(fT, fF, fCalpha, fAge, fVPD, fSW, 
                             PhysMod)
    return(state)
}
