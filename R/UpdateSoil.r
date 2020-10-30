UpdateSoil <-
function(state, parms, site, general.info, weather){
    GPP <- state[["GPP"]]
    NPP <- state[["NPP"]]
    Wdl <- state[["Wdl"]]
    Wds <- state[["Wds"]]
    Wdr <- state[["Wdr"]]
    YlC <- state[["YlC"]]
    YrC <- state[["YrC"]]
    OC <- state[["OC"]]
    YlN <- state[["YlN"]]
    YrN <- state[["YrN"]]
    ON <- state[["ON"]]
    Nav <- state[["Nav"]]
    Rs <- state[["Rs"]]
    difWl <- state[["difWl"]]
    difLitter <- state[["difLitter"]]
    difRoots <- state[["difRoots"]]
    dg <- state[["dg"]]
    fT <- state[["fT"]]
    fSW <- state[["fSW"]]
    klmax <- parms[["klmax"]]
    krmax <- parms[["krmax"]]
    komax <- parms[["komax"]]
    hc <- parms[["hc"]]
    qir <- parms[["qir"]]
    qil <- parms[["qil"]]
    qbc <- parms[["qbc"]]
    qh <- parms[["qh"]]
    el <- parms[["el"]]
    er <- parms[["er"]]
    Nf <- parms[["Nf"]]
    Navm <- parms[["Navm"]]
    Navx <- parms[["Navx"]]
    ## Initialise some variables when t<-0
    if(state[["t"]]==0){
        parms.soil <- general.info$parms.soil
        parms.sw.site <- parms.soil[which(parms.soil$soilclass == 
        site[["soilclass"]]), ]

        SWconst <- parms.sw.site[["SWconst"]]
        SWpower <- parms.sw.site[["SWpower"]]
        ASW <- state[["ASW"]] 
        MaxASW <- site[["MaxASW"]] #calc MaxASW from rooting depths?
        MoistRatio <- ASW/MaxASW
        Tmin <- parms[["Tmin"]]
        Tmax <- parms[["Tmax"]]
        Topt <- parms[["Topt"]]
        Tav <- weather[1,"Tmean"]

        dg<-((state[["Wsbr"]]*1000/state[["N"]])/parms[["aS"]])^(1/parms[["nS"]])
      #Current soil water modifier - update?
          fSW <- 1/(1 + ((1 - MoistRatio)/SWconst)^SWpower)
        #Es1 threshold transpiration
        #Current cumulative value of evap
        #SWconst equiv of Es2 in Almedia?
        #May need state var to keep track of Es?
        #ASW will be latest amount of rainfall throughfall
        #Minus ASW from cumulative Es?
        # 1. get current Es cumulative value
        # 2. Es minus rainfall in mm from wetting event
        # 3. invert equation A10 in paper to get equiv start time (equivalent to current cumulative rainfall value)
        #fSW <- 1/(1 + ((Es-Es1)/SWconst))
        # or just use integrated equation? 
       # currently climate data is monthly, use pseudo-daily time steps?   
       #   t_R = time to the next rainfall event,
       # For monthly or pseudo-daily time steps, t_R is days in
       # the month/number of rain days, and soil evaporation is
       # computed separately for each rainfall event and then
       # summed over all events.
        
        
        if (Tav < Tmin | Tav > Tmax) {
            fT <- 0
        }
        else {
            fT <- ((Tav - Tmin)/(Topt - Tmin)) * ((Tmax - Tav)/(Tmax - 
                                                                Topt))^((Tmax - Topt)/(Topt - Tmin))
        }
        Wds <- 0
        Wdl <- 0
        Wdr <- 0
        difLitter <- 0
        difRoots <- 0
        coarseDifRoots <- 0
        fineDifRoots <- 0
        difWl <- 0
    }
    ## Calculate the fineCoarseRootRatio based on DBH    
    fineCoarseRatio <- 0.1276+1462.2671*exp(-1.7958*dg) ## Empirical relationship for Scots Pine - GX !!! Requires a revised method !!!
    coarseDifRoots <- difRoots / (1 + fineCoarseRatio)
    fineDifRoots <- difRoots - coarseDifRoots
    ##Calculate the decomposition rate    
    kr <- krmax * fSW * fT
    kl <- klmax * fSW * fT
    ko <- komax * fSW * fT  
    ##Calculate fluxes in, out and between carbon and nitrogen pools
    ##Carbon fluxes
    YlCflx <- kl * (1 - hc) * YlC
    YrCflx <- kr * (1 - hc) * YrC
    OCflx <- ko * OC   
    ##Humification coefficients
    hl <- kl * hc * YlC
    hr <- kr * hc * YrC
    ##Humification coefficients for nitrogen pools
    hNl <- kl * hc * (YlN / qh)
    hNr <- kr * hc * (YrN / qh)
    ##Nitrogen fluxes
    YlNflx <- kl * ((1 - hc) / (1 - el)) * (YlN - el * (YlC / qbc))
    YlNflx <- max(0.,YlNflx)
    YrNflx <- kr * ((1 - hc) / (1 - er)) * (YrN - er * (YrC / qbc))
    YrNflx <- max(0.,YrNflx)
    ONflx <- ko * ON
    ONflx <- max(0.,ONflx)
    ##Calculate heterotrophic respiration
    Rs <- YlCflx + YrCflx + OCflx
    ##Now calculate carbon and nitrogen pools
    YrC <- YrC + ( (Wds + coarseDifRoots) / 2) - YrCflx - hr
    YlC <- YlC + ((difLitter + fineDifRoots + Wdl + Wdr) / 2) - YlCflx - hl
    OC <- OC + hl + hr - OCflx
    YrN <- YrN + ((Wds + coarseDifRoots) / (2 * qir)) - YrNflx - hNr
    YlN <- YlN + ((difLitter + fineDifRoots + Wdl + Wdr) / (2 * qil)) - YlNflx - hNl
    ON <- ON + hNr + hNl - ONflx
    ## Calculate the total pools
    totC <- YrC + YlC + OC
    totN <- YrN + YlN + ON
    ## Calculate the Fertility Rating
    ## First calculate available nitrogen
    Navflx <- YrNflx + YlNflx + ONflx
    ## Calculate nitrogen uptake
    Un <- difWl * Nf
    ## Update available nitrogen pool
    Nav <- Nav + Navflx - Un
    Nav <- ifelse(is.na(max( Navm, Nav )),0,max( Navm, Nav))  ## A small change to accommodate Bayesian calibration and avoid NAs
    if( Nav > Navx){
        Nleach <- Nav - Navx
        Nav <- Navx
    }else{
        Nleach <- 0
    }
    ## Now estimate fN
    fN <- (Nav - Navm) / (Navx - Navm)

    ## Calculate ecosystem scale fluxes
    ## now that we know soil respiration
    if(state[["t"]]!=0){
        NEE <- Rs - NPP
        Ra <- GPP - NPP
        Reco <- Ra + Rs
    }else{
        NEE <- 0
        Ra <- 0
        Reco <-0
    }
    
    ## Export the results
    state[c("YrC","YlC","OC","YrN","YlN","ON",
            "kl","kr","ko","hl","hr","hNl","hNr",
            "YrCflx","YlCflx","OCflx","YrNflx","YlNflx","ONflx","Navflx",
            "totC","totN","Un","Nav","Nleach","fN",
            "NEE","Reco","Ra","Rs")] <-
        c(YrC,YlC,OC,YrN,YlN,ON,
          kl,kr,ko,hl,hr,hNl,hNr,
          YrCflx,YlCflx,OCflx,YrNflx,YlNflx,ONflx,Navflx,
          totC,totN,Un,Nav,Nleach,fN,
          NEE,Reco,Ra,Rs)
    return(state)
}
