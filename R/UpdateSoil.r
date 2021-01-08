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
    klmax <- parms[["klmax"]]*12/parms[["timeStp"]] # switch from monthly rates to time-step rates
    krmax <- parms[["krmax"]]*12/parms[["timeStp"]]
    komax <- parms[["komax"]]*12/parms[["timeStp"]]
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
      #MaxASW <- site[["MaxASW"]]
      #MoistRatio <- ASW/MaxASW
      Tmin <- parms[["Tmin"]]
      Tmax <- parms[["Tmax"]]
      Topt <- parms[["Topt"]]
      Tav <- weather[1,"Tmean"]
      
      dg<-((state[["Wsbr"]]*1000/state[["N"]])/parms[["aS"]])^(1/parms[["nS"]])
      
      
      ##change moistratio and soil water growth mod if using updated sub-models
      if(parms[["waterBalanceSubMods"]]==T){
        
        #calc wp and fc
        theta_wp= parms[["theta_wp"]]*1000
        theta_fc= parms[["theta_fc"]]*1000
        MaxASW <- (theta_fc-theta_wp)
        
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
