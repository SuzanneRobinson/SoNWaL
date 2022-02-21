

##############likelihood for mixed prior regional calibrations####################



# Likelihood function
# function relies on the parameter list (paramList) defined at the submission of the mcmc,
# this list is what is submitted to the model run function to run SoNWal
# it contains elements for each parameter to be updated from the MCMC chain proposals
# additionally it contains the paramList$weather dataframe - this contains the longitudinal
# climate data and site specific data for each of the 13 regional sites
#'@param p matrix of parameter proposals for each chain, matrix: ncols = parameter number, nrows = number of internal chains
#'@return vector of likelihood values equal to number of rows in input p matrix
NLL_Reg_sitka_mixedP <- function(p) {
  #Sometimes number of rows in p matrix is less than number of chains (uncertain if algorithm design or bug)
  #this means that sometimes it comes in as a single row, basically a vector and needs transposing when converting to dataframe (or it does for the way i've done things)
  px <- if (is.null(nrow(p)) == F)
    data.frame(p) else
      data.frame(t(p))
  
  nm_all<-c(paste0("wiltPoint_Si",unique(paramList$weather$site)),
            paste0("fieldCap_Si",unique(paramList$weather$site)),
            paste0("satPoint_Si",unique(paramList$weather$site)),
            paste0("K_s_Si",unique(paramList$weather$site)),
            paste0("V_nr_Si",unique(paramList$weather$site)),
            paste0("E_S1_Si",unique(paramList$weather$site)),
            paste0("E_S2_Si",unique(paramList$weather$site)),
            paste0("shared_area_Si",unique(paramList$weather$site)),
            paste0("maxRootDepth_Si",unique(paramList$weather$site)),
            paste0("K_drain_Si",unique(paramList$weather$site)),
            paste0("startN_Si",unique(paramList$weather$site)),
            paste0("startC_Si",unique(paramList$weather$site)),
            "pFS2","pFS20","gammaF0","tgammaF","Rttover","mF","mR",
            "mS","Nf","Navm","Navx","klmax","krmax","komax","hc","qir","qil","qh","qbc","el","er","SWconst0",
            "SWpower0","sigma_zR"
            ,"aS","nS","pRx","pRn","gammaFx",
            "SLA0","SLA1","tSLA","alpha","Y","m0","MaxCond","LAIgcx","CoeffCond","BLcond",
            "k","Qa","Qb","MaxIntcptn")
  
  
  names(px) <- nm_all
  
  #create vector of -Inf vals to return if the model run function fails
  # fail <- rep(-Inf, nrow(px))
  #Get list of site names from the paramList$weather dataframe
  siteLst <- (unique(paramList$weather$site))
  #Update dataframe with each chains proposed param values repeated for each site
  siteVec <- rep(siteLst, each = nrow(px))
  px$chain <-
    rep(1:nrow(px))  #add chain number for reference and aggregating later
  px <- px[rep(seq_len(nrow(px)), length(siteLst)),]
  px$site <- siteVec
  
  
  
  #split px dataframe into list for running in parallel lapply function
  splitParams <- split(px, seq(nrow(px)))
  
  #update px dataframe with likelihood values for each site
  if(Sys.info()[1]!="Windows")
  { px$ll <-
    do.call(rbind, mcmapply(runMod_mixedP,splitParams,MoreArgs = list(paramList,nm_all),SIMPLIFY = F, mc.cores = 4))
  }else{
    px$ll <-
      do.call(rbind, mapply(runMod_mixedP,splitParams,MoreArgs = list(paramList,nm_all),SIMPLIFY = F))
  }
  #sum likelihood values by chain and return vector of likelihood values, one for each chain being run
  
  NlogLik <-
    as.vector(px %>% group_by(chain) %>% summarise(ll = sum(ll)) %>%
                pull(ll))
  
  
  
  NlogLik[is.na(NlogLik)==T]<--Inf
  return(NlogLik)
  
}





#runMod function
#'@param newParams list element containing dataframe of proposed parameter values for single chain for single site
runMod_mixedP <- function(newParams,paramListX,nm_all) {
  res <- tryCatch({
    
    #update parameter list with site specific params proposals
    siteSpecNm<-c(nm_all[grepl(paste0("\\_Si",newParams$site), nm_all)],
                  nm_all[!grepl("\\_Si", nm_all)])
    fitNm<-sub("_Si.*", "",siteSpecNm)
    
    
    paramListX[fitNm] <- newParams[siteSpecNm]
    #filter paramList$weather dataframe by site
    paramListX$weather <- filter(paramListX$weather, site == newParams$site)
    #get observed values (also contained in paramListX$weather dataframe) for site being fitted
    observed <-
      filter(paramListX$weather, is.na(mean_dbh_cm) == F) %>% group_by(Year) %>% summarise(dbh =
                                                                                             median(mean_dbh_cm),
                                                                                           dbhSD = median(dbhSD_cm))
    observed_N <-paramListX$weather$StemsPerHa[1]
    
    observed$dbhSD <-
      ifelse(observed$dbhSD == 0, 0.0001, observed$dbhSD)
    sY = min(observed$Year)
    eY = max(observed$Year)
    
    #filter climate data by planting year so model runs from planting year
    paramListX$weather <-
      filter(paramListX$weather, Year >= paramListX$weather$plantingYear[1])
    
    #run model
    output <- do.call(fr3PGDN, paramListX)
    
    #filter simulated data to match observed data format
    modelled <-
      output %>% group_by(Year) %>% summarise(dg = mean(dg)) %>% filter(Year >=
                                                                          sY & Year <= eY)
    modelled_N <-tail(output$N,1)
    
    
    #run likelihood function of observed vs simulated and get liklelihood value
    ifelse(
      any(is.na(modelled) == T),
      -Inf,
      flogL(
        data = c(observed$dbh,observed_N),
        sims = c(modelled$dg,modelled_N),
        data_s = c(observed$dbhSD,modelled_N*0.1)
      )
    )
  },
  error = function(cond) {
    #return -Inf if something goes wrong with param proposals
    return(-Inf)
  })
  
  res<-ifelse(max(output$LAI)>10,-Inf,res)
  return(res)
}


