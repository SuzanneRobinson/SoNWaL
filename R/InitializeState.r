InitializeState <-
function (stand, site, parms, general.info, weather) 
{
    names.state <- c("Year", "Month", "t", "hdom", "N", "G", 
                     "dg", "Vu", "LAI", "Wl", "Wr", "Wsbr", "rotation", "cycle", 
                     "Nharv", "rm.sprouts", "Ww", "Wb", "Wbr", "Wa", "W", 
                     "wsbrg", "ASW", "Wlitt", "RAD.day", "RAD", "CanCover", 
                     "lightIntcptn", "phi.p", "phi.pa", "fT", "fF", "fCalpha", 
                     "fN", "fAge", "fVPD", "fSW", "PhysMod", "alphaC",
                     "GPP", "NPP", "NEE", "Reco", "Ra", "Rs",
                     "RainIntcptn", "netRad", "fCg", "CanCond", "Etransp", 
                     "CanTransp", "Transp", "EvapTransp", "excessSW", "scaleSW", 
                     "pR", "pFS", "pS", "pF", "difWl", "difWr", "difWsbr", 
                     "Littfall", "difLitter", "difRoots", "TotalLitter", "Ndead", 
                     "Wdl", "Wds", "Wdr","YrC","YlC","OC","YrN","YlN","ON","kl",
                     "kr","ko","hl","hr","hNl","hNr","YrCflx","YlCflx","OCflx",
                     "YrNflx","YlNflx","ONflx","Navflx","totC","totN","Un","Nav",
                     "Nleach","SWC_nr","MaxASW_state","E_S","SWC_rz","volSWC_rz","soilRad","totalRad","potentialEvap")
    .length.state.vec <- length(names.state)
    state.vector <- rep(NA, .length.state.vec)
    names(state.vector) <- names.state
    state.init <- state.vector
    state.init[names(stand)] <- stand
    state.init["ASW"] <- site[["ASW"]]
    
    ##hydrological submodel states
    state.init["SWC_nr"] <- site[["SWC_nr"]]
    state.init["SWC_rz"] <- site[["SWC_rz"]]
    state.init["MaxASW_state"] <- site[["MaxASW"]]
    state.init["E_S"] <- 0
    state.init["volSWC_rz"] <- 0
    state.init["soilRad"] <- 100
    state.init["totalRad"] <- 100
    state.init["potentialEvap"] <- 0.01
    
    
    if (state.init[["t"]] == 0) {
        state.init <- CreateNewPlantation(state = state.init, 
                                          parms = parms)
        state.init <- UpdateSoil(state = state.init, parms = parms,
                                 site = site, general.info = general.info,
                                 weather = weather)
    }
    return(state.init)
}
