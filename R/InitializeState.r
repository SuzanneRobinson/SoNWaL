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
                     "Nleach","SWC_nr0","MaxASW_state","E_S")
    .length.state.vec <- length(names.state)
    state.vector <- rep(NA, .length.state.vec)
    names(state.vector) <- names.state
    state.init <- state.vector
    state.init[names(stand)] <- stand
    state.init["ASW"] <- site[["ASW"]]
    state.init["SWC_nr0"] <- site[["SWC_nr0"]]
    state.init["MaxASW_state"] <- site[["MaxASW"]]
    state.init["E_S"] <- site[["E_S"]]
    
    if (state.init[["t"]] == 0) {
        state.init <- CreateNewPlantation(state = state.init, 
                                          parms = parms)
        state.init <- UpdateSoil(state = state.init, parms = parms,
                                 site = site, general.info = general.info,
                                 weather = weather)
    }
    return(state.init)
}
