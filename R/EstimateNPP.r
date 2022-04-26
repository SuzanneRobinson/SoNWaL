#' EstimateNPP
#' @param state current model state
#' @param parms model parameters
#' @return updated model state
#' @export
EstimateNPP <-
function (state, parms) 
{
    alpha <- parms[["alpha"]]
    gDM_mol <- parms[["gDM_mol"]]
    Y <- parms[["Y"]]
    phi.pa <- state[["phi.pa"]]
    fT <- state[["fT"]]
    fF <- state[["fF"]]
    fCalpha <- state[["fCalpha"]]
    fN <- state[["fN"]]
    PhysMod <- state[["PhysMod"]]
    alphaC <- alpha * fT * fF * fCalpha * fN * PhysMod
    GPP <- alphaC * gDM_mol * phi.pa/100
    NPP <- Y * GPP
    state[c("alphaC", "GPP", "NPP")] <- c(alphaC, GPP, NPP)
    return(state)
}
