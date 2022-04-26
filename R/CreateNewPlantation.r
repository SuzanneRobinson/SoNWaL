#' CreateNewPlantation initiate plantation
#' @param state model state
#' @param parms model parameters
#' @return updated models state
CreateNewPlantation <-
function (state, parms, fma) 
{
    Wl.s <- parms[["Wl.s"]]
    Wr.s <- parms[["Wr.s"]]
    Wsbr.s <- parms[["Wsbr.s"]]
    if (missing(fma)) {
        N <- state[["N"]]
    }
    else {
        N <- unique(fma$Npl)
    }
    state[["N"]] <- N
    state[["Wl"]] <- N * Wl.s/1000
    state[["WlDormant"]] <- state[["Wl"]] #added a starting point for WlDormant
    state[["Wr"]] <- N * Wr.s/1000
    state[["Wsbr"]] <- N * Wsbr.s/1000
    return(state)
}
