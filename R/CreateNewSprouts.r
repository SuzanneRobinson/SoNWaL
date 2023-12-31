#' CreateNewSprouts initiate new sprouts
#' @param state model state
#' @param parms model parameters
#' @return updated models state
CreateNewSprouts <-
function (state, parms, newsprouts) 
{
    Wl.s <- parms[["Wl.s"]]
    Wsbr.s <- parms[["Wsbr.s"]]
    N <- state[["N"]]
    Nharv <- state[["Nharv"]]
    Wl <- state[["Wl"]]
    Wsbr <- state[["Wsbr"]]
    N0 <- N
    if (newsprouts == 1) {
        N <- 4.26831 * Nharv + 285.96252 * (1/parms[["timeStp"]])
    }
    else {
        N <- N + 285.96252 * (1/parms[["timeStp"]])
    }
    Wl <- Wl + (N - N0) * Wl.s/1000
    Wsbr <- Wsbr + (N - N0) * Wsbr.s/1000
    state[["N"]] <- N
    state[["Wl"]] <- Wl
    state[["Wsbr"]] <- Wsbr
    return(state)
}
