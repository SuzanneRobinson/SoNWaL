#' PredictVariablesInterest.3PG predict key 3pg variables and update states
#' @param state current state of model
#' @param parms model parameters
#' @param cod.pred whether to use allometric predictions of dbh etc.
#' @return updated model state and variables
#' @export
PredictVariablesInterest.3PG <-
function (state, parms, cod.pred, month) #leafgrow and leaffall from parms, current.month from weather.i (or from initial.month if the function is run for initialisation)
{
    leaf.grow <- parms[["leaf.grow"]]
    leaf.fall <- parms[["leaf.fall"]]
    SLA0 <- parms[["SLA0"]]
    SLA1 <- parms[["SLA1"]]
    tSLA <- parms[["tSLA"]]
    fracBB0 <- parms[["fracBB0"]]
    fracBB1 <- parms[["fracBB1"]]
    tBB <- parms[["tBB"]]
    t <- state[["t"]]
    N <- state[["N"]]
    Wl <- state[["Wl"]]
    #WlDormant <- state[["WlDormant"]]
    Wr <- state[["Wr"]]
    Wsbr <- state[["Wsbr"]]
    #Wl <- state[["Wl"]]
 
    
    SLA <- SLA1 + (SLA0 - SLA1) * exp(-log(2) * (t/tSLA)^2)
    LAI <- Wl * SLA * 0.1 

    fracBB <- fracBB1 + (fracBB0 - fracBB1) * exp(-log(2) * (t/tBB))
    Ww <- Wsbr * (1 - fracBB)
    Wb <- (Wsbr - Ww) * 0.5675
    Wbr <- (Wsbr - Ww) - Wb
    Wa <- Wl + Wsbr
    W <- Wl + Wr + Wsbr
    hdom <- EstimateH.3PG(N = N, Wa = Wa)
    wsbrg <- Wsbr * 1000/N
    if (cod.pred == "3PG") {
        rhoMax <- parms[["rhoMax"]]
        rhoMin <- parms[["rhoMin"]]
        tRho <- parms[["tRho"]]
        WoodDensity <- rhoMax + (rhoMin - rhoMax) * exp(-log(2) * 
            t/tRho)
        Vu <- Wsbr * (1 - fracBB)/WoodDensity
        dg <- (wsbrg/parms[["aS"]])^(1/parms[["nS"]])
        G <- (dg/200)^2 * pi * N
    }
    else if (cod.pred == "Allometric") {
        G <- EstimateG.3PG(N = N, Wa = Wa)
        Vub <- EstimateV.3PG(N = N, Ww = Ww)
        dg <- ObtainDg(N = N, G = G)
    }
    state[c("LAI", "Ww", "Wb", "Wbr", "Wa", "W", "hdom", "wsbrg", #return WlDormant as well
        "G", "Vu", "dg")] <- c(LAI, Ww, Wb, Wbr, Wa, W, hdom, 
        wsbrg, G, Vu, dg)
    return(state)
}
