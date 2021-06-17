RunModel <-
function (stand.init, weather, site, parms, general.info = parms.general, 
    presc = presc, cod.pred = "3PG", cod.clim = "Month") 
{


        state.init <- InitializeState(stand = stand.init, site = site, 
        parms = parms, general.info = general.info, weather = weather)
    proj.results <- list()
    state <- PredictVariablesInterest.3PG(state = state.init, 
        parms = parms, cod.pred = cod.pred, month = state.init[["initial.month"]])
    t.proj <- 0
    proj.results[[1]] <- c(t.proj = t.proj, state)
    if (!missing(presc)) {
        fma.c <- presc[presc$cycle == state[["cycle"]], ] # which cycle in presc to be using depending on current cycle state
        t.nsprouts <- unique(fma.c[, "t.nsprouts"])
        fst.row.fma.app <- min(c(which(fma.c$t > state[["t"]]), 
            nrow(fma.c)))
        fma.app <- fma.c[fst.row.fma.app:nrow(fma.c), ]
        fma.app$t[which(fma.app$t < state[["t"]])] <- state[["t"]]
    }
    if (cod.clim == "Average") {
        N = stand.init[["nyears"]] * 12
    }else if (cod.clim == "Month") {
        N = nrow(weather)
    }
    j = 1
    weather.i <- as.numeric()
    for (i in 1:N) {
        if (cod.clim == "Average") {
            if (j > 12) 
                j = 1
            weather.i <- weather[j, ]
        }else if (cod.clim == "Month") {
            weather.i <- weather[i, ]
        }
      
      timeStp<-parms[["timeStp"]]
      
        state.apar <- EstimateAPAR(state = state, weather = weather.i, 
            parms = parms, general.info = general.info)
        state.mods <- CalculateModifiers(state = state.apar, 
            weather = weather.i, site = site, parms = parms, 
            general.info = general.info)
        state.npp <- EstimateNPP(state = state.mods, parms = parms)

     #  if(parms[["timeStp"]]==52) weather.i2<-clm_df_daily%>%filter(Year==weather.i$Year & week==weather.i$week)
     #  if(parms[["timeStp"]]==12) weather.i2<-clm_df_daily%>%filter(Year==weather.i$Year & Month==weather.i$Month)
     #  if(parms[["timeStp"]]==365)weather.i2=weather.i
     #   
     #   weather.i2 <- split(weather.i2, seq(nrow(weather.i2)))
     #     
     #   state.asw<- map(weather.i2, ~UpdateASW(.x,state = state.npp, 
     #                                             site = site, parms = parms, general.info = general.info))[[length(weather.i2)]]

        state.asw<-  UpdateASW(weather=weather.i,state = state.npp, 
                               site = site, parms = parms, general.info = general.info)
        

      #  wMnth= ifelse(i!=N,weather$Month[i+1],1)
      #  if(weather$Month[i]!= wMnth || i==1){
        state.walloc <- AllocateBiomass(state = state.asw, site = site, #also requires weather.i, to identify current.month
            parms = parms, weather = weather.i)
        state.mort <- EstimateMortality(state = state.walloc, 
            parms = parms)
        state.mort[["t"]] <- state.mort[["t"]] + 1/timeStp
        
        ############ ADD THE NEW SOIL SUBMODULE HERE #################
        state.soil <- UpdateSoil(state = state.mort, parms = parms,
                                 site = site, general.info = general.info,
                                 weather = weather)
        state.end <- PredictVariablesInterest.3PG(state = state.soil, 
            parms = parms, cod.pred = cod.pred, month = weather.i[["Month"]])
      #  }

        if (!missing(presc) && nrow(fma.app) > 0 && ((abs(fma.app$t[1] - 
            state.end[["t"]]) < 1/(timeStp*2)) | (fma.app$t[1] < state.end[["t"]]))) {
            state.end <- DoThinning(state = state.end, parms = parms, 
                fma = fma.app, presc = presc)
            fma.app <- fma.app[-1, ]
            state.end <- PredictVariablesInterest.3PG(state = state.end, 
                parms = parms, cod.pred = cod.pred, month = weather.i[["Month"]])
        }
        if (!missing(presc) && state.end[["N"]] == 0) {
            if (state.end[["cycle"]] > max(presc$cycle)) {
                t.proj <- t.proj + 1/timeStp
                proj.results[[i + 1]] <- c(t.proj = t.proj, state.end)
                break
            }
            fma.c <- presc[which(presc$cycle == state.end[["cycle"]]), 
                ]
            t.nsprouts <- unique(fma.c[, "t.nsprouts"])
            fma.app <- fma.c
            it.newsprouts <- 1
            if (state.end[["rotation"]] == 1) {
                state.end <- CreateNewPlantation(state = state.end, 
                  parms = parms, fma = fma.app)
                state.end <- PredictVariablesInterest.3PG(state = state.end, 
                  parms = parms, cod.pred = cod.pred, month = weather.i[["Month"]])
            }
        }
        if (!missing(presc) && state.end[["rotation"]] > 1 && 
            state.end[["t"]] < t.nsprouts) {
            state.end <- CreateNewSprouts(state = state.end, 
                parms = parms, newsprouts = it.newsprouts)
            state.end <- PredictVariablesInterest.3PG(state = state.end, 
                parms = parms, cod.pred = cod.pred, month = weather.i[["Month"]])
            it.newsprouts <- it.newsprouts + 1
        }
        if (!missing(presc) && state.end[["rotation"]] > 1 && 
            state.end[["t"]] > (t.nsprouts - 1/24) && state.end[["rm.sprouts"]] == 
            0) {
            state.end <- RemoveSprouts(state = state.end, parms = parms, 
                fma = fma.app)
        }
        t.proj <- t.proj + 1/timeStp
        proj.results[[i + 1]] <- c(t.proj = t.proj, state.end)
        state <- state.end
        j = j + 1
    }
    proj.df <- data.frame(do.call(rbind, proj.results))
    return(proj.df)
}
