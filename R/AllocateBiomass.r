AllocateBiomass <-
function (state, site, parms, weather) #requires weather too for current month, and leaffall and leafgrow from parms, and WlDormant from state
  {
    current.month <- weather[["Month"]]
    leaf.grow <- parms[["leaf.grow"]]
    leaf.fall <- parms[["leaf.fall"]]
    #WlDormant <- state[["WlDormant"]]
    m0 <- parms[["m0"]]
    pRx <- parms[["pRx"]]
    pRn <- parms[["pRn"]]
    fN <- state[["fN"]]
    PhysMod <- state[["PhysMod"]]
    m <- m0 + (1 - m0) * fN
    pR <- pRx * pRn/(pRn + (pRx - pRn) * m * PhysMod)
    pFS2 <- parms[["pFS2"]]
    pFS20 <- parms[["pFS20"]]
    pfsPower <- log(pFS20/pFS2)/log(10)
    pfsConst <- pFS2/(2^pfsPower)
    dg <- state[["dg"]]
    pFS <- pfsConst * dg^pfsPower
    pS <- (1 - pR)/(pFS + 1)
    pF <- 1 - pR - pS
    NPP <- state[["NPP"]]
    difWl <- NPP * pF
    difWr <- NPP * pR
    difWsbr <- NPP * pS
    gammaFx <- parms[["gammaFx"]]*12/parms[["timeStp"]] #adjust to get time step rates from monthly rates?
    gammaF0 <- parms[["gammaF0"]]*12/parms[["timeStp"]]
    tgammaF <- parms[["tgammaF"]]*12/parms[["timeStp"]]
    
    
    ##Problem here for shorter time steps, as litterfall is updated based on the previous month, so weekly or daily time-steps may cause issues?
    
    #1st case: we're currently in dormancy (i.e. no W increments) and we need to check if the previous month was dormant or not, to calculate any litterfall
    if (isTRUE(isDormantSeason(current.month, leaf.grow, leaf.fall))) { #if this is true, then there can be no W increments, but there might still be litterfall
      if (current.month - 1 == 0) { #then we are in January so the previous month wast December (month = 12, not month = 0)
        if(isFALSE(isDormantSeason(12, leaf.grow, leaf.fall))) { #Thus, if we are in January and December was not a dormant month, we have litterfall
          Littfall <- state[["WlDormant"]]
        } else { #otherwise if December was also a dormant month, there is no litterfall
          Littfall <- 0
        }
      } else { #If we are not in January, any previous month can be idexed by going -1 on the current month:
        if (isFALSE(isDormantSeason(current.month - 1, leaf.grow, leaf.fall))) { #e.g. if the previous month was not dormant, we have litterfall
          Littfall <- state[["WlDormant"]]
        } else { #otherwise if the previous month was also dormant, we don't have litterfall
          Littfall <- 0
        }
      }
      #Calculate end-of-month biomass: during the dormant period there are no W changes unless all the leaves fell the current month, which would change litter W
      state[["Wl"]] <- 0 #no leaves W
      difRoots <- 0
      state[["Wr"]] <- state[["Wr"]] + difRoots #no changes in root W
      state[["Wsbr"]] <- state[["Wsbr"]] ##no changes in stem W
      difLitter <- ifelse(Littfall == 0, state[["Wlitt"]] * Littfall, Littfall)
      TotalLitter <- state[["Wlitt"]] + difLitter #In a dormant month, difLitter is calculated as above and added to the existing litter W
      #browser()
    }
    
    #2nd case: The previous month was dormant but the current one is not (so we might have W increments IF we have available NPP after re-growing leaves)
    else if ((current.month - 1 == 0 && isFALSE(isDormantSeason(current.month, leaf.grow, leaf.fall)) && #so, if we're in January and NOT in dormancy season...
      isTRUE(isDormantSeason(12, leaf.grow, leaf.fall))) ||  #...but December was dormant, or...
      (isFALSE(isDormantSeason(current.month, leaf.grow, leaf.fall)) && #...or we're in any month other than January and NOT in dormancy season, and...
       isTRUE(isDormantSeason(current.month - 1, leaf.grow, leaf.fall)))) { #...and the previous month WAS a dormant month, then:
      
      Wl <- state[["WlDormant"]] #we assume all the leaf W is re-grown this month, using the value of leaf W stored from the last foliated month (see further down this script) 
      #but we need to calculate any NPP Debt that needs to be paid off in full before growth can resume:
      NPPDebt <- NPP - Wl 
      #Now we calculate litterfall and root turnover, because we're not in dormancy:
      #For litterfall:
      t <- state[["t"]]
      #For info on the litterfall eq. See A.3 in Sands and Landsberg (2002). Parameterisation of 3-PG forplantation grown Eucalyptus globules.
      Littfall <- gammaFx * gammaF0/(gammaF0 + (gammaFx - gammaF0) * 
                                       exp(-parms[["timeStp"]] * log(1 + gammaFx/gammaF0) * t/tgammaF))
      difLitter <- Littfall * Wl
      #For root turnover:
      Wr <- state[["Wr"]]
      Rttover <- parms[["Rttover"]]
      difRoots <- Rttover * Wr
      #Now we can recalculate biomass increments for the first month after dormancy:
      if (NPPDebt > 0) { #there can be some growth
        difWl <- NPPDebt * pF
        difWr <- NPPDebt * pR
        difWsbr <- NPPDebt * pS
      } else { #still no increments
        difWl <- 0
        difWr <- 0
        difWsbr <- 0
      }
      #Now we calculate end-of-month W:
      state[["Wl"]] <- Wl
      state[["Wr"]] <- Wr + difWr - difRoots
      state[["Wsbr"]] <- state[["Wsbr"]] + difWsbr
      TotalLitter <- state[["Wlitt"]] + difLitter
      #browser()
    }
    
    #3rd case: If the species is not dormant the current month, nor was it the previous month, but it is a deciduous species
    else if (leaf.grow != 0) {
      #We calculate litterfall and root turnover, because we're not in dormancy:
      #For litterfall:
      Wl <- state[["Wl"]]
      t <- state[["t"]]
      Littfall <- gammaFx * gammaF0/(gammaF0 + (gammaFx - gammaF0) * 
                                       exp(-parms[["timeStp"]] * log(1 + gammaFx/gammaF0) * t/tgammaF))
      difLitter <- Littfall * Wl
      #For root turnover:
      Wr <- state[["Wr"]]
      Rttover <- parms[["Rttover"]]
      difRoots <- Rttover * Wr
      
      #We check if we have any NPPDebt to pay before any new growth can occur:
      NPPDebt <- NPP - Wl
      if (NPPDebt < 0) { #If we have any NPPDebt to pay off...
        #...firstly we recalculate NPPDebt with this month's NPP:
        NPPDebt <- NPPDebt + NPP
        #and we check if we have now paid off the NPPDebt to allow for some growth:
        if (NPPDebt > 0) {
          difWl <- NPPDebt * pF
          difWr <- NPPDebt * pR
          difWsbr <- NPPDebt * pS
        } else { #still no increments
          difWl <- 0
          difWr <- 0
          difWsbr <- 0
        }
      #Otherwise, if we've already paid off all our previous NPPDebt the month(s) before, growth can resume as normal:  
      } else if (NPPDebt >= 0) {
        difWl <- NPP * pF
        difWr <- NPP * pR
        difWsbr <- NPP * pS
      }
      #Now we calculate end-of-month W:
      state[["Wl"]] <- Wl + difWl - difLitter
      state[["Wr"]] <- Wr + difWr - difRoots
      state[["Wsbr"]] <- state[["Wsbr"]] + difWsbr
      TotalLitter <- state[["Wlitt"]] + difLitter
      #browser()
    }
    
    #4th case: The species is not dormant this month, nor was it the previous month, and is not deciduous:
    else {
      #We calculate litterfall and root turnover, because we're not in dormancy:
      #For litterfall:
      Wl <- state[["Wl"]]
      t <- state[["t"]]
      Littfall <- gammaFx * gammaF0/(gammaF0 + (gammaFx - gammaF0) * 
                                       exp(-parms[["timeStp"]] * log(1 + gammaFx/gammaF0) * t/tgammaF))
      difLitter <- Littfall * Wl
      state[["Wlitt"]] <- state[["Wlitt"]] + difLitter
      #For root turnover:
      Wr <- state[["Wr"]]
      Rttover <- parms[["Rttover"]]
      difRoots <- Rttover * Wr
      
      #Now we calculate end-of-month W:
      Wsbr <- state[["Wsbr"]]
      TotalLitter <- state[["Wlitt"]]
      state[["Wl"]] <- Wl + difWl - difLitter
      state[["Wr"]] <- Wr + difWr - difRoots
      state[["Wsbr"]] <- state[["Wsbr"]] + difWsbr
      TotalLitter <- state[["Wlitt"]] + difLitter
      #browser()
    }
    
    #Finally, for all cases apart from case 4 (i.e. for all cases when the species is deciduous) and case 1 (i.e. if we are in dormant season), we need to store the current month's Wl if it is the last month of the growing season:
    if (current.month + 1 == 13) { #thus, if we are in December
      if (isFALSE(isDormantSeason(current.month, leaf.grow, leaf.fall)) && isTRUE(isDormantSeason(1, leaf.grow, leaf.fall))) {
        state[["WlDormant"]] <- state[["Wl"]]
      }
    } else { #or if we are in any month other than December:
      if (isFALSE(isDormantSeason(current.month, leaf.grow, leaf.fall)) && isTRUE(isDormantSeason(current.month + 1, leaf.grow, leaf.fall))) {
        state[["WlDormant"]] <- state[["Wl"]]
      }
      #browser()
    }
    
#Now we complete the output with the state variables caclulated in the script:
    state[c("pR", "pFS", "pS", "pF", "difWl", "difWr", "difWsbr", 
            "Littfall", "difLitter", "difRoots", "TotalLitter")] <- c(pR, 
                                                                      pFS, pS, pF, difWl, difWr, difWsbr, Littfall, difLitter, 
                                                                      difRoots, TotalLitter)
    return(state)
  }
