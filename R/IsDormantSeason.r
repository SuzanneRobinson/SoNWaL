isDormantSeason <-
function(current.month, leafgrow, leaffall) {
  is.dormant <- as.logical()
  if (any(c(leafgrow, leaffall) == 0)) { #if it's a non-deciduous species
    is.dormant <- FALSE
  } else if (leafgrow > leaffall) { #else, if it's a deciduous species, and we're in the southern emishpere
    if (current.month >= leaffall && current.month < leafgrow) { #if the current month falls in [leaffall, leafgrow],
      is.dormant <- TRUE #then we are in a dormancy month
    } else is.dormant <- FALSE
  } else if (leafgrow < leaffall) { #else, if it's a deciduous species, and we're in the northern emishpere
    if (current.month < leafgrow || current.month >= leaffall) { #similarly, if the current month falls in [leafgrow, leaffall],
      is.dormant <- TRUE #then we are in a dormancy month
    } else is.dormant <- FALSE
  } 
  return(is.dormant)
}
