

#' NPPfunc
#' @description takes hazard years, NPP_values and grid square data to return good and bad hazard years
#' @param hzYrs list of hazard years
#' @param NPP_value 
#' @param hazPeriod 
#' @grid_id grid identification 
#' @param x x value
#' @param y y value
#' @param yc_value yield class
#' @export
NPPfunc<-function(hzYrs,NPP_value,hazPeriod, grid_id, x, y, yc_value){
  
  hzYY<-as.numeric(unlist(hzYrs))
  if(hazPeriod %in% hzYY){
    badNPP<-NPP_value
    goodNPP<-NA
  } else {
    goodNPP<-NPP_value
    badNPP<-NA
  }
  return(cbind(grid_id = grid_id,badNPP=badNPP, goodNPP =goodNPP, x=x, y=y, 
               NPP = NPP_value, yc_value=yc_value))
}



#' uqFunc
#' @description takes good and bad years and grid ID, calculates risk uncertainty etc.
#' @export
uqFunc<-function(badNPP,goodNPP,numDraws){
  badNPP<-as.numeric(badNPP)
  goodNPP<-as.numeric(goodNPP)
  
  nH<-round(length(na.omit(badNPP))/numDraws) #number hazard years needs to be yearly not weekly
  n <- 30#round(length(badNPP)) # total number of years should be 30?
  
  z<-as.numeric(na.omit(c(goodNPP,badNPP)))
  
  Ez <- mean(z)
  Ez_H <- mean( badNPP, na.rm=T )
  Ez_notH <- mean(goodNPP, na.rm=T )
  
  z_H<-na.omit(badNPP)
  z_nH<-na.omit(goodNPP)
  
  pH <- nH / n
  V<- Ez_notH - Ez_H
  R<- Ez_notH - Ez
  
  a <- 1 + nH
  b <- 1 + n - nH
  
  s_Ez <- sqrt( var(z ) / n )
  s_Ez_H <- sqrt( var(z_H) / nH )
  s_Ez_notH <- sqrt( var(z_nH) / (n-nH) )
  s_pH <- sqrt( a*b/(a+b+1) ) / (a+b)
  s_V <- sqrt( s_Ez_notH^2 + s_Ez_H^2 )
  s_R <- sqrt( s_Ez_notH^2 + s_Ez^2 - 2 * s_Ez_notH^2 * (1-nH/n) )
  s_R<-if(is.na(s_R)==T) 0 else s_R
  return(data.frame(pH = pH, R = R, V = V, s_R = s_R, s_V = s_V, s_pH = s_pH))
}
