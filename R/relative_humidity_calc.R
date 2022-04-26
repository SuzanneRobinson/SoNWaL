
# calc relative humidity
#' @param Tmean mean temp K
#' @param Tref reference temperature (often 273.16)
#' @param p  pressure
#' @param spec_hum specific humidity
#' @return relative humidity
relative_humidity_calc<-function(Tmean,Tref=273.16,pp,spec_hum){
  TmeanK<-Tmean+273.15
  
  RH<- 0.263*pp*spec_hum*(exp((17.67*(TmeanK-Tref))/(TmeanK-29.65)))^-1
  RH[RH>100]<-100
  RH[RH<0]<-0
  
  return(RH)
}