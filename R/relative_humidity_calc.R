relative_humidity_calc<-function(Tmean,Tref=273.16,p,q){
  TmeanK<-Tmean+273.15
  
  RH<- 0.263*p*q*(exp((17.67*(TmeanK-Tref))/(TmeanK-29.65)))^-1
  RH[RH>100]<-100
  RH[RH<0]<-0
  
  return(RH)
}