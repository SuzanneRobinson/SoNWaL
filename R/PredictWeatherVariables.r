PredictWeatherVariables <-
function (weather) 
{
  
  weather1<-weather%>%filter(is.na(RH)==T)
  weather2<-weather%>%filter(is.na(RH)==F)
  
  vpdFunc<-function(weather,RH=T){
    if(RH==F){
    VPD.Tmax <- 0.61078 * exp(17.269 * weather$Tmax/(237.3 + weather$Tmax))
    VPD.Tmin <- 0.61078 * exp(17.269 * weather$Tmin/(237.3 + weather$Tmin))
    VPD <- (VPD.Tmax - VPD.Tmin)
    weather$VPDx <- VPD.Tmax
    weather$VPDn <- VPD.Tmin
    weather$VPD <- VPD
    } 
    if(RH==T){
      VPD.Tmax <- 0.61078 * exp(17.269 * weather$Tmax/(237.3 + weather$Tmax))
      VPD.Tmin <- 0.61078 * exp(17.269 * weather$Tmin/(237.3 + weather$Tmin))
      VPD.Tmean <- 0.61078 * exp(17.269 * weather$Tmean/(237.3 + weather$Tmean))
      
      VPD <- VPD.Tmean*(1-weather$RH/100)
      weather$VPDx <- VPD.Tmax
      weather$VPDn <- VPD.Tmin
      weather$VPD <- VPD
    }
    
    return(weather)
    
  }
  
  weather<-rbind(vpdFunc(weather1,RH=F),  vpdFunc(weather2,RH=T))
    return(weather)
}
