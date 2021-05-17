PredictWeatherVariables <-
  function (weather) 
  {
    
    weather1<-weather%>%filter(is.na(RH)==T)
    weather2<-weather%>%filter(is.na(RH)==F)
    
    vpdFunc<-function(weather,RH=T){
      if(RH==F){
        VPD.Tmax <- 6.1078 * exp(17.269 * weather$Tmax/(237.3 + weather$Tmax))
        VPD.Tmin <- 6.1078 * exp(17.269 * weather$Tmin/(237.3 + weather$Tmin))
        VPD <- (VPD.Tmax - VPD.Tmin)
        weather$VPDx <- VPD.Tmax
        weather$VPDn <- VPD.Tmin
        weather$VPD <- VPD
      } 
      if(RH==T){
        VPD.Tmax <- 6.1078 * exp(17.269 * weather$Tmax/(237.3 + weather$Tmax))
        VPD.Tmin <- 6.1078 * exp(17.269 * weather$Tmin/(237.3 + weather$Tmin))
        VPD.Tmean <- 6.1078 * exp(17.269 * weather$Tmean/(237.3 + weather$Tmean))
        
        VPD <- VPD.Tmean*(1-weather$RH/100)
        weather$VPDx <- VPD.Tmax
        weather$VPDn <- VPD.Tmin
        weather$VPD <- VPD
      }
      
      return(weather)
      
    }
    
    #linear model to predict historical RH from other weather variables
   if(nrow(clm_df_full<600)) mod1<-lm(RH~SolarRad+Tmin,data=weather2)
   if(nrow(clm_df_full>600)) mod1<-lm(RH~SolarRad+Tmin+Rain,data=weather2)
    
    weather1$RH<- as.vector(predict(mod1,weather1))
    weather<-rbind(vpdFunc(weather1,RH=T),  vpdFunc(weather2,RH=T))
    return(weather)
  }


