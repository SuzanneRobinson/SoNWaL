PredictWeatherVariables <-
  function (weather,VPDconv) 
  {
    
    weather1<-weather%>%filter(is.na(RH)==T)
    weather2<-weather%>%filter(is.na(RH)==F)
    
    vpdFunc<-function(weather){

        VPD.Tmax <- 0.61078 * exp(17.269 * weather$Tmax/(237.3 + weather$Tmax))
        VPD.Tmin <- 0.61078 * exp(17.269 * weather$Tmin/(237.3 + weather$Tmin))
        VPDSat <- (VPD.Tmax - VPD.Tmin)
  
        VPD.Tmean <- 0.61078 * exp(17.269 * weather$Tmean/(237.3 + weather$Tmean))
        
        VPD <- VPD.Tmean*(1-weather$RH/100)
       # VPD<-VPDSat-VP
        weather$VPD <- VPD*0.1 #convert from Pa to mbar
      
      return(weather)
      
    }
    
    #linear model to predict historical RH from other weather variables
  mod1<-lm(RH~SolarRad+log(1+Rain)+Tmean+Tmax+Tmin,data=weather2)

      weather1$RH<- as.vector(predict(mod1,weather1))
    weather<-rbind(vpdFunc(weather1),  weather2)
    return(weather)
  }


