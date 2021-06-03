PredictWeatherVariables <-
  function (weather,VPDconv) 
  {
    
    weather1<-weather%>%filter(is.na(RH)==T)
    weather2<-weather%>%filter(is.na(RH)==F)
    
    vpdFunc<-function(weather){

       # VPD.Tmax <- 0.61078 * exp(17.269 * weather$Tmax/(237.3 + weather$Tmax))
     #   VPD.Tmin <- 0.61078 * exp(17.269 * weather$Tmin/(237.3 + weather$Tmin))
      #  VPDSat <- (VPD.Tmax - VPD.Tmin)/2
  
        VPDact <-((( (0.61078 * exp(17.269 * weather$Tmean/(237.3 + weather$Tmean)))*(1-weather$RH/100)))*0.01)

      #  Tmean = 25
      #  RH =80
      # #  VPD <- 610.7*10^(7.5*Tmean/(237.3+Tmean)) #*((100-RH)/100)*0.01
     #    VPDsat <- 0.61078 * exp(17.269 * Tmean/(237.3 + Tmean))
      #   VPDair <- 0.61078 * exp(17.269 * Tmean/(237.3 + Tmean))*(1-RH/100)
      #   VPDsat-VPDair
         
        # weather$VPD.Tmax<-VPD.Tmax
       # weather$VPD.Tmin<-VPD.Tmin
       # VPD<-VPDSat-VP
        weather$VPD <-VPDact
      
      return(weather)
      
    }
    
    #linear model to predict historical RH from other weather variables
   # cut(hr,breaks = 2)
  mod1<-lm(RH~SolarRad+log(1+Rain)+Tmean+Tmax+Tmin,data=weather2)

      weather1$RH<- as.vector(predict(mod1,weather1))
    weather<-rbind(vpdFunc(weather1),  vpdFunc(weather2))
    weather$VPD[weather$VPD<0]<-0
    return(weather)
  }


