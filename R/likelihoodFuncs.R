




## Extract simulated data for use in likelihood function for shorter time-steps
sampleOutputTS<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Rs~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=mean)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
     #  filter(df,Year==2015&Month==7)$Wr[1],
     #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1],
       aggregate(df$volSWC_rz~ df$Month+df$Year,FUN=mean)[,3]
  )
  m
  return(m)
}




## Likelihood function
NLL<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
  output<-do.call(fr3PGDN,sitka)
  #use sampleOutputTS if using smaller time-steps
  modelled <-sampleOutputTS(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
  
  NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
  
  
    },
  error=function(cond) {
    return(-Inf)
  })
  return(NlogLik)
}





## Extract simulated data for use in likelihood function for shorter time-steps
sampleOutputTS_noHyd<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Rs~ df$Month+df$Year,FUN=sum)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=mean)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       #  filter(df,Year==2015&Month==7)$Wr[1],
       #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1]  )
  m
  return(m)
}




## Likelihood function
NLL_noHYD<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputTS_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}