
## Extract simulated data for use in likelihood function
sampleOutput<-function(df,sY,eY){
  m<-c(filter(df,Year>=sY&Year<=eY)$GPP,
       filter(df,Year>=sY&Year<=eY)$NPP,
       filter(df,Year>=sY&Year<=eY)$NEE,
       filter(df,Year>=sY&Year<=eY)$Reco,
       filter(df,Year>=sY&Year<=eY)$Rs,
       filter(df,Year>=sY&Year<=eY)$Etransp,
       filter(df,Year==2015&Month==8)$LAI,
       filter(df,Year==2018&Month==8)$LAI,
       filter(df,Year==2018&Month==8)$N,
       filter(df,Year==2018&Month==8)$dg,
      # filter(df,Year==2015&Month==7)$Wr,
      # filter(df,Year==2015&Month==7)$difRoots,
       filter(df,Year==2015&Month==7)$totC,
       filter(df,Year==2015&Month==7)$totN,
       filter(df,Year>=sY&Year<=eY)$volSWC_rz
       
  )
  return(m)
}





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
NLL <- function(p){
  sitka[.GlobalEnv$nm]<-p
  output<-do.call(fr3PGDN,sitka)
  #use sampleOutputTS if using smaller time-steps
  modelled <-sampleOutput(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
  
  #Some observed SWC values are missing, so NA, thus need to remove any NA results for sum
  NlogLik  <- sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T)
  return(NlogLik)
}



## Likelihood function
NLL_weekly <- function(p){
  sitka[.GlobalEnv$nm]<-p
  output<-do.call(fr3PGDN,sitka)
  #use sampleOutputTS if using smaller time-steps
  modelled <-sampleOutputTS(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
  
  NlogLik  <- sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T)
  return(NlogLik)
}
