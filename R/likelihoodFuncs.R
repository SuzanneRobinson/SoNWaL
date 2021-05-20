




## Extract simulated data for use in likelihood function for shorter time-steps
#daily (and therefore weekly) data is in grams per meter square per day, so need to convert 3pg output from tons per hectare
sampleOutputTS<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Etransp~ df$Month+df$Year,FUN=sum)[,3],
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

## Extract simulated data for use in likelihood function for shorter time-steps
#monthly data is in tons per hectare per month which matches 3pg output so no need to convert
sampleOutputMonth<-function(df,sY,eY){
  #convert to average grams per m2 per day depending on timestep of model
  
  if(nrow(df)<600){
  df<-df%>%
    mutate(GPP=GPP/2*100/days_in_month(Month))%>%
    mutate(NPP=NPP/2*100/days_in_month(Month))%>%
    mutate(NEE=NEE/2*100/days_in_month(Month))%>%
    mutate(Rs=Rs/2*100/days_in_month(Month))%>%
    mutate(Reco=Reco/2*100/days_in_month(Month))
  } else {
     df<-df%>%
      mutate(GPP=GPP/2*100/7)%>%
      mutate(NPP=NPP/2*100/7)%>%
      mutate(NEE=NEE/2*100/7)%>%
      mutate(Rs=Rs/2*100/7)%>%
      mutate(Reco=Reco/2*100/7)
  }
  
  df<- filter(df,Year>=sY&Year<=eY)
  
  
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$EvapTransp~ df$Month+df$Year,FUN=sum)[,3],
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



## Extract simulated data for use in likelihood function for shorter time-steps
#IF USING WEEKLY/daily, YOU USE DAILY DATA WHICH IS IN gCm-2d-1, so need to convert prafor output from tons per hectare per week (timestep)
sampleOutputWeekly<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  df$week<-c(1:53)
  m<-c(aggregate(df$GPP~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NPP~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NEE~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Reco~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Rs~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Etransp~ df$week+df$Year,FUN=sum)[,3],
       filter(df,Year==2015&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$LAI[1],
       filter(df,Year==2018&Month==8)$N[1],
       filter(df,Year==2018&Month==8)$dg[1],
       #  filter(df,Year==2015&Month==7)$Wr[1],
       #  filter(df,Year==2015&Month==7)$difRoots[1],
       filter(df,Year==2015&Month==7)$totC[1],
       filter(df,Year==2015&Month==7)$totN[1],
       aggregate(df$volSWC_rz~ df$week+df$Year,FUN=mean)[,3]
  )
  m
  return(m)
}






## Extract simulated data for use in likelihood function for shorter time-steps
sampleOutputTS_noHyd<-function(df,sY,eY){
  if(nrow(df)<600){
    df<-df%>%
      mutate(GPP=GPP*100/days_in_month(Month))%>%
      mutate(NPP=NPP*100/days_in_month(Month))%>%
      mutate(NEE=NEE*100/days_in_month(Month))%>%
      mutate(Rs=Rs*100/days_in_month(Month))%>%
      mutate(Reco=Reco*100/days_in_month(Month))
  } else {
    df<-df%>%
      mutate(GPP=GPP*100/7)%>%
      mutate(NPP=NPP*100/7)%>%
      mutate(NEE=NEE*100/7)%>%
      mutate(Rs=Rs*100/7)%>%
      mutate(Reco=Reco*100/7)
  }
  df<- filter(df,Year>=sY&Year<=eY)
  
  m<-c(aggregate(df$GPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NPP~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$NEE~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Reco~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Rs~ df$Month+df$Year,FUN=mean)[,3],
       aggregate(df$Etransp~ df$Month+df$Year,FUN=sum)[,3],
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


## Extract simulated data for use in likelihood function for shorter time-steps
sampleOutputTS_noHyd_weekly<-function(df,sY,eY){
  df<- filter(df,Year>=sY&Year<=eY)
  df$week<-c(1:53)
  
  m<-c(aggregate(df$GPP~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NPP~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$NEE~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Reco~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Rs~ df$week+df$Year,FUN=mean)[,3]*100/7,
       aggregate(df$Etransp~ df$week+df$Year,FUN=sum)[,3],
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



## Extract simulated data for use in likelihood function
#2007 is removed as a year due to some odd observed data values, basically there's no GPP or NEE etc.
sampleOutputPineWeek<-function(df,sY,eY){
  df<-df[-1,]
df$week<-c(1:52)
  
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean))*100/7,
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1997&Year<=2014&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(volSWC_rz,na.rm=TRUE))%>%
           select(mean))
    
    
  )
  return(m)
}

sampleOutputPine<-function(df,sY,eY){
  #convert to match hytialla data format

    if(nrow(df)<600){
      df<-df%>%
        mutate(GPP=GPP*100/days_in_month(Month))%>%
        mutate(NPP=NPP*100/days_in_month(Month))%>%
        mutate(NEE=NEE*100/days_in_month(Month))%>%
        mutate(Reco=Reco*100/days_in_month(Month))
    } else {
      df<-df%>%
        mutate(GPP=GPP*100/7)%>%
        mutate(NPP=NPP*100/7)%>%
        mutate(NEE=NEE*100/7)%>%
        mutate(Reco=Reco*100/7)
    }
    
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1997&Year<2014)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(volSWC_rz,na.rm=TRUE))%>%
           select(mean))
    
    
  )
  return(m)
}




## Extract simulated data for use in likelihood function
#2007 is removed as a year due to some odd observed data values, basically there's no GPP or NEE etc.
sampleOutputPineWeek_noHyd<-function(df,sY,eY){
  df<-df[-1,]
  df$week<-c(1:52)
  
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean))*100/7,
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,week)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean))*100/7,
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean))
    
    
  )
  return(m)
}

sampleOutputPine_noHyd<-function(df,sY,eY){
  #convert to match hytialla data format
  if(nrow(df)<600){
    df<-df%>%
      mutate(GPP=GPP*100/days_in_month(Month))%>%
      mutate(NPP=NPP*100/days_in_month(Month))%>%
      mutate(NEE=NEE*100/days_in_month(Month))%>%
      mutate(Reco=Reco*100/days_in_month(Month))
  } else {
    df<-df%>%
      mutate(GPP=GPP*100/7)%>%
      mutate(NPP=NPP*100/7)%>%
      mutate(NEE=NEE*100/7)%>%
      mutate(Reco=Reco*100/7)
  }
  
  m<-c(
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(GPP,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(NEE,na.rm=TRUE))%>%
           select(mean)),
    pull(filter(df,Year>=sY&Year<=eY&Year!=2007)%>%
           group_by(Year,Month)%>%
           dplyr::summarise(mean=mean(Reco,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(LAI,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(dg,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totC,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year==1996)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(totN,na.rm=TRUE))%>%
           select(mean)),
    
    pull(filter(df,Year>=1995&Year<=2011)%>%
           group_by(Year)%>%
           dplyr::summarise(mean=mean(N,na.rm=TRUE))%>%
           select(mean))
    
  )
  return(m)
}


## Likelihood function
NLL<- function(p){
  sitka[.GlobalEnv$nm]<-p
  sitka$SWpower0<-round(sitka$SWpower0)
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputMonth(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


## Likelihood function
NLL_weekly<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputWeekly(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
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


## Likelihood function
NLL_noHYDWeekly<- function(p){
  sitka[.GlobalEnv$nm]<-p
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,sitka)
      #use sampleOutputTS if using smaller time-steps
      modelled <-sampleOutputTS_noHyd_weekly(output,.GlobalEnv$startYear,.GlobalEnv$endYear)
      
      NlogLik  <-   ifelse(any(is.na(modelled)==T),-Inf,sum(dnorm(.GlobalEnv$observed,mean=modelled,sd=.GlobalEnv$dev,log=T),na.rm = T))
      
      
    },
    error=function(cond) {
      return(-Inf)
    })
  return(NlogLik)
}


## Likelihood function
pineLL <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPine(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}


## Likelihood function
pineLL_weekly <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPineWeek(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}




## Likelihood function
pineLL_noHyd <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPine_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}


## Likelihood function
pineLL_weekly_noHyd <- function(p){
  pine[.GlobalEnv$nm]<-p
  
  
  NlogLik <- tryCatch(
    {
      output<-do.call(fr3PGDN,pine)
      modelled <-suppressWarnings(suppressMessages(sampleOutputPineWeek_noHyd(output,.GlobalEnv$startYear,.GlobalEnv$endYear)))
      ifelse(any(is.na(modelled)==TRUE),-Inf,sum(dnorm(x=.GlobalEnv$observedPine,sd =.GlobalEnv$devPine, mean=modelled,log=T),na.rm = T))
      
    },
    error=function(cond) {
      return(-Inf)
    })
  
  return(NlogLik)
}

